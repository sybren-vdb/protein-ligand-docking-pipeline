#!usr/bin/env python3
import argparse, os, random, sys, tempfile, time
from pathlib import Path

import requests
from datetime import datetime, date
import json
import gemmi

WATER_NAMES = {"HOH", "DOD", "WAT", "H2O"}

RCSB_ENTRY_IDS = "https://data.rcsb.org/rest/v1/holdings/current/entry_ids"
FASTA_ENTRY_URL = "https://www.rcsb.org/fasta/entry/{pdbid}/download"
MMCIF_URL = "https://files.rcsb.org/download/{pdbid}.cif"
MMCIF_ASMB_URL = "https://files.rcsb.org/download/{pdbid}-assembly{n}.cif"
SEARCH_URL = "https://search.rcsb.org/rcsbsearch/v2/query"
CORE_ENTRY_URL = "https://data.rcsb.org/rest/v1/core/entry/{pdbid}"

"""def get_all_current_pdb_ids(session: requests.Session) -> list[str]:
	Fetch all current PDB IDs (uppercase strings).
	for i in range(3):
		r = session.get(RCSB_ENTRY_IDS, timeout=30, headers={
			"User-Agent": "boltz2-retraining)
	r.raise_for_status()
	data = r.json()
	# Endpoint returns a list
	if not isinstance(data, list) or not all(isinstance(x, str) for x in data):
		raise RuntimeError("Unexcepted /entry_ids payload shape")
	return [x.upper() for x in data]"""

def get_ids_by_date(session, after: date | None = None, page: int = 10000) -> list[str]:
	ids: list[str] = []
	start = 0
	while True:
		nodes = []
		if after:
			nodes.append({
				"type": "terminal",
				"service": "text",
				"parameters": {
					"attribute": "rcsb_accession_info.initial_release_date",
					"operator": "greater_or_equal",
					"value": after.strftime("%Y-%m-%d"),
					},
				})
		if not nodes:
			query = {
			"type": "terminal", "service": "full_text", "parameters": {"value": "*"}}
		elif len(nodes) == 1:
			query = nodes[0]
		else:
			query = {"type": "group", "logical_operator": "and", "nodes": nodes}
		payload = {
			"query": query,
			"return_type": "entry",
			"request_options": {
				"paginate": {"start": start, "rows": page},
				"results_content_type": ["experimental"],
			},
		}
		r = session.post(SEARCH_URL, json=payload, timeout=45)
		if r.status_code == 400:
			raise RuntimeError(f"RCSB 400: {r.text}")
		r.raise_for_status()

		batch = r.json().get("result_set", [])
		if not batch:
			break
		ids.extend(row["identifier"].upper() for row in batch)

		if len(batch) < page:
			break
		start += page

	return ids

def download_text(session: requests.Session, url: str) -> str:
	r = session.get(url, timeout=45)
	r.raise_for_status()
	return r.text

def download_bytes(session: requests.Session, url: str) -> bytes:
	r = session.get(url, timeout=60)
	r.raise_for_status()
	return r.content

def write_fasta_for_entry(session: requests.Session, pdbid: str, out_fasta: Path):
	#Per entry FASTA
	fasta = download_text(session, FASTA_ENTRY_URL.format(pdbid=pdbid))
	#Append to fasta
	with out_fasta.open("a", encoding="utf-8") as fh:
		#Normalize header
		for line in fasta.splitlines():
			if line.startswith(">"):
				if pdbid not in line.upper():
					fh.write(f">{pdbid} " + line[1:] + "\n")
				else:
					fh.write(line + "\n")
			else:
				fh.write(line + "\n")

def strip_waters(mmcif_bytes: bytes) -> bytes:
	"""Return mmCIF bytes with all solvent waters removed."""
	with tempfile.TemporaryDirectory() as td:
		in_p = os.path.join(td, "in.cif")

		with open(in_p, "wb") as f:
			f.write(mmcif_bytes)
	
		st = gemmi.read_structure(in_p)
		for model in st:
			for chain in model:
				to_delete = []
				for i in range(len(chain)):
					res = chain[i]
					if (hasattr(res, "is_water") and res.is_water()) or (res.name.upper() in WATER_NAMES):
						to_delete.append(i)
				for i in reversed(to_delete):
					del chain[i]
		for model in st:
			for i in range(len(model)):
				if len(model[i]) == 0:
					del model[i]
		doc = st.make_mmcif_document()
		return doc.as_string().encode("utf-8")

def parse_date(s: str) -> date:
	try:
		return datetime.strptime(s, "%Y-%m-%d").date()
	except Exception as e:
		raise argparse.ArgumentTypeError(f"Invalid year and month {s!r}, expected YYYY-MM-DD") from e

def get_release_date(session: requests.Session, pdbid: str) -> date | None:
	"""
	Return the entry's initial public release date (preferred), else deposite date, else None.
	"""
	r = session.get(CORE_ENTRY_URL.format(pdbid=pdbid), timeout=30)
	r.raise_for_status()
	j = r.json()
	info = j.get("rcsb_accession_info") or {}
	for key in ("initial_release_date", "deposit_date"):
		v = info.get(key)
		if v:
			try:
				return datetime.fromisoformat(v[:10]).date()
			except Exception:
				pass
	return None

def main():
	p = argparse.ArgumentParser(
		description="Sample n random PDBs, append FASTA to one file, and write water-stripped mmCIFs."
	)
	p.add_argument("-n", type=int, required=True, help="Number of random PDBs to sample")
	p.add_argument("--out-fasta", required=True, help="Path to output FASTA file (single file)")
	p.add_argument("--out-dir", required=True, help="Directory to write mmCIFs")
	p.add_argument("--seed", type=int, default=None, help="RNG seed for reproducibility")
	p.add_argument("--biounit", type=int, default=None, help="If set, download assemblyN.cif instead of asym unit")
	p.add_argument("--no-strip", action="store_true", help="Do NOT strip waters (write raw mmCIF")
	p.add_argument("--after-date", type=parse_date, default=None, help="Only include entries with initial release >= this YYYY-MM-DD")
	p.add_argument("--skip-existing", action="store_true", help="Skip IDs whose CIF already exists in out_dir")
	args = p.parse_args()
	
	rng = random.Random(args.seed)
	out_fasta = Path(args.out_fasta)
	out_dir = Path(args.out_dir)
	out_dir.mkdir(parents=True, exist_ok=True)
		
	if args.no_strip:
		print("[INFO] --no-strip set: mmCIFs will be saved as-downloaded.")

	s = requests.Session()
	
	if args.after_date:
		all_ids = get_ids_by_date(s, after=args.after_date)
		if not all_ids:
			print("[WARN] No entries match the date window.", file=sys.stderr)
	else:	
		all_ids = get_all_current_pdb_ids(s)
	rng.shuffle(all_ids)

	picked: list[str] = []
	examined = 0
	for pdbid in all_ids:
		if args.after_date:
			d = None
			try:
				d = get_release_date(s, pdbid)
				print("[DATE]", pdbid, d if d else "None")
			except requests.HTTPError as e:
				print(f"[HTTP {e.response.status_code}] {pdbid}: {e.response.text[:200]}")
				continue
			except Exception as e:
				print(f"[ERR] {pdbid}: {e}")
				continue
			if not d or d < args.after_date:
				continue
		print(pdbid)
		picked.append(pdbid)
		if len(picked) >= args.n:
			break

	if len(picked) < args.n:
		print(f"[WARN] Only found {len(picked)} entries meeting criteria (requested {args.n}).", file=sys.stderr)

	for k, pdbid in enumerate(picked, 1):
		t0 = time.time()
		pdbid_up = pdbid.upper()

		try:
			write_fasta_for_entry(s, pdbid_up, out_fasta)
		except Exception as e:
			print(f"[WARN] FASTA failed for {pdbid_up}: {e}", file.sys.stderr)

		try:
			if args.biounit:
				url = MMCIF_ASMB_URL.format(pdbid=pdbid_up, n=args.biounit)
			else:
				url = MMCIF_URL.format(pdbid=pdbid_up)
			mmcif_bytes = download_bytes(s, url)
		except Exception as e:
			print(f"[WARN] FASTA failed for {pdbid_up}: {e}", file=sys.stderr)

		try:
			if args.biounit:
				url = MMCIF_ASMB_URL.format(pdbid=pdbid_up, n=args.biounit)
			else:
				url = MMCIF_URL.format(pdbid=pdbid_up)
			mmcif_bytes = download_bytes(s, url)
		except Exception as e:
			print(f"[WARN] mmCIF download failed for {pdbid_up}: {e}", file=sys.stderr)
			continue

		try:
			if not args.no_strip:
				mmcif_bytes = strip_waters(mmcif_bytes)
		except Exception as e:
			print(f"[WARN] Water stripping failed for {pdbid_up} {e}", file=sys.stderr)

		out_path = out_dir / f"{pdbid_up}.cif"
		if args.skip_existing and out_path.exists():
			print(f"[SKIP] {pdbid_up} (exists)")
		else:
			with open(out_path, "wb") as f:
				f.write(mmcif_bytes)
			dt = time.time() - t0
			print(f"[OK] {k:4d}/{args.n} {pdbid_up} -> {out_path} ({dt:.1f}s)")

	print(f"[DONE] Wrote mmCIFs to {out_dir} and appended FASTA to {out_fasta}")

if __name__ =="__main__":
	main()
		
