#!/usr/bin/env python3
import argparse
import os
import re
import subprocess
import sys
import tempfile
from pathlib import Path
import requests
import argparse
import gzip, io

import rdkit
from rdkit import Chem

WATER_NAMES = {"HOH", "WAT"}

def load_pdb_text(pdb_source):
	"""Accepts raw pdb text, a string path, or a Path. Transparently handles .pdb.gz."""
	if isinstance(pdb_source, str) and (
		"\n" in pdb_source or pdb_source.startswith(("ATOM ", "HETATM", "HEADER"))
			):
		return pdb_source
	p = Path(pdb_source)
	data = p.read_bytes()
	if data[:2] == b"\x1f\x8b":
		with gzip.GzipFile(fileobj=io.BytesIO(data)) as gz:
			return gz.read().decode("utf-8", errors="replace")
	return data.decode("utf-8", errors="replace")

def run(cmd, cwd=None):
	print("[RUN]", " ".join(cmd))
	res = subprocess.run(cmd, cwd=cwd, stdout=subprocess.PIPE, 
				stderr=subprocess.STDOUT, text=True)
	if res.returncode != 0:
		print(res.stdout, file=sys.stderr)
		raise RuntimeError(f"Command failed: {' '.join(cmd)}")
	return res.stdout

def ensure_obabel():
	try:
		out = run(["obabel", "-V"])
		if "Open Babel" not in out:
			raise RuntimeError
	except Exception:
		raise RuntimeError("Open Babel (obabel) not found on Path.")

def download_pdb(pdb_id: str, out_pdb: Path):
	url = f"https://files.rcsb.org/download/{pdb_id.upper()}.pdb"
	r = requests.get(url, timeout=60)
	r.raise_for_status()
	out_pdb.write_text(r.text)
	return out_pdb

# ----Minimal PDB parsing helpers (ATOM/HETATM + CONECT) ----

def parse_atom_line(line):
	"""Return dict with PDB atom fields (1-based columns)."""
	#PDB fixed width fields
	rec = line[0:6].strip()
	if rec not in ("ATOM", "HETATM"):
		return None
	
	s = line.rstrip("\n")
	if len(s) < 80:
		s = s.ljust(80)
		
	return {
		"record": rec,
		"serial": int(s[6:11]),
		"atom": s[12:16],
		"altloc": s[16],
		"resname": s[17:20].strip(),
		"chain": s[21].strip(),
		"resseq": s[22:26],
		"icode": s[26].strip() or "",
		"x": float(line[30:38]),
		"y": float(line[38:46]),
		"z": float(line[46:54]),
		"element": line[76:78].strip() if len(line) >= 78 else "",
		"line": s
	}

def parse_conect_line(line):
	if not line.startswith("CONECT"):
		return None
	#Serial numbers are right-justified in 5 char fields after "CONECT"
	fields = [line[6+i:11+i] for i in range(0, 25, 5)]
	ints = []
	for f in fields:
		f = f.strip()
		if f:
			try:
				ints.append(int(f))
			except ValueError:
				pass
	return ints

def parse_ligand_selector(ligand_id: str):
	"""
	'LIG' or A:LIG' or 'A:LIG:123' (chain:resname:resseq)
    Returns (chain or None, resname, resseq or None)
    """
	pat_full = re.compile(r"^(?P<chain>[^:]):(?P<res>\w{1,3}):(?P<id>\d+)$")
	m = pat_full.match(ligand_id)
	if m:
		return (m.group("chain"), m.group("res").upper(), int(m.group("id")))
	pat_chain = re.compile(r"^(?P<chain>[^:]):(?P<res>\w{1,3})$")
	m = pat_chain.match(ligand_id)
	if m:
		return (m.group("chain"), m.group("res").upper(), None)
	return (None, ligand_id.upper(), None)

def split_receptor_and_pick_ligand(pdb_source, ligand_id_str: str):
	"""
	Returns:
		receptor_lines: list[str] (ATOM, non-water; no HETATM)
		ligand_lines: list[str] (HETATM for the chosen residue only)
		conect_lines_filtered: list[str] (CONECT lines involving ligand_atoms_only)
	"""
	chain_sel, resname_sel, resseq_sel = parse_ligand_selector(ligand_id_str)

	atoms = []
	conects_raw = []
	pdb_raw = load_pdb_text(pdb_source)
	for ln in pdb_raw.splitlines():
		if ln.startswith(("ATOM ", "HETATM")):
			rec = parse_atom_line(ln)
			if rec:
				atoms.append(rec)
		elif ln.startswith("CONECT"):
			cs = parse_conect_line(ln)
			if cs:
				conects_raw.append((cs, ln.rstrip("\n")))
		else:
			pass

	#build resi buckets for HETATM
	het_by_residue = {}
	receptor_lines = []
	for a in atoms:

		resname = str(a["resname"]).strip().upper()
		chain   = str(a["chain"]).strip()
		resseq  = str(a["resseq"]).strip()
		icode   = str(a.get("icode","")).strip()
		
		if a["record"] == "ATOM":
			if a["resname"] not in WATER_NAMES:
				receptor_lines.append(a["line"])
		elif a["record"] == "HETATM":
			if resname in WATER_NAMES:
				continue
			key = (chain, resname, resseq, icode)
			het_by_residue.setdefault(key, []).append(
				{"line": a["line"], "serial": int(a["serial"])}
				)
			
	#filter candidate residues
	resname_sel_norm = None if resname_sel is None else str(resname_sel).strip().upper()
	chain_sel_norm   = None if chain_sel   is None else str(chain_sel).strip()
	resseq_sel_norm  = None if resseq_sel  is None else str(resseq_sel).strip()
	candidates = []
	for key, alist in het_by_residue.items():
		chain, resname, resseq, icode = key
		if resname_sel_norm is not None and resname != resname_sel_norm:
			continue
		if chain_sel_norm is not None and chain != chain_sel_norm:
			continue
		if resseq_sel_norm is not None and resseq != resseq_sel_norm:
			continue
		candidates.append((len(alist), (chain, resname, resseq, icode), alist))

	if not candidates:
		present_names = sorted({k[1] for k in het_by_residue.keys()})
		present_by_chain = sorted({(k[0], k[1]) for k in het_by_residue.keys()})
		raise ValueError(
			f"No HETATM residue found matching '{ligand_id_str}'. "
			f"Seen HET names: {present_names[:25]}; "
			f"examples (chain,resname): {present_by_chain[:10]}"
				)
	
	#if resseq not specified, we choose the largest residue instance instead
	candidates.sort(key=lambda t: t[0], reverse=True)
	_, chosen_key, chosen_atoms = candidates[0]

	#Write ligand HETATM lines (preserve original)
	ligand_lines = [a["line"] for a in chosen_atoms]
	
	#Filter CONECT to only those involving ligand atom serials
	ligand_serials = {a["serial"] for a in chosen_atoms}
	conect_lines_filtered = []
	for ints, raw in conects_raw:
		if any(s in ligand_serials for s in ints):
			conect_lines_filtered.append(raw)

	
	return receptor_lines, ligand_lines, conect_lines_filtered, chosen_key

def write_pdb(lines, out_path: Path):
	with out_path.open("w") as fh:
		for ln in lines:
			fh.write(ln if ln.endswith("\n") else ln + "\n")
		fh.write("END\n")

def maybe_write_autobox_mol2_with_rdkit(ligand_pdb: Path, out_mol2: Path) -> bool:
	"""Try RDKit to convert ligand PDB -> MOL2 (optional). Returns True if written."""
	if not HAVE_RDKIT:
		return False
	pdb_block = ligand_pdb.read_text()
	# Keep original coordinates; avoid sanitization that may fail without full chemistry info
	mol = Chem.MolFromPDBBlock(pdb_block, sanitize=False, removeHs=False)
	if mol is None:
		return False
	# RDKit's Mol2 writer is available in recent versions
	try:
		Chem.MolToMol2File(mol, str(out_mol2))
		return True
	except Exception:
		return False

#----Prep and Docking

def prepare_receptor_pdbqt(receptor_pdb: Path, out_pdbqt: Path):
	#protonate at pH 7.4 and assign Gasteiger charges via obabel
	run([
		"obabel", "-ipdb", str(receptor_pdb),
		"-opdbqt",
		"-p", "7.4",
		"--partialcharge", "gasteiger",
		"-O", str(out_pdbqt)
		])

def run_gnina(receptor_pdbqt: Path, ligands_sdf: Path, autobox_ref: Path,
	out_poses: Path, out_log: Path,
	seed: int = 0, exhaustiveness: int = 8, autobox_add: float = 4.0):
	cmd = [
		"gnina",
		"-r", str(receptor_pdbqt),
		"-l", str(ligands_sdf),
		"--autobox_ligand", str(autobox_ref),
		"--autobox_add", str(autobox_add),
		"--exhaustiveness", str(exhaustiveness),
		"--seed", str(seed),
		"--num_modes", "1",
		"--out", str(out_poses),
		"--log", str(out_log)
	]
	with out_log.open("a") as lf:
		proc = subprocess.run(cmd, stdout=lf, stderr=subprocess.STDOUT, text=True)
		if proc.returncode != 0:
			raise RuntimeError("gnina run failed; see log.")

def main():
	ap =argparse.ArgumentParser(
		description="Dock an SDF into the pocket defined by a PDB ligand using gnina AutoBox."
	)
	ap.add_argument("--pdb-id", required=True, help="RCSB PDB ID, eg 8C1Y")
	ap.add_argument("--ligand-id", required=True, help="Ligand selector: 'LIG' or 'A:LIG' or 'A:LIG:123'")
	ap.add_argument("--sdf", required=True, type=Path, help="Input molecules (SDF)")
	ap.add_argument("--out-poses", default="poses.sdf", type=Path, help="Output docked poses SDF")
	ap.add_argument("--out-log", default="dock.log", type=Path, help="gnina run log")
	ap.add_argument("--seed", type=int, default=0)
	ap.add_argument("--autobox-add", type=float, default=4.0)
	ap.add_argument("--exhaustiveness", type=int, default=8)
	ap.add_argument("--keep-temp", action="store_true")
	ap.add_argument("--autobox-format", choices=["pdb", "mol2"], help="Use ligand PDB directly (default) or try to emit MOL2 with RDKit.")
	args = ap.parse_args()

	ensure_obabel()
	
	with tempfile.TemporaryDirectory() as tdir_:
		tdir = Path(tdir_)
		pdb_raw = tdir / f"{args.pdb_id.upper()}.pdb"
		receptor_pdb = tdir / "receptor_only.pdb"
		receptor_pdbqt = tdir / "receptor.pdbqt"
		ligand_pdb = tdir / "autobox_ligand.pdb"
		ligand_mol2 = tdir / "autobox_ligand.mol2"
		
		print(f"[+] Downloading PDB {args.pdb_id} ...")
		download_pdb(args.pdb_id, pdb_raw)

		print("[+] Splitting receptor and selecting ligand:", args.ligand_id)
		rec_lines, lig_lines, conect_lines, chosen_key = split_receptor_and_pick_ligand(pdb_raw, args.ligand_id)
		#Write files
		write_pdb(rec_lines, receptor_pdb)
		write_pdb(lig_lines, ligand_pdb)
		sel_chain, sel_res, sel_resseq, sel_icode = chosen_key
		print(f"[i] Chosen ligand instance: chain={sel_chain or ''} res={sel_res} id={sel_resseq}{sel_icode or ''}")

		print("[+] Preparing receptor PDBQT via Open Babel ...")
		prepare_receptor_pdbqt(receptor_pdb, receptor_pdbqt)

		#decide autobox ref
		autobox_ref = ligand_pdb
		if args.autobox_format == "mol2":
			if not HAVE_RDKIT:
				print("[!] RDKit not available; falling back to ligand PDB for --autobox_ligand.")
			else:
				print("[+] Writing ligand MOL2 via RDKit ...")
			ok = maybe_write_autobox_mol2_with_rdkit(ligand_pdb, ligand_mol2)
			if ok:
				autobox_ref = ligand_mol2
			else:
				print("[!] RDKit MOL2 export failed; falling back to ligand PDB.")
			
		print("[+] Running gnina with AutoBox ...")
		run_gnina(
			receptor_pdbqt=receptor_pdbqt,
			ligands_sdf=args.sdf,
			autobox_ref=autobox_ref,
			out_poses=args.out_poses,
			out_log=args.out_log,
			seed=args.seed,
			exhaustiveness=args.exhaustiveness,
			autobox_add=args.autobox_add
		)

		print(f"[✓] Done.\n  Poses: {args.out_poses}\n  Log:   {args.out_log}")
		if args.keep_temp:
			keep_dir = Path.cwd() / f"gnina_tmp_{args.pdb_id}"
			keep_dir.mkdir(parents=True, exist_ok=True)
			for p in [pdb_raw, receptor_pdb, receptor_pdbqt, ligand_pdb, ligand_mol2]:
				if p.exists():
					p.replace(keep.dir / p.name)
			print(f"[i] Kept intermediates in: {keep_dir}")

if __name__ == "__main__":
	main()
