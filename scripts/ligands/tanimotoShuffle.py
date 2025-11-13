import os
import sys
import argparse
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem
import random
from collections import defaultdict, deque
from typing import Dict, List, Tuple, Optional
import base64, zlib
from pathlib import Path

ASSIGN_TAG = "__assigned_dir__"
FP_TAG_DEFAULT = "__morgan_fp__"
SIM_TAG_DEFAULT = f"__Tanimoto_sim__"
RADIUS_DEFAULT = 2
NBITS_DEFAULT = 2048

def mol_fp_from_tag(mol: Chem.Mol, 
			fp_tag: str = FP_TAG_DEFAULT,
			radius: int = RADIUS_DEFAULT,
			nBits: int = NBITS_DEFAULT,
			write_back: bool = False):
	"""Attempt to read an RDKit bitvector from SD tag "Morgan_Fingerprint. If missing, compute."""
	if mol is None:
		return None
	if mol.HasProp(fp_tag):
		tag = mol.GetProp(fp_tag)
		try:
			return DataStructs.CreateFromBitString(tag)
		except Exception:
			try:
				raw = zlib.decompress(base64.b64decode(tag.encode()))
				return DataStructs.CreateFromBinary(raw)
			except Exception as e:
				print(f"[ERROR] {context}: {e}")
				raise

	#recompute
	bv = AllChem.GetMorganFingerprintAsBitVect(mol, radius=radius, nBits=nBits)
	if write_back:
		mol.SetProp(fp_tag, bv.ToBitString())
	return bv

def read_sdf_ids_fps(sdf_path: str,
			id_prop: Optional[str] = None,
			fp_tag: str = FP_TAG_DEFAULT,
			radius: int = RADIUS_DEFAULT,
			nBits: int = NBITS_DEFAULT,
			write_back: bool = False) -> Tuple[List[str], List[Chem.Mol], List[DataStructs.ExplicitBitVect]]:
	"""Read an SDF and return parallel lists: ids, mols, fps."""
	supp = Chem.SDMolSupplier(sdf_path, removeHs=False)
	ids, mols, fps = [], [], []
	idx = 0
	for mol in supp:
		if mol is None:
			idx += 1
			continue
		if id_prop and mol.HasProp(id_prop):
			mid = mol.GetProp(id_prop)
		elif mol.HasProp('_Name') and mol.GetProp('_Name'):
			mid = mol.GetProp('_Name')
		else:
			mid = f"mol_{idx}"
		bv = mol_fp_from_tag(mol, fp_tag=fp_tag, radius=radius, nBits=nBits, write_back=write_back)
		if bv is None:
			idx += 1
			continue
		ids.append(mid)
		mols.append(mol)
		fps.append(bv)
		idx += 1
	return ids, mols, fps

	#similarity helper
def max_sim_to_set(fp, fps_list):
	if not fps_list:
		return 0.0
	return max(DataStructs.BulkTanimotoSimilarity(fp, fps_list))

	
	#The Shuffler: uses Greedy logic with built-in failsafes

def shuffle_core(ids: List[str],
			fps: List[DataStructs.ExplicitBitVect],
			inc_fps_by_dir: Dict[str, List[DataStructs.ExplicitBitVect]],
			target_counts: Dict[str, int],
			thr: float = 0.3,
			enforce_newcomer_diversity: bool = False,
			max_restarts: int = 10,
			rng=None) -> Optional[Dict[str, str]]:
	
	dirs = list(target_counts.keys())
	d2i = {d:i for i,d in enumerate(dirs)}
	N = len(fps)
	
	#feasibility against included sets
	feasible = []
	for i in range(N):
		feas_i = []
		for d in dirs:
			feas_i.append(max_sim_to_set(fps[i], inc_fps_by_dir.get(d, [])) < thr)
		feasible.append(feas_i)

	best_sol = None
	for _ in range(max_restarts):
		order = list(range(N))
		(rng or random).shuffle(order)
		remaining = target_counts.copy()
		placed_in_dir = {d: [] for d in dirs}
		assignment_idx = {i: None for i in range(N)}
		deferred = deque()
 
		#added helper for workflows where ligand diversity is enforced
		def ok_newcomers(i, d):
			if not enforce_newcomer_diversity or not placed_in_dir[d]:
				return True
			sims = DataStructs.BulkTanimotoSimilarity(fps[i], [fps[j] for j in placed_in_dir[d]])
			return max(sims) < thr

	
		# pass greedily
		for i in order:
			cands = [di for di, d in enumerate(dirs) if feasible[i][di] and remaining[d] > 0]
			cands.sort(key=lambda di: -remaining[dirs[di]])
			placed = False
			for di in cands:
				d = dirs[di]
				if ok_newcomers(i, d):
					assignment_idx[i] = d
					placed_in_dir[d].append(i)
					remaining[d] -= 1
					placed = True
					break
			if not placed:
				deferred.append(i)

		if not deferred:
			best_sol = {ids[i]: assignment_idx[i] for i in range(N)}
			break

		#local repair via swaps
		success = True
		tries = 0
		while deferred and tries < 5 * N:
			i = deferred.popleft()
			placed = False
			cands = [di for di, d in enumerate(dirs) if feasible[i][di] and remaining[d] > 0]
			for di in cands:
				d = dirs[di]
				if ok_newcomers(i, d):
					assignment_idx[i] = d
					placed_in_dir[d].append(i)
					remaining[d] -= 1
					placed = True
					break
			if placed:
				continue
				#now we attempt a swap
			for d in dirs:
				for j in list(placed_in_dir[d]):
					if assignment_idx[j] != d:
						continue
					for d2 in dirs:
						if d2 == d or remaining[d2] == 0:
							continue
						# i -> d (replaces j), j -> d2
						di_ok = (feasible[i][d2i[d]] and 
									(not enforce_newcomer_diversity or (not placed_in_dir[d] or 
									max(DataStructs.BulkTanimotoSimilarity(fps[i], [fps[k] for k in placed_in_dir[d] if k != j])) < thr)))
						dj_ok = (feasible[j][d2i[d]] and
                                                                       (not enforce_newcomer_diversity or (not placed_in_dir[d2] or
                                                                        max(DataStructs.BulkTanimotoSimilarity(fps[j], [fps[k] for k in placed_in_dir[d2]])) < thr)))
						if di_ok and dj_ok:
							placed_in_dir[d].remove(j)
							assignment_idx[i] = d
							placed_in_dir[d].append(i)
							
							assignment_idx[j] = d2
							placed_in_dir[d2].append(j)
							remaining[d2] -= 1
							placed = True
							break
					if placed:
						break
				if placed:
					break
			if not placed:
				success = False
			tries += 1
		if success and all(remaining[d] == 0 for d in dirs):
			best_sol = {ids[i]: assignment_idx[i] for i in range(N)}
			break
	return best_sol

	#  decoy SDF writer
def write_assigned_sdfs(out_root: Path,
				assignment: Dict[str, str],
				id_to_mol: Dict[str, Chem.Mol],
				inc_fps_by_dir: Dict[str, List[DataStructs.ExplicitBitVect]],
				fp_tag: str,
				sim_tag: str,
				thr: float):
	per_dir_mols: Dict[str, List[Chem.Mol]] = defaultdict(list)
	#attach tags/bucket
	for mid, d in assignment.items():
		mol = id_to_mol.get(mid)
		if mol is None:
			continue
		#similarity vs included set for assigned directory
		fp = fp_from_tag_or_compute(mol, fp_tag=fp_tag, radius=RADIUS_DEFAULT, nBits=NBITS_DEFAULT, write_back=False)
		sim = max_sim_to_set(fp, inc_fps_by_dir.get(d, [])) if fp is not None else None
		try:
			mol.SetProp(ASSIGN_TAG, d)
		except Exception:
			pass
		if sim is not None:
			try:
				mol.SetProp(sim_tag, f"{sim:.4f}")
			except Exception as e:
				print(f"[ERROR] sim_tag={sim_tag}: {e}", file=sys.stderr)
				raise
		# ensure FP tag exists (bitstring)
		if not (mol.HasProp(fp_tag)):
			fp2 = fp or fp_from_tag_or_compute(mol, fp_tag=fp_tag, radius=RADIUS_DEFAULT, nBits=NBITS_DEFAULT, write_back=False)
			if fp2 is not None:
				try:
					mol.SetProp(fp_tag, fp2.ToBitString())
				except Exception as e:
					print(f"[ERROR] fp_tag={fp_tag}: {e}", file=sys.stderr)
					raise
		per_dir_mols[d].append(mol)

	#write SDFs into {dir}_decoy/decoys.sdf
	for d, mols in per_dir_mols.items():
		out_dir = out_root.parent / f"{out_root.name}"
		pass

#------------Main------------

def main():
	ap = argparse.ArgumentParser(description="Shuffle binders across protein directories and write decoys to <dir>_decoy/decoys.sdf")
	ap.add_argument("--binders-name", default="binders.sdf", help="File name of binders SDF inside each protein directory.")
	ap.add_argument("--included-name", default="included.sdf", help="Exact filename searched at the top of each folder")
	ap.add_argument("--thr", type=float, default=0.3, help="Tanimoto threshold; require < thr vs included set")
	ap.add_argument("--enforce-newcomer-diversity", action="store_true", help="Also require newcomers in same dir to be mutually < thr")
	ap.add_argument("--fp-tag", default=FP_TAG_DEFAULT, help="SD tag where Morgan FP bitstring is stored.")
	ap.add_argument("--sim-tag", default=SIM_TAG_DEFAULT, help="SD tag to store similarity vs included set.")
	ap.add_argument("--id-prop", default=None, help="SD property to use as molecule ID (else _Name or auto)")
	ap.add_argument("--radius", type=int, default=RADIUS_DEFAULT, help="Morgan FP radius")
	ap.add_argument("--nbits", type=int, default=NBITS_DEFAULT, help="Morgan FP nBits")
	ap.add_argument("--seed", type=int, default=None, help="Random seed for reproducibility.")
	ap.add_argument("--write-back-fp", action="store_true", help="If a mol lacks FP tag, compute and write it before output")
	ap.add_argument("--scan-recursive", "--use-all-sdfs", action="store_true", help="If set, read all *.sdf files under each first-level subdirectory (recursively) as the binder pool. "
                     "Skips included.sdf, decoys.sdf, and anything under *_decoy/ folders.")
	ap.add_argument("--binders-glob", default="*.sdf", help="Pattern used when --scan-recursive is set (default: *.sdf).")
	ap.add_argument("--exclude-names", nargs="*", default=["included.sdf", "decoys.sdf"], help="Basenames to ignore when --scan-recursive is set.")
	ap.add_argument("--decoys-per-binder", type=int, default=3, metavar="K", help="Cap per-target decoys to K × (#positive binders in that target). Default: 3.")
	ap.add_argument("--included-glob", default=None, help="If set, use this glob (e.g., '*.sdf') recursively to collect included files; overrides --included-name.")
	args = ap.parse_args()

	rng = random.Random(args.seed) if args.seed is not None else random

	root = Path(".").resolve()
	subdirs = sorted([p for p in root.iterdir() if p.is_dir() and not p.name.endswith("_decoy")])

	#Discover protein dirs that contain needed files
	protein_dirs = []
	for d in subdirs:
		protein_dirs.append(d)

	if not protein_dirs:
		print("[ERROR] No protein directories with binders found.", file=sys.stderr)
		sys.exit(1)

	#Build included FP sets per directory
	inc_fps_by_dir: Dict[str, List[DataStructs.ExplicitBitVect]] = {}
	for d in protein_dirs:
		dfps: List[DataStructs.ExplicitBitVect] = []
		
		if getattr(args, "included_glob", None):
			#collect all matching SDFs under this dir (recursively)
			for inc in d.rglob(args.included_glob):
				if not inc.is_file():
					continue
				#skip excludes
				if any(part.endswith("_decoy") for part in inc.parts):
					continue
				if any(inc.name == ex for ex in getattr(args, "exclude_names", [])):
					continue
				_, _, fps = read_sdf_ids_fps(
					inc,
					id_prop=args.id_prop,
					fp_tag=args.fp_tag, radius=args.radius, nBits=args.nbits,
					write_back=False
				)
				dfps.extend(fps)
		else:
			#Fallback: single filename
			included_path = d / (args.included_name if hasattr(args, "included_name") else "included.sdf")
			if included_path.exists():
				_, _, fps = read_sdf_ids_fps(
					included_path,
					id_prop=args.id_prop,
					fp_tag=args.fp_tag, radius=args.radius, nBits=args.nbits,
					write_back=False,
				)
				dfps.extend(fps)

		inc_fps_by_dir[d.name] = dfps	

	#load all binders from all directories
	pool_ids: List[str] = []
	pool_fps: List[DataStructs.ExplicitBitVect] = []
	id_to_mol: Dict[str, Chem.Mol] = {}
	target_counts: Dict[str, int] = {}
	
	global_seen_ids = set()

	def binder_files_in_dir(d: Path) -> List[Path]:
		if not args.scan_recursive:
			binders_name = getattr(args, "binders_name", "binders.sdf")
			p = d / binders_name
			return [p] if p.exists() else []
		files = []
		for p in d.rglob(args.binders_glob):
			if not p.is_file() or any(p.name == ex for ex in args.exclude_names) or any(part.endswith("_decoy") for part in p.parts):
				continue
			files.append(p)

	for d in protein_dirs:
		files = binder_files_in_dir(d)
		if not files:
			print(f"[WARN] {d.name}: no binder SDFs found; skipping directory.")
			target_counts[d.name] = 0
			continue

		binder_ids_by_dir.setdefault(d.name, [])
		per_dir_count = 0
		for f in files:
			ids, mols, fps = read_sdf_ids_fps(f, id_prop=args.id_prop, 
								fp_tag=args.fp_tag, radius=args.radius, nBits=args.nbits,
								write_back=args.write_back_fp)
			binder_ids_by_dir[d.name] = ids

			#prefix IDs with source file so that they are unique and traceable
        # enforce uniqueness across dirs
			for mid, mol, fp in zip(ids, mols, fps):
				gid = f"{d.name}/{f.stem}/{mid}"
				if gid in global_seen_ids:
					k = 1
					while f"{gid}__gpud{k}" in global_seen_ids:
						k += 1
					gid = f"{gid}__gdup{k}"
					#reflect new id on molecule name
					try:
						mol.SetProp("_Name", gid)
					except Exception:
						pass
				global_seen_ids.add(gid)
				pool_ids.append(gid)
				pool_fps.append(fp)
				id_to_mol[gid] = mol
				per_dir_count += 1

		target_counts[d.name] = len(ids)
		print(f"[INFO] {d.name}: {len(ids)} binders; included=(len{inc_fps_by_dir[d.name]})")
	num_binders = {name: len(binder_ids_by_dir.get(name, [])) for name in binder_uds_by_dir}
	if sum(target_counts.values()) != len(pool_ids):
		print("[ERROR] Pool size mismatch.", file=sys.stderr)
		sys.exit(1)

	#we shuffle
	assignment = shuffle_core(pool_ids, pool_fps, inc_fps_by_dir, target_counts,
					thr=args.thr, enforce_newcomer_diversity=args.enforce_newcomer_diversity,
					max_restarts=10,
					rng=rng)
	if assignment is None:
		print("[ERROR] No feasible assignment found under current constraints.", file=sys.stderr)
		sys.exit(2)
	
	#prepare outputs
	assigned_by_dir: Dict[str, List[str]] = defaultdict(list)
	for mid, dname in assignment.items():
		assigned_by_dir[dname].append(mid)
	
	k = args.decoy_per_binder
	cap = k * max(1, num_binders.get(dname, 0))
	if cap and len(chosen_decoys) > cap:
		chosen_decoys = chosen_decoys[:cap]
	#Write {dir}_decoy/decoys.sdf for each dir
	for d in protein_dirs:
		out_dir = d.parent / f"{d.name}_decoy"
		out_dir.mkdir(parents=True, exist_ok=True)
		out_sdf = out_dir / "decoys.sdf"
		w = Chem.SDWriter(str(out_sdf))
		count = 0
		for mid in assigned_by_dir.get(d.name, []):
			mol = id_to_mol.get(mid)
			if mol is None:
				continue
			#set assignment + similarity tag
			try:
				mol.SetProp(ASSIGN_TAG, d.name)
			except Exception:
				pass
			fp = mol_fp_from_tag(mol, fp_tag=args.fp_tag, radius=args.radius, nBits=args.nbits, write_back=False)
			sim = max_sim_to_set(fp, inc_fps_by_dir.get(d.name, [])) if fp is not None else None
			if sim is not None:
				try:
					mol.SetProp(args.sim_tag, f"{sim:.4f}")
				except Exception:
					pass
			#ensure FP tag present if requested
			if args.write_back_fp and not mol.HasProp(args.fp_tag) and fp is not None:
				try:
					mol.SetProp(args.fp_tag, fp.ToBitString())
				except Exception:
					pass
			w.write(mol)
			count += 1
		w.close()
		print(f"[OK] Wrote {count} decoys to {out_sdf}")

	print("[DONE] All decoy sets written.")

if __name__ == "__main__":
	main()
