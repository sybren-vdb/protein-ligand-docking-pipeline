#!/usr/bin/env python3
"""Convert a Gnina multi-pose SDF into per-pose mmcif receptor/ligand complexes.

Example:
pymol -cq sdf_poses_to_mmcif.py -- \
    --pdb-id 7PAV \
    --ligand-resn 6IO \
    --sdf 7PAV_6IO_poses.sdf \
    --outdir ./out \
    --remove-waters \
    --protein-only

Notes:
	Requires PyMol (headless is ok). Run with 'pymol -cq ... --<args>.
	By default, we fetch the receptor from RCSB as mmCIF..
	We strip only the ligand residue specified (eg 6IO) from the receptor before merging in each SDF pose.
"""

import os
import argparse
from pathlib import Path

from pymol import cmd
import pymol2
import sys
import traceback


def parse_args():
	p = argparse.ArgumentParser(description="Build per-pose mmCIFs from a Gnina SDF.")
	p.add_argument("--pdb-id", required=True, help="RCSB ID of the receptor (e.g., 7PAV)")
	p.add_argument("--sdf", required=True, help="Path to Gnina output SDF (multi-pose)")
	p.add_argument("--ligand-resn", required=True, help="Residue name to remove from receptor and assign to poses.")
	p.add_argument("--outdir", default="mmcif_out", help="Output directory for mmcif files")
	p.add_argument("--prefix", default=None, help="Filename prefix (default: <PDB>_<RESN>)")
	p.add_argument("--protein-only", action="store_true", help="Keep only polymeric protein from receptor")
	p.add_argument("--remove-waters", action="store_true", help="Remove solvent from receptor")
	p.add_argument("--keep-het", action="store_true", help="Keep receptor HET groups (default: keep; combine with --protein-only to drop them)")
	p.add_argument("--pose-resi-start", type=int, default=1001, help="Residue index to assign to ligand poses")
	p.add_argument("--biounit", type=int, default=None, help="Experimental: expand BIOMT to specified assembly number (PyMOL-biomt).")
	return p.parse_args()

def _die(msg):
	print(f"[FATAL] {msg}", file=sys.stderr, flush=True)

def _check(sel, tag):
	n = cmd.count_atoms(sel)
	print(f"[CHECK] {tag}: {n}", flush=True)
	if n == 0:
		_die(f"Empty selection: {tag}")
	return n

def ensure_outdir(path):
	Path(path).mkdir(parents=True, exist_ok=True)

def fetch_receptor(pdb_id, obj_name="receptor", biounit=None):
	#fetch as mmCIF
	cmd.set("fetch_path", os.getcwd())
	#turn off async to ensure data is available
	cmd.fetch(pdb_id, name=obj_name, type="cif", async_=0)
	#Optionally expand bio assembly
	if biounit is not None:
		try:
			cmd.biomt(obj_name, biounit)
		except Exception:
			print(f"Warning: Failed to apply BIOMT assembly {biounit} for {obj_name}")

def clean_receptor(obj_name, ligand_resn, protein_only=False, remove_waters=False, keep_het=True):
	#Remove specified ligand
	cmd.remove(f"{obj_name} and resn {ligand_resn}")
	if protein_only:
		#keep only polymeric protein atoms
		cmd.select("prot_sel", f"{obj_name} and polymer.protein")
		cmd.create(f"{obj_bname}_prot", "prot_sel")
		cmd.delete(obj_name)
		cmd.set_name(f"{obj_name}_prot", obj_name)
		cmd.delete("prot_sel")
		#if protein_only, HET groups are removed by definition
	else:
		if not keep_het:
			cmd.remove(f"{obj_name} and hetatm and not resn {ligand_resn}")

	if remove_waters:
		cmd.remove(f"{obj_name} and solvent")

	#Standardize
	cmd.alter(obj_name, "chain=chain if chain else 'A'")
	cmd.sort(obj_name)
	
def load_poses(sdf_path, obj_name="poses"):
	#load sdf
	existing = set(cmd.get_names("objects"))
	if obj_name in existing:
		cmd.delete(obj_name)
	cmd.load(sdf_path, obj_name)
	n_states = cmd.count_states(obj_name)
	if n_states < 1:
		raise RuntimeError("No poses found in SDF.")
	return n_states

def prepare_pose_object(src_obj, state_idx, ligand_resn, resi):
	#Create a new single state object for this pose
	pose_obj = f"lig_{state_idx:03d}"
	cmd.create(pose_obj, src_obj, state_idx, 1)
	#normalize resi naming
	cmd.alter(pose_obj, f"resn='{ligand_resn}'")
	cmd.alter(pose_obj, f"resi='{resi}'")
	cmd.alter(pose_obj, "chain='L'")
	cmd.sort(pose_obj)
	return pose_obj

def save_complex_mmcif(receptor_obj, pose_obj, out_path, state_idx):
	cmd.save(out_path, f"({receptor_obj}) or ({pose_obj})", state=state_idx,format="cif")

def main():
	args = parse_args()
	ensure_outdir(args.outdir)
	prefix = args.prefix or f"{args.pdb_id}_{args.ligand_resn}"

	pm = pymol2.PyMOL()
	pm.start()

	try:
		globals()["cmd"] = pm.cmd
	
		receptor_obj = "receptor"
		fetch_receptor(args.pdb_id, receptor_obj, biounit=args.biounit)
		clean_receptor(
			receptor_obj,
			ligand_resn=args.ligand_resn,
			protein_only=args.protein_only,
			remove_waters=args.remove_waters,
			keep_het=args.keep_het
			)

		_check("receptor", "receptor (post-load/cleanup)")
	
		print("[INFO] SDF path:", os.path.abspath(args.sdf), flush=True)
		n_states= load_poses(args.sdf, obj_name="poses")
		print(f"Found {n_states} pose(s) in {args.sdf}")

		#Iterate poses + write mmCIFs
		for i in range(1, n_states + 1):
			lig_obj = prepare_pose_object(
					"poses", i, args.ligand_resn, args.pose_resi_start + (i - 1)
					)
			sel = f"({receptor_obj}) or ({lig_obj})"
			cmd.save(out_path := os.path.join(args.outdir, f"{prefix}_pose{i:02d}.cif"),
					sel, state=1, format="cif")
			print(f"[WRITE] {out_path}", flush=True)
			cmd.delete(lig_obj)

		#Cleanup
		cmd.delete("poses")
		print("Done.")
	
	finally:
		pm.stop()


if __name__ == "__main__":
	try:
		main()
	except SystemExit:
		raise
	except Exception:
		print("[EXCEPTION] Uncaught exception:", flush=True)
		traceback.print_exc()
		sys.exit(1)
