#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
From RCSB PDB, get receptors -> find in ChEMBL -> scrape all bound ligands (minus RCSB ligand) -> write to individual SDFs

Usage:

What it does:
1) From PDBe, gets PDB entry receptor and converts to ChEMBL target ID
2) For each binding ligand, fetches CanonicalSMILES and InchiKey.
3) Finds ChEMBL molecule record and retrieves existing activity data.
4) Collects pChEMBL values as well as any other available binding data
5) Writes ligands to individual .sdf with metadata stored as SD tags.
"""

import argparse
import json
import os
import re
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Set
import time
from collections import defaultdict
 
from functools import lru_cache
import requests


from rdkit import Chem
from rdkit.Chem import AllChem, rdDepictor, DataStructs

_ACCESSION_RE = re.compile(r"^[A-NR-Z0-9]{6,10}$")

PDBe_UNIPROT_MAP = "https://www.ebi.ac.uk/pdbe/api/mappings/uniprot/{pdbid}"
PDBe_LIGANDS_URL = "https://www.ebi.ac.uk/pdbe/api/pdb/entry/ligand_monomers/{pdbid}"
PDBe_CHEMCOMP_URL = "https://www.ebi.ac.uk/pdbe/api/pdb/chemical_components/{ccd}"
RCSB_CHEMCOMP_URL = "https://data.rcsb.org/rest/v1/core/chemcomp/{ccd}"
RCSB_LIGAND_IDEAL_SDF = "https://files.rcsb.org/ligands/download/{ccd}_ideal.sdf"
RCSB_LIGAND_MODEL_SDF = "https://files.rcsb.org/ligands/download/{ccd}_model.sdf"

CHEMBL_TARGET_BY_ID = "https://www.ebi.ac.uk/chembl/api/data/target/{tid}.json"
CHEMBL_TARGET_REL = "https://www.ebi.ac.uk/chembl/api/data/target_relation.json?{q}&limit=1000"
CHEMBL_TARGET_VIA_UNIPROT = "https://www.ebi.ac.uk/chembl/api/data/target.json?target_components__accession={acc}&limit=200"
CHEMBL_TARGET_ACTIVITIES = "https://www.ebi.ac.uk/chembl/api/data/activity.json?target_chembl_id={tid}&limit=1000&offset={offset}"
CHEMBL_MOL_RECORD = "https://www.ebi.ac.uk/chembl/api/data/molecule.json?molecule_chembl_id={mid}"
CHEMBL_MOL_FROM_INCHIKEY = "https://www.ebi.ac.uk/chembl/api/data/molecule.json?molecule_structures__standard_inchi_key={ikey}"

PUBCHEM_SMILES_CID = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/{smiles}/cids/JSON"
PUBCHEM_CID_PROPS = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/property/CanonicalSMILES,InChIKey/JSON"


def get_json(url: str, retries: int = 3, sleep: float = 0.6) -> Optional[dict]:
	for i in range(retries):
		try:
			r = requests.get(url, timeout=30)
			if r.status_code == 200:
				return r.json()
			if r.status_code in (429, 500, 502, 503, 504):
				time.sleep(sleep * (i + 1))
			else:
				return None
		except requests.RequestException:
			time.sleep(sleep * (i + 1))
	return None

#PDBe helper functions

def access_pdb_uniprot(pdbid: str) -> Set[str]:
	"""
	Return a set of UniProt accessions for a given PDB ID.
	This uses the keys in PDBe's UniProt Mappings directly,
	which are already accession ids in most cases.
	"""
	data = get_json(PDBe_UNIPROT_MAP.format(pdbid=pdbid.lower())) or {}
	entry = data.get(pdbid.lower(), {})
	u = entry.get("UniProt", {}) or {}

	if all(isinstance(v, dict) and "mappings" in v for v in u.values()):
		return set(u.keys())

	accs = set()
	for _chain, items in u.items():
		if isinstance(items, list):
			for it in items:
				if isinstance(it, dict):
					acc = it.get("accession") or it.get("uniprot_acc") or it.get("identifier")
					if acc:
						accs.add(acc)
	if acc.isempty():
		print("Access_PDB_UniProt found no accessions.")
	return accs

def pdb_ligand_inchikey(pdbid: str) -> Set[str]:
	data = get_json(PDBe_LIGANDS_URL.format(pdbid=pdbid.lower()))
	if not data or pdbid.lower() not in data:
		return set()
	inchis: Set[str] = set()
	for lig in data[pdbid.lower()]:
		smiles = lig.get("smiles")
		if not smiles:
			continue
	#get standardized InchiKey
		cid_data = get_json(PUBCHEM_SMILES_CID.format(smiles=requests.utils.quote(smiles, safe="")))
		cids = (cid_data or {}).get("IdentifierList", {}).get("CID", [])
		if not cids:
			continue
		props = get_json(PUBCHEM_CID_PROPS.format(cid=cids[0])) or {}
		try:
			props0 = props["PropertyTable"]["Properties"][0]
			key = props0.get("InChiKey")
			if key:
				inchis.add(key)
		except Exception:
			pass
	return inchis

	#ChEMBL helper functions

def chembl_targets(acc: str) -> List[str]:
	"""Return ChEMBL ID."""
	data = get_json(CHEMBL_TARGET_VIA_UNIPROT.format(acc=acc)) or {}
	out = []
	for t in data.get("targets", []):
		tid = t.get("target_chembl_id")
		if tid:
			out.append(tid)
	return out

def fetch_target_activity_data(target_id: str) -> List[dict]:
	"""Return paginated activity rows for target."""
	rows = []
	offset = 0
	while True:
		data = get_json(CHEMBL_TARGET_ACTIVITIES.format(tid=target_id, offset=offset))
		if not data:
			break
		acts = data.get("activities", [])
		rows.extend(acts)
		if len(acts) < 1000:
			break
		offset += 1000
		time.sleep(0.1)
	return rows

def chembl_molecule_records(mid:str) -> Optional[dict]:
	data = get_json(CHEMBL_MOL_RECORD.format(mid=mid)) or {}
	mols = data.get("molecules") or []
	return mols[0] if mols else None

def chembl_targets_ids_for_subset_complex(uniprots_set: set, acc: str, min_components: int = 1) -> list:
	"""Keep targets where:
		- target_type is a protein complex or protein complex group
		- component UniProt accessions are a non-empty subset
		- and size >= min_components (default 1; set 2 to exclude single subunits)
	"""
	out = []
	for t in chembl_targets_raw_for_uniprot(acc):
		ttype = (t.get("target_type") or "").upper().strip()

		comps = target_component_accessions_from_full(t)
		n = len(comps)
		
		if ttype == "SINGLE PROTEIN":
			if n == 1 and next(iter(comps), None) in uniprots_set and n >= min_components:
				tid = t.get("target_chembl_id")
				if tid:
					out.append(tid)
			continue

		if ttype in ("PROTEIN COMPLEX", "PROTEIN COMPLEX GROUP"):
			if comps and comps.issubset(uniprots_set) and n >= min_components:
				tid = t.get("target_chembl_id")
				if tid:
					out.append(tid)
			continue

	return out

def chembl_targets_raw_for_uniprot(acc: str) -> list:
	"""Return raw ChEMBL target objects for UniProt accession, sans filtering."""
	data = get_json(CHEMBL_TARGET_VIA_UNIPROT.format(acc=acc)) or {}
	return data.get("targets", [])
	
def target_component_accessions_from_full(target_full: dict) -> set:
	"""Extract UniProt accessions from a ChEMBL target."""
	accs: set[str] = set()
	for comp in (target_full.get("target_components") or []):
		acc_single = comp.get("accession")
		if isinstance(acc_single, str) and acc_single:
			accs.add(acc_single)

		#for nested instances under component or target_component
		for k in ("component", "target_component"):
			nested = comp.get(k) or {}
			acc = nested.get("accession")
			if isinstance(acc, str) and acc:
				accs.add(acc)
	return accs

def chembl_targets_ids_for_strict_complex(uniprots_set: set, acc: str) -> list:
	"""Return target IDs where 
		- If the PDB has a single UniProt, allow SINGLE PROTEIN targets with that accession
		- target_type is PROTEIN COMPLEX or PROTEIN COMPLEX GROUP
		- component UniProt accessions exactly match uniprots_set
	"""

	pdb_set = {norm_acc(u) for u in uniprots_set}
	out: list[str] = []

	for t in chembl_targets_raw_for_uniprot(acc):
		tid = t.get("target_chembl_id")
		if not tid:
			continue

		full = fetch_target_by_id(tid)
		if not isinstance(full, dict) or not full:
			continue
		
		ttype = (full.get("target_type") or "").upper().strip()
		comps = {norm_acc(a) for a in target_component_accessions_from_full(full)}
		
		if ttype == "SINGLE PROTEIN" and comps =={acc_norm}:
			out.append(tid)
			continue

		if ttype in ("PROTEIN COMPLEX", "PROTEIN COMPLEX GROUP") and comps == pdb_set:
			out.append(tid)

	return out

#Helpers for selectivity group expansion

def target_relations(q: str) -> list[dict]:
	data = get_json(CHEMBL_TARGET_REL.format(q=q)) or {}
	return data.get("target_relations", []) or []

def selectivity_group_target(tid: str) -> str | None:
	#find SELECTIVITY OF via SUBSET OF operation in ChEMBL
	rows = target_relations(f"target_chembl_id={tid}&relationship=SUBSET%20OF")
	for r in rows:
		gid = r.get("related_target_chembl_id")
		if gid:
			return gid
	return None

def get_selectivity_group_members(gid: str) -> list[str]:
	rows = target_relations(f"related_target_chembl_id={gid}&relationship=SUPERSET%20OF")
	return [r.get("target_chembl_id") for r in rows if r.get("target_chembl_id")]

def fetch_target_by_id(tid: str) -> dict:
	data = get_json(CHEMBL_TARGET_BY_ID.format(tid=tid)) or {}
	ts = data.get("targets") or []
	return ts[0] if ts else {}


def norm_acc(x: str) -> str:
	return (x or "").upper().split("-")[0]


#Helpers for Tanimoto similarity

def pdb_ligand_smiles(pdbid: str, ligand_id: str) -> str:
	"""
    Resolve a PDB ligand's SMILES. Order:
      1) PDBe ligand_monomers (per-entry)
      2) RCSB chemcomp (descriptors)
      3) PDBe chemical components
      4) RCSB ligand SDF (ideal → model) → RDKit
    Returns '' if unresolved.
    """
	if not pdbid or not ligand_id:
		return ""
	pid = pdbid.lower().strip()
	ccd = ligand_id.upper().strip()

	#PDBe per-entry approach
	try:
		data = get_json(PDBe_LIGANDS_URL.format(pdbid=pdbid.lower())) or {}
		ligs = (data.get(pid) or [])
		for lig in ligs:
			if (lig.get("chem_comp_id") or "").upper() == ccd:
				smi = lig.get("smiles") or lig.get("modified_smiles") or ""
				if smi:
					return smi
	except Exception:
		pass

	#RCSB chemcomp
	try:	
		r = requests.get(RCSB_CHEMCOMP_URL.format(ccd=ccd), timeout=15)
		if r.ok:
			j = r.json() or {}
			desc = j.get("rcsb_chem_comp_descriptor") or {}
			smi = desc.get("smiles") or desc.get("smiles_stereo") or desc.get("smiles_derived") or ""
			if not smi:
				for ent in (j.get("pdbx_chem_comp_descriptor") or []):
					if "smiles" in ent.get("type"or "").lower():
						smi = ent.get("descriptor") or ""
						if smi:
							break
			if smi:
				return smi
	
	except Exception:
		pass

	#PDBe chem components
	try:
		comp = get_json(PDBe_CHEMCOMP_URL.format(ccd=ccd.lower())) or {}
		rows = comp.get(ccd.lower())
		if rows and isinstance(rows, list):
			smi = rows[0].get("smiles") or rows[0].get("modified_smiles") or ""
			if smi:
				return smi
	except Exception:
		pass

	#Model PDB -> SMILES
	try:
		r = requests.get(RCSB_LIGAND_MODEL_PDB_URL.format(ccd=ccd), timeout=15)
		if r.ok and r.text:
			mol = Chem.MolFromPDBBlock(r.text, removeHs=False)
			if mol:
				return Chem.MolToSmiles(mol)
	except Exception:
		pass

	#if nothing found
	return ""

def tanimoto_smiles(a: str, b: str) -> float:
	"""Morgan FP (radius=2, 2048 bits) Tanimoto similarity for two SMILES; returns 0.0 if either fails."""
	try:
		ma = Chem.MolFromSmiles(a) if a else None
		mb = Chem.MolFromSmiles(b) if b else None
		if not ma or not mb:
			return 0.0
		fa = AllChem.GetMorganFingerprintAsBitVect(ma, radius=2, nBits=2048)
		fb = AllChem.GetMorganFingerprintAsBitVect(mb, radius=2, nBits=2048)
		return float(DataStructs.TanimotoSimilarity(fa, fb))
	except Exception:
		return 0.0

def smiles_from_sdf_text(sdf_text: str) -> str:
    try:
        suppl = Chem.SDMolSupplier()
        suppl.SetData(sdf_text, removeHs=False)
        mols = [m for m in suppl if m]
        if mols:
            return Chem.MolToSmiles(mols[0])
    except Exception:
        pass
    return ""

#Data retrieval helpers

def summary_of_activity_data(acts: List[dict]) -> Tuple[List[float], Dict[str, dict]]:
	"""Return a tuple of pChEMBL values and a summary by type of data (ie IC50)."""
	pvals: List[float] = []
	buckets: Dict[str, List[Tuple[float, str]]] = defaultdict(list)

	for a in acts:
		p = a.get("pchembl_value")
		if p:
			try:
				pvals.append(float(p))
			except ValueError:
				pass

		stype = (a.get("standard_type") or a.get("type") or "").upper().replace(" ", "")
		sval = a.get("standard_value")
		sunit = a.get("standard_units")
		try:
			val = float(sval) if sval is not None else None
		except ValueError:
			val = None
		if stype and (val is not None):
			buckets[stype].append((val, sunit))

	def summ(vals_units):
		vals = [v for v, _ in vals_units]
		if not vals:
			return {}
		vals.sort()
		n = len(vals)
		mid = n // 2
		median = vals[mid] if n % 2 else 0.5 * (vals[mid - 1] + vals[mid])
		units = [u for v, u in vals_units if u]
		return {"n": n, "min": vals[0], "median": median, "max": vals[-1], "units": units}
	
	summary: Dict[str, dict] = {}
	for k, v in buckets.items():
		summary[k] = summ(v)

	return pvals, summary

	#RDkit writers

def write_sdf(smiles: str, out_path: Path, props: Dict[str, str]) -> bool:
	try:
		if not smiles:
			return false

		mol = Chem.MolFromSmiles(smiles)
		if mol is None:
			return False
		
		mol = Chem.AddHs(mol)
		try:
			params = AllChem.ETKDGv3() if hasattr(AllChem, "ETKDGv3") else AllChem.ETKDG()
			if AllChem.EmbedMolecule(mol, params) >= 0:
				AllChem.UFFOptimizeMolecule(mol)
			else:
				rdDepictor.Compute2DCoords(mol)
		except Exception:
			redDepictor.compute2DCoords(mol)

		for k, v in (props or {}).items():
			if v is not None:
				mol.SetProp(str(k), str(v))

		w = Chem.SDWriter(str(out_path))
		try:
			w.write(mol)
		finally:
			w.close()
		return True
	except Exception:
		return False


	#main workflow

def dbg_targets_for_acc(acc: str):
	print(f"[DBG] Raw ChEMBL targets for {acc}:")
	ts = chembl_targets_raw_for_uniprot(acc)
	print(f"      count={len(ts)}")
	for t in ts[:20]:
		tid = t.get("target_chembl_id")
		ttype = (t.get("target_type") or "").upper().strip()
		comps = sorted(target_component_accessions_from_full(t))
		print(f"	{tid} [{ttype}] comps={comps}")

def main():
	ap = argparse.ArgumentParser(description="PDB -> ChEMBL target ligands (excluding the PDB's ligand) -> SDFs")
	ap.add_argument("--pdb", required=True, help="RCSB PDB ID, eg 8C1Y")
	ap.add_argument("--out", required=True, help="Output directory for SDF files")
	ap.add_argument("--sleep", type=float, default=0.2, help="Sleep between API calls (s)")
	ap.add_argument("--subset-complex", action="store_true", help="Keep complexes whose components are a subset of the PDB UniProt set.")
	ap.add_argument("--min-components", type=int, default=1, help="Minimum number of UniProt components required when using --subset-complex (default: 1).")
	ap.add_argument("--strict-complex", action="store_true", help="Only keep ChEMBL targets that are PROTEIN COMPLEX (/GROUP) whose UniProt components exactly match the PDB's UniProt set.")
	ap.add_argument("--dry-run", action="store_true", help="Print chosen targets and exit before fetching activities.")
	ap.add_argument("--by-selectivity-group", action="store_true", help="After selecting targets (broad/subset/strict), expand to all targets in the same selectivity group(s).")
	ap.add_argument("--tanimoto-percentage", type=float, default=None, help="Keep only ligands with Tanimoto similarity >= this threshold (0-1) to the provided PDB ligand.")
	ap.add_argument("--tanimoto-ligand", type=str, default=None, help="PDB ligand identifier to compare against (i.e. K, ATP, EDO). Requires --tanimoto-percentage.")
	args = ap.parse_args()

	pdbid = args.pdb.strip()
	outdir = Path(args.out)
	outdir.mkdir(parents=True, exist_ok=True)

	print(f"[INFO] Mapping PDB {pdbid} to UniProt...")

	uniprots = access_pdb_uniprot(pdbid)
	if not uniprots:
		print("[WARN] No UniProt mappings found.")
		return
	print(f"	UniProt: {', '.join(sorted(uniprots))}")


	print("[INFO] Resolving ChEMBL targets ...")
	#Call debugger
	for acc in sorted(uniprots):
		dbg_targets_for_acc(acc)
	
	target_ids = set()
	if args.strict_complex:
		#precise match
		for acc in sorted(uniprots):
			target_ids.update(chembl_targets_ids_for_strict_complex(uniprots, acc))
	elif args.subset_complex:
		#subset match
		for acc in sorted(uniprots):
			target_ids.update(chembl_targets_ids_for_subset_complex(uniprots, acc, args.min_components))
	else:
		for acc in sorted(uniprots):
			tids = chembl_targets_raw_for_uniprot(acc)
			target_ids.update(tids)
	print(f"		Targets: {', '.join(sorted(target_ids))}")
	
	if args.by_selectivity_group and target_ids:
		print("[INFO] Expanding by selectivity group(s) ... ")
		expanded = set(target_ids)
		seen_groups = set()
		for tid in sorted(target_ids):
			gid = selectivity_group_target(tid)
			if gid and gid not in seen_groups:
				seen_groups.add(gid)
				expanded.update(get_selectivity_group_members(gid))
				time.sleep(args.sleep * 0.5)
		target_ids = expanded
		print(f"	After selectivity-group expansion: {len(target_ids)} targets")
	
	if args.dry_run:
		print("[DRY] Selected target IDs:")
		for tid in sorted(target_ids):
			if args.by_selectivity_group:
				full = fetch_target_by_id(tid)
				ttype = (full.get("target_type") or "").upper().strip()
				sg = selectivity_group_target(full)
				print(f"	{tid} [{ttype}] selgrp={sg}")
			else:
				print(" ", tid)
		return

	print("[INFO] Excluding PDB ligands ... ")
	pdb_inchikeys = pdb_ligand_inchikey(pdbid)
	print(f"	Excluding {len(pdb_inchikeys)} ligands by InChiKey.")

	#Collect molecules from activities for target
	print("[INFO] Fetching activities for targets and collecting ligands ... ")
	mol_to_acts: Dict[str, List[dict]] = defaultdict(list)
	target_to_mols: Dict[str, str] = defaultdict(set)

	for tid in sorted(target_ids):
		acts = fetch_target_activity_data(tid)
		for a in acts:
			mid = a.get("molecule_chembl_id")
			if mid:
				mol_to_acts[mid].append(a)
				target_to_mols[tid] .add(mid)
		print(f"	{tid}: +{len(acts)} activity rows, {len(mol_to_acts)} molecules so far.")
		time.sleep(args.sleep)

	molecule_cache: Dict[str, str] = {}
	smiles_cache: Dict[str, str] = {}
	inchikey_cache: Dict[str, str] = {}

	for mid in mol_to_acts.keys():
		if mid in molecule_cache:
			continue
		rec = chembl_molecule_records(mid)
		if not rec:
			continue
		molecule_cache[mid] = rec
		struct = rec.get("molecule_structures") or {}
		smiles_cache[mid] = struct.get("canonical_smiles") or struct.get("standard_smiles") or ""
		

	#Filter PDB ligands
	keep_mols: Dict[str, List[dict]] = {}
	for mid, acts in mol_to_acts.items():
		rec = molecule_cache.get(mid)
		if not rec:
			continue
		inchi_key = inchikey_cache.get(mid, "")
		if inchi_key and inchi_key in pdb_inchikeys:
			continue
		if not smiles_cache.get(mid):
			continue
		keep_mols[mid] = acts
		time.sleep(0.05)

	for tid in list(target_to_mols.keys()):
		target_to_mols[tid] = {mid for mid in target_to_mols[tid] if mid in keep_mols}
		if not target_to_mols[tid]:
			del target_to_mols[tid]

	print(f"[INFO] {len(keep_mols)} molecules remain after exclusion.")
	

	if args.tanimoto_percentage is not None:
		thr = float(args.tanimoto_percentage)
		if not (0.0 <= thr <= 1.0):
			print(f"[WARN] --only-tanimotos must be between 0 and 1; got {thr}. Ignoring.")
		elif not args.tanimoto_ligand:
			print("[WARN] --only-tanimotos provided without --tanimoto-ligand; ignoring.")
		else:
			ref_smiles = pdb_ligand_smiles(pdbid, args.tanimoto_ligand)
			if not ref_smiles:
				print(f"[WARN] Couldn’t find SMILES for ligand '{args.tanimoto_ligand}' in {pdbid}; ignoring Tanimoto filter.")	
			else:
				ref_mol = Chem.MolFromSmiles(ref_smiles)
				if not ref_mol:
					print(f"[WARN] Bad reference SMILES for '{args.tanimoto_ligand}'; skipping Tanimoto filter.")
				else:
					ref_fp = AllChem.GetMorganFingerprintAsBitVect(ref_mol, radius=2, nBits=2048)
					
					mids_list = []
					fp_list = []
					for mid in keep_mols.keys():
						smi = smiles_cache.get(mid, "")
						m = Chem.MolFromSmiles(smi) if smi else None
						if not m:
							continue
						fp_list.append(AllChem.GetMorganFingerprintAsBitVect(m, radius=2, nBits=2048))
						mids_list.append(mid)

					if not fp_list:
						print("[WARN] No candidate fingerprints; skipping Tanimoto filter.")
					else:
						sims = DataStructs.BulkTanimotoSimilarity(ref_fp, fp_list)
						kept_by_similarity: Dict[str, List[dict]] = {}
						tanimaps: Dict[str, float] = {}
						fp_maps: Dict[str, str] = {}
						for mid, sim, fp in zip(mids_list, sims, fp_list):
							if sim >= thr:
								kept_by_similarity[mid] = keep_mols[mid]
								tanimaps[mid] = sim
								fp_maps[mid] = fp.ToBitString()

						removed = len(keep_mols) - len(kept_by_similarity)
						keep_mols = kept_by_similarity
						print(f"[INFO] Tanimoto filter ≥ {thr:.2f} vs {args.tanimoto_ligand}: kept {len(keep_mols)}, removed {removed}.")


						#prune target_to_mols again to match filtered set
						for tid in list(target_to_mols.keys()):
							target_to_mols[tid] = {mid for mid in target_to_mols[tid] if mid in keep_mols}
							if not target_to_mols[tid]:
								del target_to_mols[tid]

						#stash similarity and the Morgan Fingerprint  as SD tags
						for mid, acts in keep_mols.items():
							if acts:
								acts[0]["__tanimoto_sim__"] = tanimaps.get(mid)
								acts[0]["__morgan_fp__"] = fp_maps.get(mid)				
	
	#Write SDFs
	index: List[dict] = []

	for i, (tid, mol_ids) in enumerate(sorted(target_to_mols.items()), 1):
		sdf_path = outdir / f"{pdbid}_{tid}.sdf"
		n_written = 0
		w = Chem.SDWriter(str(sdf_path))
		try:
			for mid in sorted(mol_ids):
				rec = chembl_molecule_records(mid)
				if not rec:
					continue

				# Structural properties
				struct = rec.get("molecule_structures") or {}
				can_smiles = struct.get("canonical_smiles") or struct.get("standard_smiles") or ""
				inchi_key = struct.get("standard_inchi_key") or ""
				pref_name = rec.get("pref_name") or ""
				if not can_smiles:
					continue

				# Use per-molecule activities, not the old 'acts' loop variable
				acts_for_mid = keep_mols.get(mid, [])  # or mol_to_acts[mid]
				pvals, act_summary = summary_of_activity_data(acts_for_mid)

				mol = Chem.MolFromSmiles(can_smiles)
				if mol is None:
					continue
				mol = Chem.AddHs(mol)
				try:
					params = AllChem.ETKDGv3() if hasattr(AllChem, "ETKDGv3") else AllChem.ETKDG()
					if AllChem.EmbedMolecule(mol, params) >= 0:
						AllChem.UFFOptimizeMolecule(mol)
				except Exception:
					pass

				# SD tags: similarity + fingerprint (bitstring)
				sim_val, fp_val = None, None
				rec_keep = keep_mols.get(mid)
				if rec_keep and isinstance(rec_keep[0], dict):
					sim_val = rec_keep[0].get("__tanimoto_sim__")
					fp_val  = rec_keep[0].get("__morgan_fp__")  # should already be a bitstring

				# set props
				if sim_val is not None and args.tanimoto_ligand:
					mol.SetProp(f"Tanimoto_to_{args.tanimoto_ligand.upper()}", f"{sim_val:.3f}")
				if fp_val:
					mol.SetProp("__morgan_fp__", str(fp_val))

				mol.SetProp("PDB_ID", pdbid)
				mol.SetProp("Molecule_ChEMBL_ID", mid)
				mol.SetProp("Preferred Name", pref_name)
				mol.SetProp("InChIKey", inchi_key)  # standard capitalization
				if pvals:
					mol.SetProp("pChEMBL_values", ",".join(f"{v:.3f}" for v in pvals))
				mol.SetProp("Activity_Summary_JSON", json.dumps(act_summary, sort_keys=True))

				w.write(mol)
				n_written += 1

				index.append({
					"pdb_id": pdbid,
					"molecule_chembl_id": mid,
					"preferred_name": pref_name,
					"canonical_smiles": can_smiles,
					"inchi_key": inchi_key,
					"pchembl_values": pvals,
					"activity_summary": act_summary,
					"sdf_file": sdf_path.name,
						})
		finally:
			w.close()

		if n_written:
			print(f"[{i}/{len(target_to_mols)}] Wrote {sdf_path.name} with {n_written} molecules")
		else:
			if sdf_path.exists():
				sdf_path.unlink()
			print(f"[{i}/{len(target_to_mols)}] Skipped {tid} (no valid molecules)")
			time.sleep(0.05)

	with open(outdir / f"{pdbid}_ligand_index.json", "w") as f:
		json.dump(index, f, indent=2, sort_keys=True)
	print(f"[DONE] {len(index)} molecules processed. Index written to {outdir / (pdbid + '_ligand_index.json')}")
	

if __name__ == "__main__":
	main()


