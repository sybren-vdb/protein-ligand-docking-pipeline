# Protein–Ligand Docking Pipeline (Gnina + ChEMBL + RCSB)

This repository contains a collection of scripts for building 
protein–ligand docking datasets and benchmarking pipelines.

It is centered around Gnina-based docking, with utilities for:

- Fetching and preprocessing protein structures from the PDB (mmCIF / assemblies)
- Scraping ligands and activity data from ChEMBL and PDBe/RCSB
- Generating ligand decoys via Tanimoto similarity (RDKit)
- Running Gnina docking (optionally inside Apptainer)
- Converting Gnina pose SDFs into receptor–ligand mmCIF complexes
- Exporting top-k docking poses and affinities as JSONL for ML workflows

## Layout

- `scripts/structures/`
  - `sample_pdbs_to_cif.py`: sample and download mmCIF or assembly files from RCSB.
  - `sdf_poses_to_mmcif.py`: use PyMOL to convert multi-pose SDFs into per-pose mmCIF complexes.
  - `merge_to_assembly.py`: merge Gnina output back into biological assemblies using Gemmi.
  - `merge_lig_rec_blocks.py`: merge receptor and ligand CIF blocks into a single model.
  - `gnina_merge_one.py`: single-step receptor/ligand merge, strip waters/hets, write mmCIF.

- `scripts/ligands/`
  - `scrapeChemBLMols.py`: given a PDB ID, find corresponding ChEMBL targets and scrape ligands + activities into SDF.
  - `tanimotoShuffle.py`: compute RDKit Morgan fingerprints, Tanimoto similarities, and shuffle / assign ligands for decoys.

- `scripts/docking/`
  - `gninaDocking.py`: high-level docking script; prepares inputs, calls `gnina`, and manages intermediate files.
  - `gninaOnly.py`: thin wrapper around `gnina` (optionally via Apptainer) for direct docking runs.
  - `gnina_topk_to_affinity_map.py`: parse Gnina SDF outputs, select top-k poses per ligand, convert pK ↔ ΔG, and write a JSONL
    mapping suitable for training downstream models (e.g., Boltz2).

## Dependencies

This code uses:

- Python 3.9+
- [RDKit](https://www.rdkit.org/)
- [Gemmi](https://gemmi.readthedocs.io/)
- [PyMOL](https://pymol.org/) / `pymol2` for some structure operations
- `requests` for HTTP API calls (RCSB, PDBe, ChEMBL)
- External tools:
  - [`gnina`](https://github.com/gnina/gnina)
  - [`obabel` / Open Babel] for some format conversions

A minimal `environment.yml` is provided as a starting point.

## Example usage

1. Fetch structures and ligands
   - Use `sample_pdbs_to_cif.py` to download mmCIF / assemblies.
   - Use `scrapeChemBLMols.py` to fetch ligands and activity data from ChEMBL.

2. Generate decoys
   - Use `tanimotoShuffle.py` to construct decoy / shuffled ligand sets.

3. Run docking
   - Use `gninaDocking.py` or `gninaOnly.py` to run Gnina on the prepared receptor/ligand pairs.

4. Postprocess results
   - Use `sdf_poses_to_mmcif.py`, `merge_to_assembly.py`, or `gnina_merge_one.py` to build receptor–ligand complexes.
   - Use `gnina_topk_to_affinity_map.py` to export top-k poses and affinity labels as JSONL.

Each script has `--help` describing its arguments and options.
