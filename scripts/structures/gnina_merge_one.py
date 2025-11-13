#!/usr/bin/env python3
import os, argparse
import gemmi

WATERS = {'HOH', 'WAT', 'H2O', 'DOD'}

def count_atoms_in_block(block: gemmi.cif.Block) -> int:
    try:
        st = gemmi.make_structure_from_block(block)
        return st[0].atom_count() if len(st) > 0 else 0
    except Exception:
        return 0

def guess_blocks(doc: gemmi.cif.Document):
    # Prefer names when available; otherwise fall back to size (max = receptor, min = ligand)
    rec_i = lig_i = None
    rec_atoms, lig_atoms = -1, 10**9
    for i, blk in enumerate(doc):
        n = count_atoms_in_block(blk)
        if n == 0:
            continue
        name = (blk.name or "").lower()
        # strong hints
        if "receptor" in name and n > rec_atoms:
            rec_i, rec_atoms = i, n
            continue
        if ("lig" in name or "pose" in name) and n < lig_atoms:
            lig_i, lig_atoms = i, n
            continue
        # generic size heuristics
        if n > rec_atoms:
            rec_i, rec_atoms = i, n
        if n < lig_atoms:
            lig_i, lig_atoms = i, n

    if rec_i is None or lig_i is None or rec_i == lig_i:
        raise RuntimeError("receptor/ligand block not found")
    return rec_i, lig_i

def structure_from_block(block: gemmi.cif.Block) -> gemmi.Structure:
    st = gemmi.make_structure_from_block(block)
    if len(st) == 0:
        raise RuntimeError("block has no models/atoms")
    return st

def strip_experimental_hets(model: gemmi.Model, keep_water=True):
    # Remove non-polymer residues (keep waters if requested)
    for ch in list(model):
        for i in range(len(ch.residues)-1, -1, -1):
            r = ch.residues[i]
            if r.het_flag:
                if keep_water and r.name in WATERS:
                    continue
                del ch.residues[i]

def add_ligand_as_chain(model: gemmi.Model, lig_st: gemmi.Structure, chain_id='Z'):
    # Pick an unused chain ID
    used = {c.name for c in model}
    cid = chain_id if chain_id not in used else None
    if cid is None:
        for c in "ZYXWVUTSRQPONMLKJIHGFEDCBA":
            if c not in used:
                cid = c
                break
    ch = gemmi.Chain(cid)

    lig_model = lig_st[0]
    seq = 1
    for lch in lig_model:
        for lres in lch:
            new_res = gemmi.Residue()
            new_res.name = lres.name or "LIG"
            new_res.het_flag = True
            new_res.seqid = gemmi.SeqId(seq); seq += 1
            for lat in lres:
                a = gemmi.Atom()
                a.name = lat.name
                a.element = lat.element
                a.pos = gemmi.Position(lat.pos.x, lat.pos.y, lat.pos.z)
                a.occ = lat.occ
                a.b_iso = lat.b_iso
                new_res.add_atom(a)
            ch.add_residue(new_res)

    model.add_chain(ch)

def write_single_block_cif(structure: gemmi.Structure, out_path: str, block_name: str):
    doc = structure.make_mmcif_document()
    if block_name:
        doc[0].name = block_name
    with open(out_path, "w") as f:
        f.write(doc.as_string())

def build_complex(gnina_cif_path: str, assembly_cif_path: str, out_path: str):
    gdoc = gemmi.cif.read(gnina_cif_path)
    rec_i, lig_i = guess_blocks(gdoc)

    # Receptor: assembly if provided, else receptor block from Gnina CIF
    if assembly_cif_path and os.path.isfile(assembly_cif_path):
        st = gemmi.read_structure(assembly_cif_path)
        if len(st) == 0:
            raise RuntimeError("assembly structure has no models")
    else:
        st = structure_from_block(gdoc[rec_i])

    m = st[0]
    strip_experimental_hets(m, keep_water=True)

    # Ligand: build from ligand block Structure
    lig_st = structure_from_block(gdoc[lig_i])
    if lig_st[0].atom_count() == 0:
        raise RuntimeError("ligand block has no atoms")
    add_ligand_as_chain(m, lig_st, chain_id='Z')

    base = os.path.splitext(os.path.basename(gnina_cif_path))[0]
    write_single_block_cif(st, out_path, block_name=f"data_{base}")

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gnina", required=True, help="input Gnina mmCIF (multi-block)")
    ap.add_argument("--assembly", default="", help="RCSB assembly mmCIF for the receptor (optional)")
    ap.add_argument("--out", required=True, help="output mmCIF (single block, single model)")
    args = ap.parse_args()

    try:
        build_complex(args.gnina, args.assembly, args.out)
        print(f"[OK] {args.gnina} -> {args.out}")
    except Exception as e:
        print(f"[SKIP] {args.gnina}\n{type(e).__name__}: {e}")

if __name__ == "__main__":
    main()
