#!/usr/bin/env python3
# merge_gnina_to_complex_v2.py
import os, re, glob, argparse
import gemmi
import multiprocessing as mp

def structure_from_block(block):
    try:
        s = gemmi.make_structure_from_block(block)
        # normalize: ensure at least one model
        if len(s) == 0:
            s.add_model("1")
        return s
    except Exception:
        return None

def count_atoms(struct):
    if not struct or len(struct) == 0:
        return 0
    n = 0
    for ch in struct[0]:
        for r in ch:
            n += len(r)
    return n

def choose_blocks(doc):
    # Build structures for all blocks; keep only those with atoms
    candidates = []
    for b in doc:
        s = structure_from_block(b)
        n = count_atoms(s)
        if n > 0:
            candidates.append((b, s, n))
    if len(candidates) < 2:
        return None, None  # can't find both receptor and ligand
    # receptor = block with most atoms; ligand = any other non-empty block
    candidates.sort(key=lambda x: x[2], reverse=True)
    rec_b, rec_s, _ = candidates[0]
    # pick the smallest non-empty as ligand (works for 2-block Gnina mmCIFs)
    lig_b, lig_s, _ = candidates[-1]
    return (rec_b, rec_s), (lig_b, lig_s)

def unique_chain_name(model, base='L'):
    used = {c.name for c in model}
    if base not in used:
        return base
    i = 1
    while True:
        name = f'{base}{i}'
        if name not in used:
            return name
        i += 1

def strip_experimental_hets(model, keep_water=True):
    water_names = {'HOH', 'WAT', 'H2O'}
    # iterate over a snapshot of chains so we can modify in place
    for ch in list(model):
        # collect residue indices to delete
        to_del = []
        for i, res in enumerate(ch):
            if res.het_flag == 'H':
                if keep_water and res.name in water_names:
                    continue
                to_del.append(i)
        # delete from the end to preserve indices
        for i in reversed(to_del):
            del ch[i]

def copy_residue_atoms(src_res, dst_res):
    for a in src_res:
        na = gemmi.Atom()
        na.name = a.name
        na.element = a.element
        na.occ = a.occ
        na.b_iso = a.b_iso
        na.altloc = a.altloc
        na.pos = a.pos
        dst_res.add_atom(na)

def add_structure_as_chain(src_struct, dst_model, chain_name):
    ch_out = gemmi.Chain(chain_name)
    dst_model.add_chain(ch_out)
    for src_ch in src_struct[0]:
        for src_res in src_ch:
            r = gemmi.Residue()
            r.name = src_res.name
            # safest way to copy seqid across versions
            r.seqid = gemmi.SeqId(str(src_res.seqid))
            r.het_flag = src_res.het_flag
            copy_residue_atoms(src_res, r)
            ch_out.add_residue(r)
    return ch_out

def load_receptor_container(asm_path, rec_struct):
    # Prefer assembly1 if present; else use receptor struct as is
    if asm_path and os.path.exists(asm_path):
        try:
            s = gemmi.read_structure(asm_path)
            if len(s) == 0:
                s.add_model("1")
            return s
        except Exception:
            pass
    return rec_struct

def merge_one(gnina_cif, assemblies_dir, out_dir, strip_hets=True, debug=False):
    try:
        doc = gemmi.cif.read_file(gnina_cif)
    except Exception as e:
        print(f"[ERROR] {os.path.basename(gnina_cif)}: CIF read failed: {e}")
        return False

    rec, lig = choose_blocks(doc)
    if not rec or not lig:
        print(f"[SKIP] {os.path.basename(gnina_cif)}: receptor/ligand block not found")
        return False

    rec_block, rec_struct = rec
    lig_block, lig_struct = lig

    # PDB id from filename for assembly lookup
    m = re.search(r'([0-9][A-Z0-9]{3})', os.path.basename(gnina_cif).upper())
    pdbid = m.group(1) if m else None
    asm_path = os.path.join(assemblies_dir, f"{pdbid}_assembly1.cif") if pdbid else None

    # receptor container
    try:
        out_struct = load_receptor_container(asm_path, rec_struct)
    except Exception as e:
        print(f"[ERROR] {os.path.basename(gnina_cif)}: receptor load failed: {e}")
        return False

    # ensure a model
    if len(out_struct) == 0:
        out_struct.add_model("1")
    model = out_struct[0]

    # remove experimental ligands if requested
    if strip_hets:
        strip_experimental_hets(model, keep_water=True)

    # add the Gnina ligand as a fresh chain
    lig_chain = unique_chain_name(model, 'L')
    add_structure_as_chain(lig_struct, model, lig_chain)

    # minimal housekeeping
    out_struct.setup_entities()

    # final sanity
    tot_atoms = sum(len(r) for c in model for r in c)
    if tot_atoms == 0:
        print(f"[ERROR] {os.path.basename(gnina_cif)}: merged structure has 0 atoms")
        return False

    os.makedirs(out_dir, exist_ok=True)
    base = os.path.splitext(os.path.basename(gnina_cif))[0]
    out_path = os.path.join(out_dir, f"{base}_merged.cif")
    with open(out_path, 'w') as fh:
        fh.write(out_struct.make_mmcif_document().as_string())
    if debug:
        print(f"[OK] {os.path.basename(gnina_cif)} -> {os.path.basename(out_path)}  atoms={tot_atoms}")
    return True

def _init_worker():
    signal.signal(signal.SIGINT, signal.SIG_IGN)

def _worker_task(args):
    path, assemblies_dir, out_dir, keep_water, debug = args
    try:
        ok, msg = merge_one(
            path,
            assemblies_dir,
            out_dir,
            keep_water=keep_water,
            debug=debug
        )
        return msg
    except Exception as e:
        return f"[ERROR] {path}: {e}"



def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--src', required=True)
    ap.add_argument('--out', required=True)
    ap.add_argument('--assemblies', required=True)
    ap.add_argument('--glob', default='*.cif')
    ap.add_argument('--no-strip-hets', action='store_true',
                    help="Do not strip experimental HET ligands from receptor")
    ap.add_argument('--debug', action='store_true')
    ap.add_argument("-j", "--workers", type=int, default=1)
    args = ap.parse_args()

    os.makedirs(args.out, exist_ok=True)
    if not paths:
        print(f"[WARN] No input files matched.", file=sys.stderr)
        sys.exit(0)

    tasks = [(path, args.assemblies, args.out, args.keep_water, args.debug)
        for path in paths]

    if args.workers <= 1:
        wrote = errors = 0
        for t in tasks:
            msg = _worker_task(t)
            print(msg)
            if msg.startswith("[OK]"):
                wrote += 1
            elif msg.startswith("[ERROR]"):
                errors += 1
        print(f"Done. Wrote {wrote} complexes; {errors} errors.")
        return

    wrote = errors = 0
    ctx = mp.get_context("spawn")
    with ctx.Pool(processes=args.workers,
                  maxtasksperchild=25,
                  initializer=_init_worker) as pool:
        for msg in pool.imap_unordered(_worker_task, tasks, chunksize=4):
            print(msg)
            if msg.startswith("[OK]"):
                wrote += 1
            elif msg.startswith("[ERROR]"):
                errors += 1
    print(f"Done. Wrote {wrote} complexes; {errors} errors.")

if __name__ == '__main__':
    main()
