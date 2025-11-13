#!/usr/bin/env python
import argparse, math, json
from pathlib import Path
from collections import defaultdict
from rdkit import Chem

def pK_to_dG_kcal(pK: float, T: float = 298.15) -> float:
    R = 1.98720425864083e-3  # kcal/(mol·K)
    return -2.303 * R * T * pK  # ≈ -1.364 * pK at 298 K

def parse_one_sdf(sdf_path: Path):
    suppl = Chem.SDMolSupplier(str(sdf_path), removeHs=False)
    rows = []
    for i, mol in enumerate(suppl):
        if mol is None: 
            continue
        props = mol.GetPropsAsDict()
        # Robust parse; tolerate strings
        def fget(name):
            v = props.get(name, None)
            if v is None: return None
            try: return float(v)
            except: return None

        cnn_aff_pK = fget("CNNaffinity")
        cnn_score  = fget("CNNscore")
        vina_aff   = fget("minimizedAffinity")  # optional tie-breaker (kcal/mol)

        # heuristics to form a "record id" key from filename; customize if needed
        # e.g., 8ABC_D_lig.sdf  -> record_id = "8ABC_D"
        stem = sdf_path.stem
        record_id = stem.split(".")[0]  # edit if your naming differs

        rows.append({
            "record_id": record_id,
            "sdf": str(sdf_path),
            "pose_idx": i,
            "cnn_affinity_pK": cnn_aff_pK,
            "cnn_score": cnn_score,
            "vina": vina_aff,
        })
    return rows

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("sdf_glob", help="Glob like /path/to/gnina_out/*.sdf")
    ap.add_argument("--topk", type=int, default=3, help="poses per record")
    ap.add_argument("--unit", choices=["pK","dG"], default="pK", help="affinity unit to export")
    ap.add_argument("--out_jsonl", default="affinity_map.jsonl", help="output JSONL")
    args = ap.parse_args()

    paths = sorted(Path().glob(args.sdf_glob)) if "*" in args.sdf_glob else [Path(args.sdf_glob)]
    buckets = defaultdict(list)
    for p in paths:
        for row in parse_one_sdf(p):
            buckets[row["record_id"]].append(row)

    out = []
    for rid, rows in buckets.items():
        # sort by CNNscore desc; tie-break by better (more negative) vina if present
        rows = [r for r in rows if r["cnn_score"] is not None and r["cnn_affinity_pK"] is not None]
        if not rows:
            continue
        rows.sort(key=lambda r: (-r["cnn_score"], r["vina"] if r["vina"] is not None else 1e9))
        top = rows[:args.topk]
        for rank, r in enumerate(top, 1):
            aff = r["cnn_affinity_pK"] if args.unit == "pK" else pK_to_dG_kcal(r["cnn_affinity_pK"])
            out.append({
                "id": rid,                 # must match boltz2 record.id
                "affinity": aff,           # value you’ll inject into record.affinity
                "unit": args.unit,         # pK or dG
                "pose_weight": r["cnn_score"],  # optional: use as sample weight
                "sdf": r["sdf"],
                "pose_idx": r["pose_idx"],
                "pose_rank": rank
            })

    with open(args.out_jsonl, "w") as f:
        for row in out:
            f.write(json.dumps(row) + "\n")

    print(f"Wrote {len(out)} entries to {args.out_jsonl}")

if __name__ == "__main__":
    main()
