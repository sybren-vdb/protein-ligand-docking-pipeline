import argparse, subprocess, sys
from pathlib import Path

def main():
    ap = argparse.ArgumentParser("Thin wrapper around gnina (optionally via apptainer)")
    ap.add_argument("-r","--receptor_pdbqt", required=True, type=Path)
    ap.add_argument("-l","--ligands_sdf", required=True, type=Path)
    ap.add_argument("--autobox_ref", required=True, type=Path,
                    help="Reference ligand (PDB/MOL2/SDF accepted by gnina)")
    ap.add_argument("--out_poses", required=True, type=Path)
    ap.add_argument("--out_log", required=True, type=Path)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--exhaustiveness", type=int, default=8)
    ap.add_argument("--autobox_add", type=float, default=4.0)
    ap.add_argument("--sif", type=Path, help="Apptainer SIF to run gnina inside")
    ap.add_argument("--nv", action="store_true", help="Enable GPU with --nv for Apptainer")
    args = ap.parse_args()

    base_cmd = ["gnina",
        "-r", str(args.receptor_pdbqt),
        "-l", str(args.ligands_sdf),
        "--autobox_ligand", str(args.autobox_ref),
        "--autobox_add", str(args.autobox_add),
        "--exhaustiveness", str(args.exhaustiveness),
        "--seed", str(args.seed),
        "--out", str(args.out_poses),
        "--log", str(args.out_log)
    ]

    if args.sif:
        cmd = ["apptainer", "exec"] + (["--nv"] if args.nv else []) + [str(args.sif)] + base_cmd
    else:
        cmd = base_cmd

    args.out_log.parent.mkdir(parents=True, exist_ok=True)
    with open(args.out_log, "a") as lf:
        proc = subprocess.run(cmd, stdout=lf, stderr=subprocess.STDOUT, text=True)
        if proc.returncode != 0:
            sys.exit(proc.returncode)

if __name__ == "__main__":
    main()
