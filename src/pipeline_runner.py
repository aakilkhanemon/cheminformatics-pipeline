import argparse
import os
from src.curator_2d import (
    step1_deduplicate,
    step2_sanitize,
    step3_desalt,
    step4_standardize,
    step5_bioavailability_filter,
    step6_pains_filter,
)
from src.expander_chem import step7_physiological_ionization, step8_enumerate_stereoisomers
from src.builder_3d import steps9_to_14_3d_refinement


def main():
    parser = argparse.ArgumentParser(description="MarDockPrep: 14-Step Marine Cheminformatics Preparation Pipeline")
    parser.add_argument("--input", required=True, help="Path to raw input SDF file")
    parser.add_argument("--outdir", default="data/processed", help="Directory for step-wise outputs")
    parser.add_argument("--filter-mode", choices=["ro5", "bro5"], default="bro5", help="Physicochemical filter regime")
    parser.add_argument("--max-isomers", type=int, default=32, help="Max stereoisomers enumerated per parent")
    parser.add_argument("--num-confs", type=int, default=20, help="Initial ETKDGv3 conformers per isomer")
    parser.add_argument("--rmsd-cutoff", type=float, default=0.5, help="Heavy-atom RMSD pruning cutoff (Angstroms)")
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    os.makedirs("results/tables", exist_ok=True)

    s1 = step1_deduplicate(args.input, f"{args.outdir}/step01_dedup.sdf", "results/tables/step01_dedup.csv")
    s2 = step2_sanitize(f"{args.outdir}/step01_dedup.sdf", f"{args.outdir}/step02_sanitized.sdf", "results/tables/step02_sanitized.csv")
    s3 = step3_desalt(f"{args.outdir}/step02_sanitized.sdf", f"{args.outdir}/step03_desalted.sdf", "results/tables/step03_desalted.csv")
    s4 = step4_standardize(f"{args.outdir}/step03_desalted.sdf", f"{args.outdir}/step04_standardized.sdf", "results/tables/step04_standardized.csv")
    s5 = step5_bioavailability_filter(f"{args.outdir}/step04_standardized.sdf", f"{args.outdir}/step05_filtered.sdf", "results/tables/step05_filtered.csv", mode=args.filter_mode)
    s6 = step6_pains_filter(f"{args.outdir}/step05_filtered.sdf", f"{args.outdir}/step06_pains.sdf", "results/tables/step06_pains.csv")
    s7 = step7_physiological_ionization(f"{args.outdir}/step06_pains.sdf", f"{args.outdir}/step07_ionized.sdf", "results/tables/step07_ionized.csv")
    s8 = step8_enumerate_stereoisomers(f"{args.outdir}/step07_ionized.sdf", f"{args.outdir}/step08_isomers.sdf", "results/tables/step08_isomers.csv", max_isomers=args.max_isomers)
    s14 = steps9_to_14_3d_refinement(
        f"{args.outdir}/step08_isomers.sdf",
        f"{args.outdir}/mardockprep_final_docking_ready.sdf",
        "results/tables/step14_final_library.csv",
        num_confs=args.num_confs,
        rmsd_cutoff=args.rmsd_cutoff,
    )

    print(f"Pipeline Complete | Raw -> Final Docking Conformers: {len(s1)} -> {len(s14)}")


if __name__ == "__main__":
    main()
