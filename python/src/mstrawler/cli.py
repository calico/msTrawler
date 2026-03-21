"""msTrawler command-line interface."""

from __future__ import annotations

import argparse
import os
import sys

import pandas as pd


def main():
    parser = argparse.ArgumentParser(
        prog="mstrawler",
        description="msTrawler: Multiplexed proteomics analysis pipeline",
    )
    subparsers = parser.add_subparsers(dest="command", help="Available commands")

    # --- run command ---
    run_parser = subparsers.add_parser("run", help="Run the full analysis pipeline")
    run_parser.add_argument("--psm", required=True, help="Path to Proteome Discoverer PSM file")
    run_parser.add_argument("--protein", default="", help="Path to Proteome Discoverer Protein file (optional)")
    run_parser.add_argument("--sample", required=True, help="Path to sample file (CSV)")
    run_parser.add_argument("--covariate", required=True, help="Path to covariate file (CSV)")
    run_parser.add_argument("--output", "-o", default=".", help="Output directory (default: current directory)")
    run_parser.add_argument("--backend", choices=["python", "r"], default="python", help="Model fitting backend (default: python)")
    run_parser.add_argument("--seed", type=int, default=None, help="Random seed for reproducibility")
    run_parser.add_argument("--lod", type=float, default=0.01, help="Limit of detection threshold (default: 0.01)")
    run_parser.add_argument("--min-above", type=int, default=4, help="Min channels above LOD per scan (default: 4)")
    run_parser.add_argument("--col-adjust", type=float, default=0.5, help="SD percentile for normalization (default: 0.5)")
    run_parser.add_argument("--outlier-cutoff", type=float, default=3.0, help="Outlier residual threshold (default: 3.0)")
    run_parser.add_argument("--scale-sn", type=float, default=1.0, help="SNR scaling factor (default: 1.0)")
    run_parser.add_argument("--delimiter", default="\t", help="PSM file delimiter (default: tab)")

    # --- convert command ---
    convert_parser = subparsers.add_parser("convert", help="Convert Proteome Discoverer files only")
    convert_parser.add_argument("--psm", required=True, help="Path to PSM file")
    convert_parser.add_argument("--protein", default="", help="Path to Protein file (optional)")
    convert_parser.add_argument("--delimiter", default="\t", help="File delimiter (default: tab)")

    # --- check command ---
    check_parser = subparsers.add_parser("check", help="Check backend availability")

    args = parser.parse_args()

    if args.command is None:
        parser.print_help()
        sys.exit(1)

    if args.command == "run":
        _cmd_run(args)
    elif args.command == "convert":
        _cmd_convert(args)
    elif args.command == "check":
        _cmd_check()


def _cmd_run(args):
    from .file_converter import convert_pd_file
    from .pipeline import ms_trawl

    print(f"=== msTrawler Analysis Pipeline (backend={args.backend}) ===")
    print()

    # Convert PSM file
    print(f"Converting: {args.psm}")
    df = convert_pd_file(args.psm, args.protein, args.delimiter)
    print(f"  {len(df)} scans, {df['Protein.ID'].nunique()} proteins")
    print()

    # Load metadata
    print(f"Sample file: {args.sample}")
    sample_file = pd.read_csv(args.sample)
    print(f"Covariate file: {args.covariate}")
    covariate_file = pd.read_csv(args.covariate)
    print()

    # Run pipeline
    os.makedirs(args.output, exist_ok=True)
    orig_dir = os.getcwd()
    os.chdir(args.output)
    try:
        ms_trawl(
            df=df,
            sample_file=sample_file,
            covariate_file=covariate_file,
            scale_sn=args.scale_sn,
            lod=args.lod,
            min_above=args.min_above,
            col_adjust=args.col_adjust,
            outlier_cutoff=args.outlier_cutoff,
            seed=args.seed,
            backend=args.backend,
        )
    finally:
        os.chdir(orig_dir)

    # List outputs
    print()
    print(f"=== Output written to: {args.output} ===")
    for f in sorted(os.listdir(args.output)):
        if f.endswith(".csv"):
            print(f"  {f}")


def _cmd_convert(args):
    from .file_converter import convert_pd_file

    print(f"Converting: {args.psm}")
    df = convert_pd_file(args.psm, args.protein, args.delimiter)
    print(f"  {len(df)} scans, {df['Protein.ID'].nunique()} proteins")
    print(f"  Output: {args.psm}_converted.csv")


def _cmd_check():
    from .model_fitting import check_backend

    for backend in ["python", "r"]:
        status = check_backend(backend)
        icon = "OK" if status["available"] else "NOT AVAILABLE"
        print(f"  [{icon}] {backend}: {status['details']}")


if __name__ == "__main__":
    main()
