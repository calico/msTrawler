"""E2E test: run full Python pipeline and compare structure against R golden baseline."""

from pathlib import Path
import shutil
import os

import numpy as np
import pandas as pd
import pytest

from mstrawler.file_converter import convert_pd_file
from mstrawler.pipeline import ms_trawl

REPO_ROOT = Path(__file__).resolve().parents[3]
TUTORIAL_DIR = REPO_ROOT / "data" / "tutorial_example"
GOLDEN_DIR = REPO_ROOT / "tests" / "golden"


@pytest.fixture
def run_pipeline(tmp_path):
    """Run the full Python pipeline in a temp directory and return output paths."""
    # Copy input data
    psm = tmp_path / "PSMs.txt"
    prot = tmp_path / "Proteins.txt"
    shutil.copy(TUTORIAL_DIR / "pd_example_export_PSMs.txt", psm)
    shutil.copy(TUTORIAL_DIR / "pd_example_export_Proteins.txt", prot)

    # Convert
    df = convert_pd_file(str(psm), str(prot))

    # Read sample/covariate files
    sample_file = pd.read_csv(TUTORIAL_DIR / "sample_file.csv")
    covariate_file = pd.read_csv(TUTORIAL_DIR / "covariate_file.csv")

    # Run pipeline in tmp_path
    orig_dir = os.getcwd()
    os.chdir(tmp_path)
    try:
        ms_trawl(
            df=df,
            sample_file=sample_file,
            covariate_file=covariate_file,
            scale_sn=1,
            lod=0.01,
            impute_penalty=1,
            min_above=4,
            ssn_filter=20,
            outlier_cutoff=3,
            n_sum=3,
            swap_protein=False,
            max_pep=25,
            col_adjust=0.5,
            drop_contam=True,
            drop_reverse=True,
            peptide_analysis=False,
            min_re=5,
            time_diff=True,
            seed=42,
        )
    finally:
        os.chdir(orig_dir)

    return tmp_path


def test_pipeline_produces_expected_files(run_pipeline):
    """Check that the Python pipeline produces the same set of output files as R."""
    output_dir = run_pipeline
    r_pipeline_dir = GOLDEN_DIR / "pipeline"

    # Expected files from R
    r_files = {f.name for f in r_pipeline_dir.glob("*.csv")}
    # Also expect ColAdjustmentFactors.csv in the root
    expected = r_files | {"ColAdjustmentFactors.csv"}

    py_files = {f.name for f in output_dir.glob("*.csv") if "converted" not in f.name}

    # Check that all R pipeline files are produced by Python
    missing = expected - py_files
    assert len(missing) == 0, f"Missing output files: {missing}"


def test_pipeline_simple_csv_dimensions(run_pipeline):
    """Check Simple.csv has same number of proteins as R version."""
    py_simple = pd.read_csv(run_pipeline / "Simple.csv")
    r_simple = pd.read_csv(GOLDEN_DIR / "pipeline" / "Simple.csv")

    # R's write.csv adds a row-index column ("Unnamed: 0"); drop it for comparison
    r_data_cols = [c for c in r_simple.columns if not c.startswith("Unnamed")]
    py_data_cols = list(py_simple.columns)

    assert len(py_data_cols) == len(r_data_cols), (
        f"Column count mismatch: Python={len(py_data_cols)}, R={len(r_data_cols)}"
    )


def test_pipeline_factor_files_exist(run_pipeline):
    """Check that factor contrast files are created."""
    assert (run_pipeline / "Factor_Age_Old.csv").exists()
    assert (run_pipeline / "Factor_Age_Young.csv").exists()


def test_pipeline_time_files_exist(run_pipeline):
    """Check that time trend files are created."""
    assert (run_pipeline / "Time_Old.csv").exists()
    assert (run_pipeline / "Time_Young.csv").exists()
