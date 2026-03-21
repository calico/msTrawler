"""E2E test: compare Python file converter output against R golden baseline."""

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from mstrawler.file_converter import convert_pd_file

REPO_ROOT = Path(__file__).resolve().parents[3]
TUTORIAL_DIR = REPO_ROOT / "data" / "tutorial_example"
GOLDEN_DIR = REPO_ROOT / "tests" / "golden"

TOLERANCE = 1e-6


@pytest.fixture
def golden_converted_df():
    return pd.read_csv(GOLDEN_DIR / "converted_df.csv")


@pytest.fixture
def python_converted_df(tmp_path):
    # Run converter with tutorial data, using tmp_path to avoid polluting repo
    import shutil

    psm_file = tmp_path / "pd_example_export_PSMs.txt"
    protein_file = tmp_path / "pd_example_export_Proteins.txt"
    shutil.copy(TUTORIAL_DIR / "pd_example_export_PSMs.txt", psm_file)
    shutil.copy(TUTORIAL_DIR / "pd_example_export_Proteins.txt", protein_file)

    return convert_pd_file(str(psm_file), str(protein_file))


def test_dimensions_match(golden_converted_df, python_converted_df):
    assert golden_converted_df.shape == python_converted_df.shape, (
        f"Shape mismatch: R={golden_converted_df.shape}, Python={python_converted_df.shape}"
    )


def test_columns_match(golden_converted_df, python_converted_df):
    assert list(golden_converted_df.columns) == list(python_converted_df.columns)


def test_numeric_values_match(golden_converted_df, python_converted_df):
    for col in golden_converted_df.columns:
        g = golden_converted_df[col]
        p = python_converted_df[col]

        if pd.api.types.is_numeric_dtype(g):
            # Compare non-NA values
            g_vals = g.to_numpy(dtype=float)
            p_vals = p.to_numpy(dtype=float)

            both_valid = ~(np.isnan(g_vals) | np.isnan(p_vals))
            if both_valid.sum() == 0:
                continue

            max_diff = np.max(np.abs(g_vals[both_valid] - p_vals[both_valid]))
            denom = np.maximum(np.abs(g_vals[both_valid]), np.abs(p_vals[both_valid]))
            denom = np.maximum(denom, 1e-10)
            max_rel_diff = np.max(np.abs(g_vals[both_valid] - p_vals[both_valid]) / denom)

            assert max_diff <= TOLERANCE or max_rel_diff <= TOLERANCE, (
                f"Column '{col}': max_abs_diff={max_diff:.2e}, max_rel_diff={max_rel_diff:.2e}"
            )


def test_string_values_match(golden_converted_df, python_converted_df):
    for col in ["Protein.ID", "Peptide", "Plex"]:
        g = golden_converted_df[col].astype(str).tolist()
        p = python_converted_df[col].astype(str).tolist()
        assert g == p, f"Column '{col}' string mismatch at first diff"
