"""E2E tests for preprocessing functions against R golden baseline.

Note: LOD imputation uses random draws, so R and Python produce different imputed
values even with the same seed (different RNG algorithms). We test:
1. LOD: same censoring pattern (tooFew flags match)
2. LOD: same dimensions and non-imputed values match
3. geoNorm: exact match (deterministic)
"""

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from mstrawler.file_converter import convert_pd_file
from mstrawler.preprocessing import tmt_lod, geo_norm

REPO_ROOT = Path(__file__).resolve().parents[3]
TUTORIAL_DIR = REPO_ROOT / "data" / "tutorial_example"
GOLDEN_DIR = REPO_ROOT / "tests" / "golden"

TOLERANCE = 1e-6


@pytest.fixture
def converted_df(tmp_path):
    import shutil

    psm = tmp_path / "PSMs.txt"
    prot = tmp_path / "Proteins.txt"
    shutil.copy(TUTORIAL_DIR / "pd_example_export_PSMs.txt", psm)
    shutil.copy(TUTORIAL_DIR / "pd_example_export_Proteins.txt", prot)
    return convert_pd_file(str(psm), str(prot))


@pytest.fixture
def matrices(converted_df):
    sn_cols = [c for c in converted_df.columns if c.endswith(".Sn")]
    i_cols = [c for c in converted_df.columns if "Adjusted.Intensity" in c]
    sn_mat = converted_df[sn_cols].to_numpy(dtype=float)
    i_mat = converted_df[i_cols].to_numpy(dtype=float)
    return sn_mat, i_mat, converted_df


# ---------- LOD Tests ----------


class TestTmtLOD:
    def test_too_few_flags_match(self, matrices):
        """The tooFew boolean vector should match R (deterministic, no randomness)."""
        sn_mat, i_mat, _ = matrices
        golden_toofew = pd.read_csv(GOLDEN_DIR / "lod_toofew.csv")

        _, _, too_few = tmt_lod(sn_mat, i_mat, lod=0.01, min_above=4, scale_sn=1, impute_penalty=1)

        r_flags = golden_toofew["tooFew"].to_numpy().astype(bool)
        assert np.array_equal(too_few, r_flags), (
            f"tooFew mismatch: Python has {too_few.sum()} flagged, R has {r_flags.sum()}"
        )

    def test_dimensions_match(self, matrices):
        sn_mat, i_mat, _ = matrices
        golden_i = pd.read_csv(GOLDEN_DIR / "lod_intensities.csv")
        golden_sn = pd.read_csv(GOLDEN_DIR / "lod_snr.csv")

        new_i, new_sn, _ = tmt_lod(
            sn_mat, i_mat, lod=0.01, min_above=4, scale_sn=1, impute_penalty=1
        )

        assert new_i.shape == golden_i.shape, f"Intensity shape: {new_i.shape} vs {golden_i.shape}"
        assert new_sn.shape == golden_sn.shape, f"SNR shape: {new_sn.shape} vs {golden_sn.shape}"

    def test_non_censored_values_unchanged(self, matrices):
        """Values above LOD should not be modified by imputation."""
        sn_mat, i_mat, _ = matrices

        # Determine which entries are NOT censored (same logic as R)
        ssi = np.nansum(i_mat, axis=1)
        pi_mat = i_mat / ssi[:, np.newaxis]
        not_censored = (pi_mat >= 0.01) & ~np.isnan(pi_mat)

        new_i, new_sn, _ = tmt_lod(
            sn_mat, i_mat, lod=0.01, min_above=4, scale_sn=1, impute_penalty=1
        )

        # Non-censored intensities should be identical
        np.testing.assert_allclose(
            new_i[not_censored], i_mat[not_censored], rtol=TOLERANCE,
            err_msg="Non-censored intensities were modified!"
        )

        # Non-censored SNR should be identical
        np.testing.assert_allclose(
            new_sn[not_censored], sn_mat[not_censored], rtol=TOLERANCE,
            err_msg="Non-censored SNR values were modified!"
        )


# ---------- Normalization Tests ----------


class TestGeoNorm:
    def test_normalization_factors_match(self, matrices):
        """Normalization factors should match R exactly (deterministic)."""
        sn_mat, i_mat, df = matrices
        golden_factors = pd.read_csv(GOLDEN_DIR / "geonorm_factors.csv")

        # Replicate R preprocessing: LOD → remove tooFew → then normalize
        # Use R golden LOD output to avoid RNG divergence
        golden_i = pd.read_csv(GOLDEN_DIR / "lod_intensities.csv").to_numpy(dtype=float)
        golden_toofew = pd.read_csv(GOLDEN_DIR / "lod_toofew.csv")["tooFew"].to_numpy().astype(bool)

        clean_i = golden_i[~golden_toofew]
        clean_plex = df["Plex"].to_numpy()[~golden_toofew]

        norm_bool = np.ones(clean_i.shape[0])
        n_plex = len(np.unique(clean_plex))
        ratios = np.ones(n_plex * clean_i.shape[1])

        norm_i, norm_factors = geo_norm(clean_i, norm_bool, clean_plex, ratios)

        # Compare factors
        r_factors = golden_factors.to_numpy(dtype=float)
        np.testing.assert_allclose(
            norm_factors, r_factors, rtol=TOLERANCE,
            err_msg="Normalization factors don't match R"
        )

    def test_normalized_intensities_match(self, matrices):
        """Normalized intensity matrix should match R exactly."""
        sn_mat, i_mat, df = matrices
        golden_norm = pd.read_csv(GOLDEN_DIR / "geonorm_intensities.csv")

        golden_i = pd.read_csv(GOLDEN_DIR / "lod_intensities.csv").to_numpy(dtype=float)
        golden_toofew = pd.read_csv(GOLDEN_DIR / "lod_toofew.csv")["tooFew"].to_numpy().astype(bool)

        clean_i = golden_i[~golden_toofew]
        clean_plex = df["Plex"].to_numpy()[~golden_toofew]

        norm_bool = np.ones(clean_i.shape[0])
        n_plex = len(np.unique(clean_plex))
        ratios = np.ones(n_plex * clean_i.shape[1])

        norm_i, _ = geo_norm(clean_i, norm_bool, clean_plex, ratios)

        r_norm = golden_norm.to_numpy(dtype=float)
        np.testing.assert_allclose(
            norm_i, r_norm, rtol=TOLERANCE,
            err_msg="Normalized intensities don't match R"
        )
