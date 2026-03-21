"""Restructure data from wide to long format for modeling."""

from __future__ import annotations

import numpy as np
import pandas as pd
from numpy.typing import NDArray


def restructure_for_modeling(
    df: pd.DataFrame,
    norm_i: NDArray,
    sn_mat: NDArray,
    sample_file: pd.DataFrame,
    bridge_mod: bool,
    scale_sn: float,
) -> tuple[pd.DataFrame, list[str]]:
    """Melt intensity/SNR matrices and merge with sample covariates.

    Returns:
        Tuple of (ready_df, sample_names).
    """
    id_vars = df[["PA.Gene.Symbol", "Protein.ID", "Peptide", "Plex"]].reset_index(drop=True)
    n_rows = len(id_vars)

    # Build wide DataFrames with Scan column
    i_cols = [c for c in df.columns if "Adjusted.Intensity" in c]
    sn_cols = [c for c in df.columns if c.endswith(".Sn")]

    norm_df = pd.DataFrame(norm_i, columns=i_cols).reset_index(drop=True)
    sn_df = pd.DataFrame(sn_mat, columns=sn_cols).reset_index(drop=True)

    scan_col = pd.Series(range(1, n_rows + 1), name="Scan")

    # Melt intensities
    i_wide = pd.concat([scan_col, id_vars, norm_df], axis=1)
    melt_i = i_wide.melt(
        id_vars=["Scan", "PA.Gene.Symbol", "Protein.ID", "Peptide", "Plex"],
        var_name="variable",
        value_name="value",
    )

    # Melt SNR
    sn_wide = pd.concat([scan_col, id_vars, sn_df], axis=1)
    melt_sn = sn_wide.melt(
        id_vars=["Scan", "PA.Gene.Symbol", "Protein.ID", "Peptide", "Plex"],
        var_name="Channel",
        value_name="SNR",
    )

    # Combine
    melt_sn["lIntensity"] = np.log2(melt_i["value"].values)
    melt_sn["SampleID"] = melt_sn["Plex"].astype(str) + "_" + melt_sn["Channel"].astype(str)

    # Merge with sample covariates
    cov_names = list(sample_file.columns[1:])
    if cov_names:
        # Match SampleIDs
        valid = melt_sn["SampleID"].isin(sample_file["SampleID"])
        melt_sn = melt_sn[valid].copy()

        id_match = sample_file.set_index("SampleID").loc[melt_sn["SampleID"].values]
        for col in cov_names:
            melt_sn[col] = id_match[col].values

    ready_df = melt_sn.copy()

    # Remove NaN lIntensity
    ready_df = ready_df[~ready_df["lIntensity"].isna()].copy()

    # Bridge or scan-mean normalization
    if not bridge_mod:
        scan_mean = ready_df.groupby("Scan")["lIntensity"].transform("mean")
        ready_df["lIntensity"] = ready_df["lIntensity"] - scan_mean
    else:
        bridge_ids = sample_file.loc[sample_file["Bridge"] == 1, "SampleID"].values
        ready_df["Bridge"] = (ready_df["SampleID"].isin(bridge_ids)).astype(int)

    # Sample names and weights
    sample_names = sorted(ready_df["SampleID"].unique().tolist())
    ready_df["Plex"] = pd.Categorical(ready_df["Plex"])
    ready_df["techVar"] = 1.0 / (scale_sn * ready_df["SNR"])

    return ready_df, sample_names
