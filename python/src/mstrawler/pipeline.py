"""Main msTrawler pipeline orchestrator."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from .covariate_setup import CovariateConfig, setup_covariates
from .file_converter import convert_pd_file
from .model_fitting import adapt_model, check_backend
from .preprocessing import find_outliers, geo_norm, tmt_lod
from .restructure import restructure_for_modeling
from .results import init_result_tables, write_result_tables

# Markers for outlier/overflow filtering
_MARKER_OUTLIER = "OUTLIER_REMOVE_AT_ONCE!"
_MARKER_TOO_MANY = "TOO_MANY!!!"


def ms_trawl(
    df: pd.DataFrame,
    sample_file: pd.DataFrame | None = None,
    covariate_file: pd.DataFrame | None = None,
    scale_sn: float = 1.0,
    lod: float = 0.01,
    impute_penalty: float = 1.0,
    min_above: int = 4,
    ssn_filter: float = 20.0,
    outlier_cutoff: float = 3.0,
    n_sum: int = 3,
    swap_protein: bool = False,
    max_pep: int = 25,
    col_adjust: float | np.ndarray | None = 0.5,
    col_ratios: np.ndarray | None = None,
    drop_contam: bool = True,
    drop_reverse: bool = True,
    peptide_analysis: bool = False,
    min_re: int = 5,
    time_diff: bool = True,
    seed: int | None = None,
    backend: str = "python",
) -> None:
    """Run the full msTrawler analysis pipeline.

    Args:
        backend: Model fitting backend. Options:
            - "python" (default): pure Python using statsmodels WLS + Wald F-tests.
              No R dependency. Best for production/deployment.
            - "r": uses R's lme4 + pbkrtest via rpy2 for mixed-effects models
              and Kenward-Roger tests. Requires R + rpy2 + R packages.
              Best for small-sample experiments and publication-quality statistics.

    Outputs CSV files to the current working directory.
    """
    # Validate backend
    status = check_backend(backend)
    if not status["available"]:
        raise RuntimeError(
            f"Backend '{backend}' is not available: {status['details']}"
        )
    rng = np.random.default_rng(seed)

    # Validate covariates
    if covariate_file is not None:
        _validate_covariates(covariate_file, sample_file)

    # Clean input data
    df, sn_mat, i_mat = _clean_input_data(
        df, col_adjust, drop_reverse, drop_contam, peptide_analysis, ssn_filter
    )

    # Preprocessing pipeline (LOD + outliers + normalization)
    df, sn_mat, norm_i, u_plex = _preprocess_pipeline(
        df, sn_mat, i_mat, lod, min_above, scale_sn, impute_penalty,
        outlier_cutoff, n_sum, swap_protein, max_pep, peptide_analysis,
        col_adjust, col_ratios, rng,
    )

    # Set up covariates
    cov = setup_covariates(sample_file, covariate_file, u_plex)

    # Restructure for modeling
    ready_df, sample_names = restructure_for_modeling(
        df, norm_i, sn_mat, cov.sample_file, cov.bridge_mod, scale_sn,
    )

    # Initialize result tables
    tab_res = init_result_tables(
        ready_df, sample_names, cov.n_cat, cov.factor_names,
        cov.level_names, cov.n_levels, cov.n_cont, cov.cont_names,
        cov.time_level_n, cov.time_vars, cov.t_cat_levels,
    )
    table_list = tab_res["table_list"]
    time_tables = tab_res["time_tables"]
    u_prot = tab_res["u_prot"]

    # Run per-protein modeling loop
    table_list, time_tables = _run_protein_models(
        ready_df, u_prot, table_list, time_tables,
        cov, sample_names, scale_sn, min_re, time_diff,
        backend=backend,
    )

    # Write results with FDR
    write_result_tables(
        table_list, cov.table_names, time_tables, cov.time_table_names,
        apply_fdr=True, suffix=None,
    )


def _validate_covariates(covariate_file, sample_file):
    """Validate covariate and sample file consistency."""
    co_vector = covariate_file["Covariate"].values
    samp_vector = sample_file.columns
    n_matches = sum(c in samp_vector for c in co_vector)
    if n_matches < len(co_vector):
        raise ValueError("At least one covariate name does not match across files.")

    # Check for missing covariates (ignoring bridge)
    bridge_cols = [c for c in sample_file.columns if c.upper() == "BRIDGE"]
    if bridge_cols:
        non_bridge = sample_file[sample_file[bridge_cols[0]] != 1]
        if non_bridge.isna().sum().sum() > 0:
            raise ValueError("Missing values are not allowed in covariates")
    elif sample_file.isna().sum().sum() > 0:
        raise ValueError("Missing values are not allowed in covariates")


def _clean_input_data(df, col_adjust, drop_reverse, drop_contam, peptide_analysis, ssn_filter):
    """Filter, order, and extract matrices."""
    if isinstance(col_adjust, np.ndarray) and len(col_adjust) > 1:
        df = df.copy()
        df["colAdjust"] = col_adjust

    if drop_reverse:
        df = df[~df["Protein.ID"].astype(str).str.contains("##", na=False)]
    if drop_contam:
        df = df[~df["Protein.ID"].astype(str).str.contains("contam", na=False)]

    df = df.sort_values(["Protein.ID", "Plex", "Peptide"]).reset_index(drop=True)

    if peptide_analysis:
        df["PA.Gene.Symbol"] = df["PA.Gene.Symbol"].astype(str) + "___" + df["Protein.ID"].astype(str)
        df["Protein.ID"] = df["Peptide"]

    sn_cols = [c for c in df.columns if c.endswith(".Sn")]
    i_cols = [c for c in df.columns if "Adjusted.Intensity" in c]
    sn_mat = df[sn_cols].values.astype(float)
    i_mat = df[i_cols].values.astype(float)

    if sn_mat.shape[1] != i_mat.shape[1]:
        raise ValueError("Number of SNR columns doesn't match intensity columns.")

    # SSN filter
    if ssn_filter is not None:
        ssn = np.sum(sn_mat, axis=1)
        keep = ssn >= ssn_filter
        df = df[keep].reset_index(drop=True)
        sn_mat = sn_mat[keep]
        i_mat = i_mat[keep]

    # Remove rows with < 2 non-zero intensities
    n_nonzero = np.sum(i_mat > 0, axis=1)
    keep = n_nonzero >= 2
    df = df[keep].reset_index(drop=True)
    sn_mat = sn_mat[keep]
    i_mat = i_mat[keep]

    return df, sn_mat, i_mat


def _preprocess_pipeline(
    df, sn_mat, i_mat, lod, min_above, scale_sn, impute_penalty,
    outlier_cutoff, n_sum, swap_protein, max_pep, peptide_analysis,
    col_adjust, col_ratios, rng,
):
    """LOD imputation, outlier detection, and normalization."""
    # LOD
    new_i, new_sn, too_few = tmt_lod(sn_mat, i_mat, lod, min_above, scale_sn, impute_penalty, rng)
    if too_few.sum() > 0:
        keep = ~too_few
        df = df[keep].reset_index(drop=True)
        new_i = new_i[keep]
        new_sn = new_sn[keep]

    i_mat = new_i
    sn_mat = new_sn

    # Outlier detection
    df = df.copy()
    df["Protein.ID"] = df["Protein.ID"].astype(str)
    df["Peptide"] = df["Peptide"].astype(str)
    u_plex = df["Plex"].unique()
    u_prot = df["Protein.ID"].unique()

    remove_mask = np.zeros(len(df), dtype=bool)

    for plex in u_plex:
        for prot in u_prot:
            idx = np.where((df["Protein.ID"].values == prot) & (df["Plex"].values == plex))[0]
            if len(idx) == 0:
                continue

            sub_i = i_mat[idx]
            sub_sn = sn_mat[idx]

            if outlier_cutoff is not None and len(idx) >= 2:
                outlier_mat = find_outliers(sub_i, sub_sn, outlier_cutoff, scale_sn)
                outlier_rows = np.sum(outlier_mat, axis=1) > 0
            else:
                outlier_rows = np.zeros(len(idx), dtype=bool)

            outlier_idx = idx[outlier_rows]
            if len(outlier_idx) > 0:
                if swap_protein:
                    for oi in outlier_idx:
                        df.iloc[oi, df.columns.get_loc("Protein.ID")] = (
                            df.iloc[oi]["Protein.ID"] + "___" + df.iloc[oi]["Peptide"]
                        )
                else:
                    remove_mask[outlier_idx] = True

            # Aggregation check
            n_after = int((~outlier_rows).sum())
            if (n_after < n_sum and n_after > 1) or (peptide_analysis and n_after > 1):
                good_idx = idx[~outlier_rows]
                keep_idx = good_idx[0]
                sn_cols = [c for c in df.columns if c.endswith(".Sn")]
                i_cols = [c for c in df.columns if "Adjusted.Intensity" in c]
                df.iloc[keep_idx, df.columns.get_loc("Peptide")] = df.iloc[idx[0]]["Peptide"] + "_SUM"
                for k, sc in enumerate(sn_cols):
                    df.iat[keep_idx, df.columns.get_loc(sc)] = sub_sn[~outlier_rows, k].sum()
                for k, ic in enumerate(i_cols):
                    df.iat[keep_idx, df.columns.get_loc(ic)] = sub_i[~outlier_rows, k].sum()
                sn_mat[keep_idx] = sub_sn[~outlier_rows].sum(axis=0)
                i_mat[keep_idx] = sub_i[~outlier_rows].sum(axis=0)
                for ri in idx:
                    if ri != keep_idx:
                        remove_mask[ri] = True
                continue

            # Max pep limit
            if len(idx) - len(outlier_idx) > max_pep:
                sub_ssn = np.sum(sub_sn, axis=1)
                good_mask = ~outlier_rows
                good_indices = np.where(good_mask)[0]
                sorted_good = good_indices[np.argsort(-sub_ssn[good_indices])]
                to_remove = sorted_good[max_pep:]
                for ri in to_remove:
                    remove_mask[idx[ri]] = True

    if remove_mask.sum() > 0:
        keep = ~remove_mask
        df = df[keep].reset_index(drop=True)
        i_mat = i_mat[keep]
        sn_mat = sn_mat[keep]

    # Normalization
    if col_adjust is None:
        norm_i = i_mat
    else:
        if isinstance(col_adjust, (int, float)) and not isinstance(col_adjust, np.ndarray):
            norm_bool = np.zeros(len(df))
            for p in u_plex:
                plex_mask = df["Plex"].values == p
                plex_idx = np.where(plex_mask)[0]
                if len(plex_idx) <= 1:
                    raise ValueError("Plex with <= 1 observation")
                pep_sd = np.array([
                    np.nanstd(np.log(i_mat[pi])) for pi in plex_idx
                ])
                med_sd = np.nanquantile(pep_sd, col_adjust)
                norm_bool[plex_idx[pep_sd < med_sd]] = 1
        else:
            norm_bool = df["colAdjust"].values if "colAdjust" in df.columns else np.ones(len(df))

        n_plex = len(np.unique(df["Plex"].values))
        if col_ratios is None or len(col_ratios) < 2:
            ratios = np.ones(n_plex * i_mat.shape[1])
        else:
            ratios = col_ratios

        norm_i, norm_facs = geo_norm(i_mat, norm_bool, df["Plex"].values, ratios)
        pd.DataFrame(norm_facs).to_csv("ColAdjustmentFactors.csv", index=False)

    return df, sn_mat, norm_i, u_plex


def _run_protein_models(
    ready_df, u_prot, table_list, time_tables,
    cov: CovariateConfig, sample_names, scale_sn, min_re, time_diff,
    backend: str = "python",
):
    """Per-protein modeling loop."""
    for prot_idx, prot_id in enumerate(u_prot):
        prot_dat = ready_df[ready_df["Protein.ID"] == prot_id].copy()

        # Split bridge/non-bridge
        if cov.bridge_mod:
            bridge_dat = prot_dat[prot_dat["Bridge"] == 1]
            not_bridge = prot_dat[prot_dat["Bridge"] == 0].copy()
        else:
            bridge_dat = None
            not_bridge = prot_dat.copy()

        # Sample estimates (covType = "None")
        res = adapt_model(
            not_bridge, "None", cov.fixed_str, bridge_dat, cov.t_parm,
            None, sample_names, min_re, cov.rand_id, False, False,
            backend=backend,
        )
        if res is not None:
            estimates, _, model_type = res
            n_data_cols = len(table_list[0].columns) - 3
            if len(estimates) >= n_data_cols:
                table_list[0].iloc[prot_idx, 3:] = estimates[:n_data_cols]
            table_list[0].iloc[prot_idx, 2] = model_type

        # Factor models
        tab_index = 1
        if cov.n_cat > 0 and cov.level_names is not None:
            for i, fn in enumerate(cov.factor_names):
                if fn == "":
                    continue
                levels = cov.level_names[i]
                for j in range(len(levels)):
                    # Set reference level
                    non_ref = [lv for k, lv in enumerate(levels) if k != j]
                    full_cols = [f"{fn}{lv}" for lv in non_ref]

                    res = adapt_model(
                        not_bridge, "Factor", cov.fixed_str, bridge_dat, cov.t_parm,
                        None, full_cols, min_re, cov.rand_id, False, True,
                        backend=backend,
                    )
                    if res is not None and tab_index < len(table_list):
                        estimates, time_res, model_type = res
                        n_data_cols = len(table_list[tab_index].columns) - 3
                        if len(estimates) >= n_data_cols:
                            table_list[tab_index].iloc[prot_idx, 3:] = estimates[:n_data_cols]
                        table_list[tab_index].iloc[prot_idx, 2] = model_type

                    tab_index += 1

        # Continuous models
        if cov.n_cont > 0:
            cont_cols = [cn for cn in cov.cont_names if cn != ""]
            res = adapt_model(
                not_bridge, "Continuous", cov.fixed_str, bridge_dat, cov.t_parm,
                None, cont_cols, min_re, cov.rand_id, False, True,
                backend=backend,
            )
            if res is not None:
                estimates, _, model_type = res
                for k, cn in enumerate(cont_cols):
                    if tab_index < len(table_list):
                        if isinstance(estimates, np.ndarray) and estimates.ndim == 2 and k < len(estimates):
                            table_list[tab_index].iloc[prot_idx, 3:6] = estimates[k]
                        table_list[tab_index].iloc[prot_idx, 2] = model_type
                    tab_index += 1

        # Time models
        if cov.t_parm == "Time" and time_tables:
            lht_list = [["X_" + tv for tv in cov.time_vars]] if cov.time_vars else None
            res = adapt_model(
                not_bridge, "Time", cov.fixed_str, bridge_dat, cov.t_parm,
                lht_list, cov.time_vars or [], min_re, cov.rand_id, False, True,
                backend=backend,
            )
            if res is not None:
                _, time_res, model_type = res
                if time_res is not None:
                    n_data = len(time_tables[0].columns) - 3
                    time_tables[0].iloc[prot_idx, 3:3 + min(len(time_res), n_data)] = time_res[:n_data]
                time_tables[0].iloc[prot_idx, 2] = model_type

        print(f"Protein {prot_idx + 1} Complete")

    return table_list, time_tables
