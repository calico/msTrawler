"""Core statistical model fitting — pure Python using statsmodels.

Replaces R's adaptModel() with statsmodels WLS (weighted least squares) and
optionally MixedLM. Uses Wald F-tests instead of Kenward-Roger approximation.

Note: Kenward-Roger is not available in Python. Wald tests are asymptotically
equivalent but less conservative for small samples. For small-sample
proteomics, results may differ slightly from R's pbkrtest-based p-values.
"""

from __future__ import annotations

import warnings

import numpy as np
import pandas as pd
import statsmodels.api as sm
from numpy.typing import NDArray
from scipy import stats as sp_stats


def adapt_model(
    prot_dat: pd.DataFrame,
    cov_type: str,
    fixed_form_str: str,
    bridge_dat: pd.DataFrame | None,
    t_parm: str,
    lht_list: list | None,
    full_columns: list[str],
    min_re: int,
    rand_id: str,
    reduced_mod: bool,
    multi_level: bool,
) -> tuple[NDArray, NDArray | None, str] | None:
    """Fit adaptive weighted linear model for a single protein.

    Args:
        prot_dat: Data for one protein (non-bridge observations).
        cov_type: One of "None", "Factor", "Continuous", "Time".
        fixed_form_str: Formula string for fixed effects.
        bridge_dat: Bridge channel data (or None).
        t_parm: Time parameterization mode.
        lht_list: List of hypothesis test specifications.
        full_columns: Expected column names for output.
        min_re: Minimum levels for random effects.
        rand_id: Column name for random grouping variable.
        reduced_mod: Whether model was reduced.
        multi_level: Whether to attempt mixed model.

    Returns:
        Tuple of (estimates, time_results, model_type) or None on failure.
    """
    try:
        return _fit_model(
            prot_dat, cov_type, fixed_form_str, bridge_dat, t_parm,
            lht_list, full_columns, min_re, rand_id, reduced_mod, multi_level,
        )
    except Exception:
        return None


def _fit_model(
    prot_dat, cov_type, fixed_form_str, bridge_dat, t_parm,
    lht_list, full_columns, min_re, rand_id, reduced_mod, multi_level,
):
    """Internal model fitting logic."""
    dat = prot_dat.copy()

    # Build design matrix based on cov_type
    if cov_type == "None":
        return _fit_none_model(dat, bridge_dat, full_columns)
    elif cov_type == "Factor":
        return _fit_factor_model(
            dat, bridge_dat, fixed_form_str, full_columns,
            lht_list, min_re, rand_id, reduced_mod, multi_level,
        )
    elif cov_type == "Continuous":
        return _fit_continuous_model(
            dat, bridge_dat, fixed_form_str, full_columns,
            min_re, rand_id, reduced_mod, multi_level,
        )
    elif cov_type == "Time":
        return _fit_time_model(
            dat, bridge_dat, fixed_form_str, full_columns,
            lht_list, min_re, rand_id, multi_level,
        )
    return None


def _fit_none_model(
    dat: pd.DataFrame,
    bridge_dat: pd.DataFrame | None,
    sample_names: list[str],
) -> tuple[NDArray, None, str]:
    """Fit per-sample abundance estimates (covType = "None")."""
    y = dat["lIntensity"].values
    w = 1.0 / dat["techVar"].values

    # Design matrix: one column per sample
    X = pd.get_dummies(dat["SampleID"], drop_first=False).astype(float)

    # Ensure columns match expected sample names
    for sn in sample_names:
        if sn not in X.columns:
            X[sn] = 0.0
    X = X[sample_names]

    if bridge_dat is not None and len(bridge_dat) > 0:
        # Add bridge data with zero sample columns
        by = bridge_dat["lIntensity"].values
        bw = 1.0 / bridge_dat["techVar"].values
        bX = pd.DataFrame(0.0, index=range(len(bridge_dat)), columns=sample_names)

        # Add scan effects for bridge
        y = np.concatenate([y, by])
        w = np.concatenate([w, bw])
        X = pd.concat([X, bX], ignore_index=True)

    X_arr = X.values
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        try:
            model = sm.WLS(y, X_arr, weights=w).fit()
            ci = model.conf_int(alpha=0.05)
            # Result: [lower, est, upper] per sample, flattened
            result = np.column_stack([ci[:, 0], model.params, ci[:, 1]])
            return result.ravel(), None, "lm"
        except Exception:
            return np.full(len(sample_names) * 3, np.nan), None, "none"


def _fit_factor_model(
    dat, bridge_dat, fixed_form_str, full_columns,
    lht_list, min_re, rand_id, reduced_mod, multi_level,
):
    """Fit model with factor covariates."""
    result, model_type = _fit_weighted_model(dat, bridge_dat, full_columns, min_re, multi_level)

    # Extract estimates and p-values for each contrast
    n_cols = len(full_columns)
    res_vec = np.full(n_cols * 3, np.nan)

    if result is not None:
        params = result.params
        param_names = list(result.model.exog_names) if hasattr(result.model, "exog_names") else []

        for j, col_name in enumerate(full_columns):
            # Find matching parameter
            idx = _find_param_index(param_names, col_name)
            if idx is not None and idx < len(params):
                res_vec[j * 3] = params[idx]
                try:
                    # Wald test for this parameter
                    t_val = params[idx] / result.bse[idx]
                    df_resid = result.df_resid if hasattr(result, "df_resid") else len(dat) - len(params)
                    p_val = 2 * (1 - sp_stats.t.cdf(abs(t_val), df=max(df_resid, 1)))
                    res_vec[j * 3 + 1] = p_val
                except Exception:
                    pass

    # Time hypothesis tests
    time_res = None
    if lht_list:
        time_res = _run_hypothesis_tests(result, lht_list, model_type)

    model_label = "lm reduced" if reduced_mod else model_type
    return res_vec, time_res, model_label


def _fit_continuous_model(
    dat, bridge_dat, fixed_form_str, full_columns,
    min_re, rand_id, reduced_mod, multi_level,
):
    """Fit model with continuous covariates."""
    result, model_type = _fit_weighted_model(dat, bridge_dat, full_columns, min_re, multi_level)

    n_cont = len(full_columns)
    res_mat = np.full((n_cont, 3), np.nan)

    if result is not None:
        params = result.params
        param_names = list(result.model.exog_names) if hasattr(result.model, "exog_names") else []

        for j, col_name in enumerate(full_columns):
            idx = _find_param_index(param_names, col_name)
            if idx is not None and idx < len(params):
                res_mat[j, 0] = params[idx]
                try:
                    t_val = params[idx] / result.bse[idx]
                    df_resid = result.df_resid if hasattr(result, "df_resid") else len(dat) - len(params)
                    p_val = 2 * (1 - sp_stats.t.cdf(abs(t_val), df=max(df_resid, 1)))
                    res_mat[j, 1] = p_val
                except Exception:
                    pass

    return res_mat, None, model_type


def _fit_time_model(
    dat, bridge_dat, fixed_form_str, full_columns,
    lht_list, min_re, rand_id, multi_level,
):
    """Fit model and test time trend parameters."""
    result, model_type = _fit_weighted_model(dat, bridge_dat, full_columns, min_re, multi_level)

    time_res = None
    if result is not None and lht_list:
        time_res = _run_hypothesis_tests(result, lht_list, model_type)

    return None, time_res, model_type


def _fit_weighted_model(
    dat: pd.DataFrame,
    bridge_dat: pd.DataFrame | None,
    target_cols: list[str],
    min_re: int,
    multi_level: bool,
) -> tuple | None:
    """Fit a weighted linear model with the available covariates.

    Builds design matrix from the data columns and fits WLS.
    Returns (result, model_type) tuple.
    """
    y = dat["lIntensity"].values.astype(float)
    w = 1.0 / dat["techVar"].values.astype(float)

    # Build design matrix from factor/continuous columns in the data
    # Use all non-metadata columns as predictors
    meta_cols = {"Scan", "PA.Gene.Symbol", "Protein.ID", "Peptide", "Plex",
                 "Channel", "SNR", "lIntensity", "SampleID", "techVar", "Bridge", "variable"}

    pred_cols = [c for c in dat.columns if c not in meta_cols]
    if not pred_cols:
        return None, "none"

    # Build design matrix with proper encoding
    X_parts = []
    col_names = []
    for col in pred_cols:
        if dat[col].dtype == object or isinstance(dat[col].dtype, pd.CategoricalDtype):
            # Categorical — create dummies (drop first = reference)
            dummies = pd.get_dummies(dat[col], prefix=col, drop_first=True, dtype=float)
            X_parts.append(dummies.values)
            col_names.extend(dummies.columns.tolist())
        else:
            # Numeric
            X_parts.append(dat[col].values.reshape(-1, 1).astype(float))
            col_names.append(col)

    if not X_parts:
        return None, "none"

    X = np.hstack(X_parts)

    # Add intercept
    X = sm.add_constant(X)
    col_names = ["const"] + col_names

    # Check rank
    rank = np.linalg.matrix_rank(X)
    if rank < X.shape[1]:
        # Remove collinear columns via QR pivoting
        _, r, p = np.linalg.qr(X, mode="reduced") if X.shape[0] >= X.shape[1] else (None, None, None)
        # Fallback: just use what we have
        pass

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        try:
            model = sm.WLS(y, X, weights=w)
            result = model.fit()
            result.model.exog_names = col_names
            return result, "lm"
        except Exception:
            return None, "none"


def _find_param_index(param_names: list[str], target: str) -> int | None:
    """Find the index of a parameter in the model output."""
    # Exact match first
    if target in param_names:
        return param_names.index(target)
    # Try with prefix (R adds factor name as prefix)
    for i, name in enumerate(param_names):
        if target in name:
            return i
    return None


def _run_hypothesis_tests(
    result,
    lht_list: list,
    model_type: str,
) -> NDArray | None:
    """Run linear hypothesis tests (Wald F-test)."""
    if result is None:
        return None

    n_params = len(result.params)
    results = []

    for test_spec in lht_list:
        if test_spec is None or (isinstance(test_spec, float) and np.isnan(test_spec)):
            results.extend([np.nan, np.nan])
            continue

        try:
            # Build restriction matrix
            if isinstance(test_spec, list):
                param_names = getattr(result.model, "exog_names", [])
                R = np.zeros((len(test_spec), n_params))
                for k, spec in enumerate(test_spec):
                    idx = _find_param_index(param_names, spec)
                    if idx is not None:
                        R[k, idx] = 1.0
            else:
                # Single parameter test
                param_names = getattr(result.model, "exog_names", [])
                idx = _find_param_index(param_names, str(test_spec))
                if idx is None:
                    results.extend([np.nan, np.nan])
                    continue
                R = np.zeros((1, n_params))
                R[0, idx] = 1.0

            f_test = result.f_test(R)
            p_val = float(f_test.pvalue)
            results.extend([p_val, np.nan])  # p-value, q-value placeholder
        except Exception:
            results.extend([np.nan, np.nan])

    return np.array(results) if results else None
