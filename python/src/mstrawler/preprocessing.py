"""Preprocessing functions: LOD imputation, normalization, and outlier detection."""

from __future__ import annotations

import numpy as np
import pandas as pd
import statsmodels.api as sm
from numpy.typing import NDArray


def tmt_lod(
    sn_mat: NDArray,
    i_mat: NDArray,
    lod: float = 0.01,
    min_above: int = 3,
    scale_sn: float = 1.0,
    impute_penalty: float = 1.0,
    rng: np.random.Generator | None = None,
) -> tuple[NDArray, NDArray, NDArray]:
    """Impute below-LOD values from a log-normal distribution.

    Args:
        sn_mat: Signal-to-noise matrix (scans x channels).
        i_mat: Intensity matrix (scans x channels).
        lod: Limit of detection threshold.
        min_above: Minimum channels above LOD to keep a scan.
        scale_sn: SNR scaling factor.
        impute_penalty: Penalty for imputed values.
        rng: Random number generator (for reproducibility).

    Returns:
        Tuple of (imputed_intensities, imputed_snr, too_few_flags).
    """
    if rng is None:
        rng = np.random.default_rng()

    sn = sn_mat.copy().astype(float)
    im = i_mat.copy().astype(float)

    # Row-wise sums
    ssn = np.nansum(sn, axis=1)
    ssi = np.nansum(im, axis=1)

    # Proportion matrix
    with np.errstate(divide="ignore", invalid="ignore"):
        pi_mat = im / ssi[:, np.newaxis]

    # Determine censored entries
    if lod < 1:
        censored = (pi_mat < lod) | np.isnan(pi_mat)
    else:
        censored = (sn < lod) & ~np.isnan(pi_mat)

    # Flag scans with too few above-LOD channels
    n_above = np.sum(~censored, axis=1)
    too_few = n_above < min_above

    # Compute imputation parameters
    n_cols = im.shape[1]
    if lod < 1:
        e_miss = ssi * lod / 2
        e_mat = np.tile(e_miss[:, np.newaxis], (1, n_cols))
        sd_miss = ssn * lod / (2 * impute_penalty)
        sd_mat = np.tile(sd_miss[:, np.newaxis], (1, n_cols))
    else:
        # Find rank of SNR nearest to LOD
        ranked_sn = np.sort(sn.ravel())
        sn_index = np.max(np.where(ranked_sn > lod))
        ranked_i = np.sort(im.ravel())
        lod_i = ranked_i[sn_index]

        e_miss_val = lod_i / 2
        e_mat = np.full_like(im, e_miss_val)
        sd_miss_val = lod / 2
        sd_mat = np.full_like(sn, sd_miss_val)

    # Generate imputed values for censored entries
    sd_vec = sd_mat[censored]
    e_vec = e_mat[censored]
    with np.errstate(divide="ignore", invalid="ignore"):
        imp_vec = 2 ** rng.normal(np.log2(e_vec), 1.0 / np.sqrt(scale_sn * sd_vec))

    # Replace censored values
    new_sn = sn.copy()
    new_sn[censored] = sd_vec

    new_i = im.copy()
    new_i[censored] = imp_vec

    return new_i, new_sn, too_few


def geo_norm(
    mat: NDArray,
    norm_index: NDArray,
    plex: NDArray,
    ratios: NDArray,
) -> tuple[NDArray, NDArray]:
    """Geometric mean normalization across channels within each plex.

    Args:
        mat: Reporter ion intensity matrix (scans x channels).
        norm_index: Boolean/int vector (1 = use for normalization).
        plex: Vector of plex membership per scan.
        ratios: Expected channel ratios (typically all 1s).

    Returns:
        Tuple of (normalized_matrix, normalization_factors).
    """
    mat = mat.copy().astype(float)
    mat[mat == 0] = np.nan

    l_mat = np.log2(mat)
    u_plex = np.unique(plex)
    n_channels = mat.shape[1]

    # Compute channel means per plex using normalization rows
    ch_means = np.zeros((n_channels, len(u_plex)))
    for x, p in enumerate(u_plex):
        mask = (norm_index == 1) & (plex == p)
        ch_means[:, x] = np.nanmean(l_mat[mask, :], axis=0)

    # Reshape ratios to match
    log_ratios = np.log2(ratios.reshape(n_channels, len(u_plex), order="F"))
    grand_mean = np.nanmean(ch_means - log_ratios)
    norm_factors = ch_means - log_ratios - grand_mean

    # Apply normalization
    new_mat = mat.copy()
    for x, p in enumerate(u_plex):
        mask = plex == p
        new_mat[mask, :] = 2 ** (l_mat[mask, :] - norm_factors[:, x])

    return new_mat, norm_factors


def find_outliers(
    prot_mat: NDArray,
    sn_mat: NDArray,
    outlier_cutoff: float = 3.0,
    scale_sn: float = 1.0,
) -> NDArray:
    """Detect outlier scans using consensus standardized residuals.

    Fits two weighted linear models (channel model and intercept model)
    and flags entries where |standardized residual| > cutoff in either model.

    Args:
        prot_mat: Intensity matrix for one protein (scans x channels).
        sn_mat: SNR matrix for one protein (scans x channels).
        outlier_cutoff: Standardized residual threshold (default: 3).
        scale_sn: SNR scaling factor.

    Returns:
        Boolean matrix (scans x channels) where 1 = outlier.
    """
    n_scans, n_channels = prot_mat.shape

    # Log-transform and center by row mean
    with np.errstate(divide="ignore", invalid="ignore"):
        log_mat = np.log2(prot_mat.astype(float))
    row_means = np.nanmean(log_mat, axis=1)
    log_mat = log_mat - row_means[:, np.newaxis]

    # Melt to long format
    scan_ids = np.repeat(np.arange(n_scans), n_channels)
    channel_ids = np.tile(np.arange(n_channels), n_scans)
    values = log_mat.ravel()
    weights = (scale_sn * sn_mat.ravel()).astype(float)

    # Remove NaN entries
    valid = ~np.isnan(values)

    # Model 1: value ~ channel (no intercept)
    outlier_mat1 = _fit_and_flag(
        scan_ids, channel_ids, values, weights, valid,
        n_scans, n_channels, outlier_cutoff, include_channel=True,
    )

    # Model 2: value ~ 1 (intercept only)
    outlier_mat2 = _fit_and_flag(
        scan_ids, channel_ids, values, weights, valid,
        n_scans, n_channels, outlier_cutoff, include_channel=False,
    )

    # Consensus: outlier in either model
    return np.maximum(outlier_mat1, outlier_mat2)


def _fit_and_flag(
    scan_ids: NDArray,
    channel_ids: NDArray,
    values: NDArray,
    weights: NDArray,
    valid: NDArray,
    n_scans: int,
    n_channels: int,
    cutoff: float,
    include_channel: bool,
) -> NDArray:
    """Fit a WLS model and return outlier flag matrix."""
    y = values[valid]
    w = weights[valid]

    if include_channel:
        # Design matrix: one dummy per channel (no intercept)
        X = np.zeros((valid.sum(), n_channels))
        ch = channel_ids[valid]
        for j in range(n_channels):
            X[ch == j, j] = 1.0
    else:
        # Intercept only
        X = np.ones((valid.sum(), 1))

    try:
        model = sm.WLS(y, X, weights=w).fit()
        resid = model.get_influence().resid_studentized_internal
    except Exception:
        return np.zeros((n_scans, n_channels), dtype=int)

    # Map residuals back to full matrix
    full_resid = np.full(len(values), np.nan)
    full_resid[valid] = resid
    resid_mat = full_resid.reshape(n_scans, n_channels)

    return (np.abs(resid_mat) > cutoff).astype(int)
