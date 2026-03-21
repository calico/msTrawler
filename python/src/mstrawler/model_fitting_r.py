"""R-backed model fitting via rpy2 — uses lme4, pbkrtest, and car.

Requires: R installed, plus R packages lme4, pbkrtest, car, reshape2.
Install the optional dependency: pip install mstrawler[r]
"""

from __future__ import annotations

import warnings

import numpy as np
import pandas as pd
from numpy.typing import NDArray

try:
    import rpy2.robjects as ro
    from rpy2.robjects import pandas2ri, Formula
    from rpy2.robjects.packages import importr
    from rpy2.robjects.conversion import localconverter

    _RPY2_AVAILABLE = True
except ImportError:
    _RPY2_AVAILABLE = False


def _check_rpy2():
    if not _RPY2_AVAILABLE:
        raise ImportError(
            "rpy2 is required for the 'r' backend. "
            "Install with: pip install mstrawler[r]\n"
            "Also requires R with packages: lme4, pbkrtest, car"
        )


def _ensure_r_packages():
    """Load required R packages, raising a clear error if missing."""
    _check_rpy2()
    try:
        importr("lme4")
        importr("pbkrtest")
        importr("car")
    except Exception as e:
        raise ImportError(
            f"Required R package not found: {e}\n"
            "Install in R: install.packages(c('lme4', 'pbkrtest', 'car'))"
        )


def _pd_to_r(df: pd.DataFrame) -> ro.DataFrame:
    """Convert pandas DataFrame to R DataFrame."""
    with localconverter(ro.default_converter + pandas2ri.converter):
        return ro.conversion.py2rpy(df)


def _r_to_py(obj):
    """Convert R object to Python."""
    with localconverter(ro.default_converter + pandas2ri.converter):
        return ro.conversion.rpy2py(obj)


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
    """Fit model using R's lme4/pbkrtest via rpy2.

    Same interface as model_fitting_python.adapt_model, but calls
    the original R adaptModel() function for full statistical fidelity
    including mixed-effects models and Kenward-Roger tests.
    """
    _check_rpy2()
    _ensure_r_packages()

    try:
        return _fit_via_r(
            prot_dat, cov_type, fixed_form_str, bridge_dat, t_parm,
            lht_list, full_columns, min_re, rand_id, reduced_mod, multi_level,
        )
    except Exception:
        return None


def _fit_via_r(
    prot_dat, cov_type, fixed_form_str, bridge_dat, t_parm,
    lht_list, full_columns, min_re, rand_id, reduced_mod, multi_level,
):
    """Call R's adaptModel via rpy2."""
    base = importr("base")
    stats = importr("stats")

    # Source msTrawler's adaptModel if not already loaded
    # We load the installed package
    try:
        mstrawler_r = importr("msTrawler")
    except Exception:
        raise ImportError(
            "R package 'msTrawler' not found. Install with:\n"
            "  R CMD INSTALL /path/to/msTrawler --no-staged-install"
        )

    # Convert Python DataFrames to R
    r_prot_dat = _pd_to_r(prot_dat)

    r_bridge_dat = ro.NULL
    if bridge_dat is not None and len(bridge_dat) > 0:
        r_bridge_dat = _pd_to_r(bridge_dat)

    # Convert fixed_form_str to R formula
    r_fixed_form = Formula(fixed_form_str)

    # Convert lht_list to R list
    if lht_list is not None:
        r_lht = ro.ListVector({})
        for i, item in enumerate(lht_list):
            if item is None or (isinstance(item, float) and np.isnan(item)):
                r_lht.rx2[i + 1] = ro.NA_Real
            elif isinstance(item, list):
                r_lht.rx2[i + 1] = ro.StrVector(item)
            else:
                r_lht.rx2[i + 1] = ro.StrVector([str(item)])
    else:
        r_lht = ro.NULL

    r_full_columns = ro.StrVector(full_columns)

    # Call R's adaptModel
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        result = mstrawler_r.adaptModel(
            protDat=r_prot_dat,
            covType=cov_type,
            fixedForm=r_fixed_form,
            bridgeDat=r_bridge_dat,
            timeParm=t_parm,
            lhtList=r_lht,
            fullColumns=r_full_columns,
            minRE=min_re,
            randID=rand_id,
            reducedMod=reduced_mod,
            multiLevel=multi_level,
        )

    # Convert R result list back to Python
    # result[[1]] = estimates, result[[2]] = time results, result[[3]] = model type
    estimates = np.array(_r_to_py(result.rx2(1)))
    time_res = None
    r_time = result.rx2(2)
    if r_time != ro.NULL and r_time is not None:
        time_res = np.array(_r_to_py(r_time))
    model_type = str(result.rx2(3)[0])

    return estimates, time_res, model_type
