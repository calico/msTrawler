"""Model fitting dispatcher — routes to Python or R backend.

Usage:
    # Pure Python (default, no R dependency)
    from mstrawler.model_fitting import adapt_model
    result = adapt_model(data, ..., backend="python")

    # R backend via rpy2 (requires R + lme4 + pbkrtest)
    result = adapt_model(data, ..., backend="r")

The two backends have identical interfaces but different statistical methods:

    backend="python":
        - statsmodels WLS (weighted least squares)
        - Wald F-tests for hypothesis testing
        - No mixed-effects models
        - No R dependency

    backend="r":
        - lme4::lmer (weighted mixed-effects models)
        - pbkrtest::KRmodcomp (Kenward-Roger tests)
        - car::lht (linear hypothesis tests)
        - Requires R + rpy2 + R packages installed
"""

from __future__ import annotations

from typing import Literal

import numpy as np
import pandas as pd
from numpy.typing import NDArray


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
    backend: Literal["python", "r"] = "python",
) -> tuple[NDArray, NDArray | None, str] | None:
    """Fit adaptive weighted model for a single protein.

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
        backend: "python" for statsmodels WLS, "r" for lme4+pbkrtest via rpy2.

    Returns:
        Tuple of (estimates, time_results, model_type) or None on failure.
    """
    if backend == "r":
        from .model_fitting_r import adapt_model as _adapt_r

        return _adapt_r(
            prot_dat, cov_type, fixed_form_str, bridge_dat, t_parm,
            lht_list, full_columns, min_re, rand_id, reduced_mod, multi_level,
        )
    elif backend == "python":
        from .model_fitting_python import adapt_model as _adapt_py

        return _adapt_py(
            prot_dat, cov_type, fixed_form_str, bridge_dat, t_parm,
            lht_list, full_columns, min_re, rand_id, reduced_mod, multi_level,
        )
    else:
        raise ValueError(f"Unknown backend: {backend!r}. Use 'python' or 'r'.")


def check_backend(backend: str) -> dict:
    """Check if a backend is available and return status info.

    Returns:
        Dict with 'available' (bool), 'backend' (str), 'details' (str).
    """
    if backend == "python":
        return {
            "available": True,
            "backend": "python",
            "details": "statsmodels WLS + Wald F-tests (always available)",
        }
    elif backend == "r":
        try:
            from .model_fitting_r import _check_rpy2, _ensure_r_packages

            _check_rpy2()
            _ensure_r_packages()
            return {
                "available": True,
                "backend": "r",
                "details": "lme4 + pbkrtest + car via rpy2 (R installed and packages found)",
            }
        except ImportError as e:
            return {
                "available": False,
                "backend": "r",
                "details": str(e),
            }
    else:
        return {
            "available": False,
            "backend": backend,
            "details": f"Unknown backend: {backend!r}",
        }
