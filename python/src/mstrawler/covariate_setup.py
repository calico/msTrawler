"""Set up covariates, factors, time variables, and model formula from sample/covariate files."""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import pandas as pd


@dataclass
class CovariateConfig:
    """All covariate-related configuration derived from sample and covariate files."""

    sample_file: pd.DataFrame
    bridge_mod: bool
    bridge_col: list[int]
    table_names: list[str]
    co_vars: list[str]
    n_cat: int
    factor_names: list[str]
    level_names: list[list[str]] | None
    n_levels: list[int] | int
    n_cont: int
    cont_names: list[str]
    time_vars: list[str] | None
    time_degree: int
    time_level_n: int
    time_table_names: list[str] | None
    t_cat_index: int
    t_cat_factor_index: int
    t_cat_name: str | None
    t_cat_levels: list[str] | None
    t_parm: str
    circadian: int
    rand_id: str
    fixed_str: str


def setup_covariates(
    sample_file: pd.DataFrame,
    covariate_file: pd.DataFrame,
    u_plex: np.ndarray,
) -> CovariateConfig:
    """Set up covariates, factors, time variables, and model formula."""
    sf = sample_file.copy()

    # Reduce plexes to those present in data
    sf_plexes = sf["SampleID"].str.extract(r"^([^_]+)_")[0].values
    used_plexes = np.intersect1d(u_plex, sf_plexes)

    # Bridge detection
    n_bridges = int(sf["Bridge"].sum())
    bridge_mod = n_bridges > 0
    if bridge_mod and n_bridges != len(used_plexes):
        raise ValueError("The number of plexes must match the number of bridge samples.")

    table_names = ["Simple.csv"]

    # Identify bridge column and build covariate list
    bridge_col = [i for i, c in enumerate(sf.columns) if c.upper() == "BRIDGE"]
    if bridge_col:
        exclude = {0} | set(bridge_col)
        co_vars = [sf.columns[i] for i in range(len(sf.columns)) if i not in exclude]
    else:
        co_vars = list(sf.columns[1:])

    if not co_vars:
        co_vars = ["1"]

    cov_type = covariate_file["Type"].tolist()
    cat_indices = [i for i, t in enumerate(cov_type) if "FACTOR" in str(t).upper()]
    n_cat = len(cat_indices)

    t_parm = ""

    # Factor setup
    if n_cat > 0:
        factor_names = [str(covariate_file["Covariate"].iloc[i]) for i in cat_indices]
        if not all(fn in sf.columns for fn in factor_names):
            raise ValueError("Covariate names in covariateFile don't match sampleFile columns.")

        level_names = []
        n_levels_list = []
        for fn in factor_names:
            if bridge_col:
                real_idx = sf[sf.iloc[:, bridge_col[0]] == 0].index
            else:
                real_idx = sf.index
            levels = sf.loc[real_idx, fn].unique().astype(str).tolist()
            level_names.append(levels)
            n_levels_list.append(len(levels))
            table_names.extend([f"Factor_{fn}_{lv}.csv" for lv in levels])
        n_levels = n_levels_list
    else:
        factor_names = [""]
        level_names = None
        n_levels = 0

    # Continuous covariates
    cont_indices = [i for i, t in enumerate(cov_type) if "CONTINUOUS" in str(t).upper()]
    n_cont = len(cont_indices)
    if n_cont > 0:
        cont_names = [str(covariate_file["Covariate"].iloc[i]) for i in cont_indices]
        table_names.extend([f"Continuous_{cn}.csv" for cn in cont_names])
        for cn in cont_names:
            sf[cn] = sf[cn].astype(float) - sf[cn].astype(float).mean()
    else:
        cont_names = [""]

    # Time covariates
    time_indices = [i for i, t in enumerate(cov_type) if "TIME" in str(t).upper()]
    if time_indices:
        ti = time_indices[0]
        orig_time_name = str(covariate_file["Covariate"].iloc[ti])
        sf = sf.rename(columns={orig_time_name: "Time"})
        covariate_file = covariate_file.copy()
        covariate_file.loc[covariate_file.index[ti], "Covariate"] = "Time"

        if not pd.api.types.is_numeric_dtype(sf["Time"]):
            raise ValueError("Time variable is non-numeric")

        # Remove Time from coVars
        if ti < len(co_vars):
            co_vars = [v for i, v in enumerate(co_vars) if i != ti]

        time_vars = []
        time_degree = int(covariate_file["TimeDegree"].iloc[ti])
        if time_degree >= 1:
            time_vars.append("Time")
        if time_degree >= 2:
            time_vars.append("Time2")
            sf["Time2"] = sf["Time"].astype(float) ** 2
        if time_degree == 3:
            time_vars.append("Time3")
            sf["Time3"] = sf["Time"].astype(float) ** 3

        # Time category
        t_cat_mask = covariate_file["TimeCategory"] == 1
        t_cat_indices = covariate_file.index[t_cat_mask].tolist()
        if len(t_cat_indices) > 1:
            raise ValueError("Only one time category is allowed")

        if len(t_cat_indices) == 1:
            t_cat_idx = t_cat_indices[0]
            t_cat_name_val = str(covariate_file["Covariate"].iloc[t_cat_idx])
            t_cat_factor_idx = factor_names.index(t_cat_name_val) if t_cat_name_val in factor_names else -1
            if t_cat_factor_idx >= 0:
                t_cat_index = t_cat_idx + 1  # 1-based like R
                t_cat_factor_index = t_cat_factor_idx + 1  # 1-based like R
                t_cat_name = t_cat_name_val
                time_level_n = n_levels[t_cat_factor_idx] if isinstance(n_levels, list) else 0
                t_cat_levels = level_names[t_cat_factor_idx] if level_names else None
                time_table_names = [f"Time_{lv}.csv" for lv in t_cat_levels] if t_cat_levels else None
                t_parm = "Category"
            else:
                t_cat_index = 0
                t_cat_factor_index = 0
                t_cat_name = None
                time_level_n = 1
                time_table_names = ["Time.csv"]
                t_cat_levels = None
                t_parm = "Continuous" if n_cont > 0 else "Time"
        else:
            t_cat_index = 0
            t_cat_factor_index = 0
            t_cat_name = None
            time_level_n = 1
            time_table_names = ["Time.csv"]
            t_cat_levels = None
            t_parm = "Continuous" if n_cont > 0 else "Time"

        circadian = int(covariate_file["Circadian"].iloc[ti])
    else:
        time_degree = 0
        time_vars = None
        time_level_n = 0
        t_cat_index = 0
        t_cat_factor_index = 0
        circadian = 0
        t_parm = ""
        time_table_names = None
        t_cat_name = None
        t_cat_levels = None

    # Circadian terms
    if circadian != 0:
        sf["Sine"] = np.sin((2 * np.pi / 24) * sf["Time"].astype(float))
        sf["Cosine"] = np.cos((2 * np.pi / 24) * sf["Time"].astype(float))
        if time_vars is None:
            time_vars = []
        time_vars.extend(["Sine", "Cosine"])

    # Longitudinal ID
    id_indices = [i for i, t in enumerate(cov_type) if "ID" in str(t).upper()]
    if id_indices:
        id_name = str(covariate_file["Covariate"].iloc[id_indices[0]])
        rand_id = id_name
        co_vars = [v for v in co_vars if v != id_name]
    else:
        rand_id = "SampleID"

    # Build formula string
    all_terms = co_vars + (time_vars or [])
    fixed_str = "lIntensity ~ " + " + ".join(all_terms)

    if t_cat_index > 0 and t_cat_name and time_vars:
        interaction_terms = [f"{t_cat_name}:{tv}" for tv in time_vars]
        fixed_str += " + " + " + ".join(interaction_terms)

    return CovariateConfig(
        sample_file=sf,
        bridge_mod=bridge_mod,
        bridge_col=bridge_col,
        table_names=table_names,
        co_vars=co_vars,
        n_cat=n_cat,
        factor_names=factor_names,
        level_names=level_names,
        n_levels=n_levels,
        n_cont=n_cont,
        cont_names=cont_names,
        time_vars=time_vars,
        time_degree=time_degree,
        time_level_n=time_level_n,
        time_table_names=time_table_names,
        t_cat_index=t_cat_index,
        t_cat_factor_index=t_cat_factor_index,
        t_cat_name=t_cat_name,
        t_cat_levels=t_cat_levels,
        t_parm=t_parm,
        circadian=circadian,
        rand_id=rand_id,
        fixed_str=fixed_str,
    )
