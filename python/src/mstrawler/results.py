"""Result table initialization and output."""

from __future__ import annotations

import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests


def init_result_tables(
    ready_df: pd.DataFrame,
    sample_names: list[str],
    n_cat: int,
    factor_names: list[str],
    level_names: list[list[str]] | None,
    n_levels: list[int] | int,
    n_cont: int,
    cont_names: list[str],
    time_level_n: int,
    time_vars: list[str] | None,
    t_cat_levels: list[str] | None,
) -> dict:
    """Initialize empty results tables for protein modeling output."""
    u_prot = ready_df["Protein.ID"].dropna().unique()
    gene_map = ready_df.dropna(subset=["Protein.ID"]).groupby("Protein.ID")["PA.Gene.Symbol"].first()
    u_gene = gene_map.reindex(u_prot).values
    n_prot = len(u_prot)

    table_list = []

    # Sample estimate table
    res_cols = []
    for sn in sample_names:
        res_cols.extend([f"lower_{sn}", f"est_{sn}", f"upper_{sn}"])
    sample_tab = pd.DataFrame(
        np.nan,
        index=range(n_prot),
        columns=["Gene", "Protein", "modelFit"] + res_cols,
    )
    sample_tab["Gene"] = u_gene
    sample_tab["Protein"] = u_prot
    table_list.append(sample_tab)

    # Factor tables
    if n_cat > 0 and level_names is not None:
        n_levels_list = n_levels if isinstance(n_levels, list) else [n_levels]
        for i, fn in enumerate(factor_names):
            if fn == "":
                continue
            lvls = level_names[i]
            for j in range(len(lvls)):
                non_ref = [lv for k, lv in enumerate(lvls) if k != j]
                cols = []
                for lv in non_ref:
                    cols.extend([f"Est_{lv}", f"Pval_{lv}", f"Qval_{lv}"])
                tab = pd.DataFrame(
                    np.nan, index=range(n_prot),
                    columns=["Gene", "Protein", "modelFit"] + cols,
                )
                tab["Gene"] = u_gene
                tab["Protein"] = u_prot
                table_list.append(tab)

    # Continuous tables
    if n_cont > 0:
        for cn in cont_names:
            if cn == "":
                continue
            tab = pd.DataFrame(
                np.nan, index=range(n_prot),
                columns=["Gene", "Protein", "modelFit", f"Est_{cn}", f"Pval_{cn}", f"Qval_{cn}"],
            )
            tab["Gene"] = u_gene
            tab["Protein"] = u_prot
            table_list.append(tab)

    # Time tables
    time_tables = []
    if time_level_n > 0 and time_vars:
        if time_level_n == 1:
            cols = ["Gene", "Protein", "modelFit"] + time_vars + ["Pval_Time", "Qval_Time"]
            tab = pd.DataFrame(np.nan, index=range(n_prot), columns=cols)
            tab["Gene"] = u_gene
            tab["Protein"] = u_prot
            time_tables.append(tab)
        elif t_cat_levels:
            for l_idx in range(time_level_n):
                ref_and_others = ["REF"] + [lv for k, lv in enumerate(t_cat_levels) if k != l_idx]
                p_q_cols = []
                for lv in ref_and_others:
                    p_q_cols.extend([f"Pval_Time_{lv}", f"Qval_Time_{lv}"])
                cols = ["Gene", "Protein", "modelFit"] + time_vars + p_q_cols
                tab = pd.DataFrame(np.nan, index=range(n_prot), columns=cols)
                tab["Gene"] = u_gene
                tab["Protein"] = u_prot
                time_tables.append(tab)

    # Ensure modelFit column accepts string values in all tables
    for tab in table_list + time_tables:
        if "modelFit" in tab.columns:
            tab["modelFit"] = tab["modelFit"].astype(object)

    return {
        "table_list": table_list,
        "time_tables": time_tables,
        "u_prot": u_prot,
        "u_gene": u_gene,
        "n_prot": n_prot,
    }


def write_result_tables(
    table_list: list[pd.DataFrame],
    table_names: list[str],
    time_tables: list[pd.DataFrame],
    time_table_names: list[str] | None,
    apply_fdr: bool = True,
    suffix: str | None = None,
) -> None:
    """Write result tables to CSV, optionally applying FDR correction."""
    for i, tab in enumerate(table_list):
        q_cols = [c for c in tab.columns if "Qval" in c]
        for qc in q_cols:
            p_col = qc.replace("Qval", "Pval")
            if p_col in tab.columns:
                p_vals = tab[p_col].values.astype(float)
                # Replace 0 p-values
                p_vals[p_vals == 0] = 1e-300
                tab[p_col] = p_vals
                if apply_fdr:
                    valid = ~np.isnan(p_vals)
                    if valid.sum() > 0:
                        _, q_vals, _, _ = multipletests(p_vals[valid], method="fdr_bh")
                        tab.loc[valid, qc] = q_vals

        filename = table_names[i]
        if suffix:
            filename = filename.replace(".csv", f"___{suffix}.csv")
        tab.to_csv(filename, index=False)

    if time_tables and time_table_names:
        for i, tab in enumerate(time_tables):
            q_cols = [c for c in tab.columns if "Qval" in c]
            for qc in q_cols:
                p_col = qc.replace("Qval", "Pval")
                if p_col in tab.columns:
                    p_vals = tab[p_col].values.astype(float)
                    p_vals[p_vals == 0] = 1e-300
                    tab[p_col] = p_vals
                    if apply_fdr:
                        valid = ~np.isnan(p_vals)
                        if valid.sum() > 0:
                            _, q_vals, _, _ = multipletests(p_vals[valid], method="fdr_bh")
                            tab.loc[valid, qc] = q_vals

            filename = time_table_names[i]
            if suffix:
                filename = filename.replace(".csv", f"___{suffix}.csv")
            tab.to_csv(filename, index=False)
