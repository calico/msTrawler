"""Convert Proteome Discoverer output files to msTrawler input format."""

from __future__ import annotations

import re
from pathlib import Path

import numpy as np
import pandas as pd


def convert_pd_file(
    pd_psm_file: str | Path,
    pd_protein_file: str | Path = "",
    delimiter: str = "\t",
) -> pd.DataFrame:
    """Convert Proteome Discoverer PSM/Protein exports to msTrawler dataframe.

    Args:
        pd_psm_file: Path to PSM tab-delimited file from Proteome Discoverer.
        pd_protein_file: Optional protein file for gene name mapping.
        delimiter: File delimiter (default: tab).

    Returns:
        DataFrame with columns: PA.Gene.Symbol, Protein.ID, Peptide, Plex,
        X{channel}.Adjusted.Intensity (per channel),
        X{channel}.Sn (per channel).

    Raises:
        FileNotFoundError: If pd_psm_file does not exist.
        ValueError: If required columns are missing.
    """
    pd_psm_path = Path(pd_psm_file)
    if not pd_psm_path.exists():
        raise FileNotFoundError(f"PSM file not found: {pd_psm_file}")

    # Build protein→gene mapping if protein file provided
    protein_gene_map: dict[str, str] = {}
    pd_protein_path = Path(pd_protein_file) if pd_protein_file else None
    if pd_protein_path and pd_protein_path.exists():
        protein_df = pd.read_csv(pd_protein_path, sep=delimiter)
        # Column names may have dots (R read.csv) or spaces
        gene_col = _find_column(protein_df, "Gene.Symbol", "Gene Symbol")
        acc_col = _find_column(protein_df, "Accession")
        protein_gene_map = dict(zip(protein_df[acc_col], protein_df[gene_col]))

    # Parse the PSM file line by line (matching R behavior for large files)
    rows = []
    channels = None
    intensity_cols = None
    sn_cols = None

    with open(pd_psm_path, "r") as f:
        header_line = f.readline().strip().replace('"', "")
        if not header_line:
            raise ValueError("PSM file is empty")

        headers = header_line.split(delimiter)

        # Find column indices
        avg_sn_idx = _grep_index(headers, "Average Reporter S")
        protein_idx = _grep_indices(headers, "Protein Accessions")
        seq_idx = _exact_index(headers, "Annotated Sequence")
        ptm_idx = _exact_index(headers, "Modifications")
        plex_idx = _exact_index(headers, "File ID")
        abd_indices = [i for i, h in enumerate(headers) if "Abundance" in h]

        if avg_sn_idx is None:
            raise ValueError("Missing 'Average Reporter S/N' column")
        if not protein_idx:
            raise ValueError("Missing 'Protein Accessions' column")
        if seq_idx is None:
            raise ValueError("Missing 'Annotated Sequence' column")
        if plex_idx is None:
            raise ValueError("Missing 'File ID' column")
        if not abd_indices:
            raise ValueError("Missing Abundance columns")

        # Extract channel names: "Abundance: 126" or "Abundance 126" → "126"
        channels = [re.sub(r"Abundance:?\s*", "", headers[i]).lower() for i in abd_indices]
        intensity_cols = [f"X{ch}.Adjusted.Intensity" for ch in channels]
        sn_cols = [f"X{ch}.Sn" for ch in channels]

        # Process rows
        for line in f:
            line = line.strip()
            if not line:
                continue

            fields = line.split(delimiter)
            fields = [f.replace('"', "") for f in fields]

            # Skip if average SN is NA, empty, or 0
            avg_sn_str = fields[avg_sn_idx] if avg_sn_idx < len(fields) else ""
            if not avg_sn_str or avg_sn_str == "NA":
                continue
            try:
                avg_sn = float(avg_sn_str)
            except ValueError:
                continue
            if avg_sn == 0:
                continue

            # Extract protein accession (use first match)
            protein = fields[protein_idx[0]]
            gene = protein_gene_map.get(protein, "")

            # Clean peptide sequence (remove brackets)
            peptide = fields[seq_idx].replace("[", "").replace("]", "")
            peptide = f"{peptide} {fields[ptm_idx]}"

            plex = fields[plex_idx]

            # Extract abundance values, convert NA to 0
            intensities = np.zeros(len(abd_indices))
            for k, idx in enumerate(abd_indices):
                val = fields[idx] if idx < len(fields) else ""
                try:
                    intensities[k] = float(val)
                except (ValueError, TypeError):
                    intensities[k] = 0.0

            # Calculate per-channel SNR from average reporter SN
            noise = np.mean(intensities) / avg_sn
            if noise > 0:
                snr = intensities / noise
            else:
                snr = np.zeros(len(abd_indices))

            row = [gene, protein, peptide, plex] + intensities.tolist() + snr.tolist()
            rows.append(row)

    # Build dataframe
    columns = ["PA.Gene.Symbol", "Protein.ID", "Peptide", "Plex"] + intensity_cols + sn_cols
    df = pd.DataFrame(rows, columns=columns)

    # Match R behavior: empty strings become NaN (R's read.csv converts "" to NA)
    df["PA.Gene.Symbol"] = df["PA.Gene.Symbol"].replace("", np.nan)
    df["Protein.ID"] = df["Protein.ID"].replace("", np.nan)

    # Write converted CSV (matching R behavior)
    output_path = str(pd_psm_file) + "_converted.csv"
    df.to_csv(output_path, index=False)

    return df


def _grep_index(headers: list[str], pattern: str) -> int | None:
    """Find first index where pattern appears in header (case-sensitive substring match)."""
    for i, h in enumerate(headers):
        if pattern in h:
            return i
    return None


def _grep_indices(headers: list[str], pattern: str) -> list[int]:
    """Find all indices where pattern appears in header."""
    return [i for i, h in enumerate(headers) if pattern in h]


def _exact_index(headers: list[str], name: str) -> int | None:
    """Find index where header exactly matches name."""
    for i, h in enumerate(headers):
        if h == name:
            return i
    return None


def _find_column(df: pd.DataFrame, *candidates: str) -> str:
    """Find the first matching column name from candidates (handles dots vs spaces)."""
    for name in candidates:
        if name in df.columns:
            return name
        # Try with dots replaced by spaces and vice versa
        alt = name.replace(".", " ")
        if alt in df.columns:
            return alt
        alt = name.replace(" ", ".")
        if alt in df.columns:
            return alt
    raise ValueError(f"None of {candidates} found in columns: {list(df.columns)}")
