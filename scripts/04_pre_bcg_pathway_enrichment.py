#!/usr/bin/env python3
# ============================================================
# 04_pre_bcg_pathway_enrichment.py
#
# Pre-BCG pseudobulk differential expression + Hallmark GSEA
# for Figure 1c / Supplementary Figure 11.
#
# Historical analysis conventions recovered from prior code:
#   - raw single-nucleus RNA counts
#   - pseudobulk by sample_id × cell_type
#   - PyDESeq2 / DESeq2-style Wald test
#   - rank genes by signed Wald statistic ("stat")
#   - gseapy.prerank()
#   - Hallmark gene sets
#   - min_size = 10
#   - max_size = 500
#   - permutation_num = 2000
#   - seed = 0
#
# IMPORTANT:
# RAW_DATA_OBJECT.h5ad stores obs["sample_id"] as numeric sample
# numbers (1...32). These are mapped directly to
# metadata/sample_metadata.csv using the sample_number column.
#
# Positive log2FoldChange / Wald statistic / NES =
# higher or enriched in early-recurrence patients.
#
# Run from repository root:
#
#   python scripts\04_pre_bcg_pathway_enrichment.py
#
# Or explicitly:
#
#   python scripts\04_pre_bcg_pathway_enrichment.py `
#       --input-h5ad data\raw\RAW_DATA_OBJECT.h5ad
#
# If raw counts are in a layer rather than .X:
#
#   python scripts\04_pre_bcg_pathway_enrichment.py `
#       --input-h5ad data\raw\RAW_DATA_OBJECT.h5ad `
#       --counts-layer counts
#
# ============================================================

from __future__ import annotations

import argparse
import json
import os
import platform
import re
import sys
from pathlib import Path

import anndata as ad
import gseapy as gp
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy
from scipy import sparse

from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats


# ============================================================
# 1. ANALYSIS PARAMETERS
# ============================================================

MIN_CELLS_PER_SAMPLE_CT = 1

# Recovered from the historical Python pseudobulk workflow.
MIN_TOTAL_COUNTS_GENE = 10
MIN_SAMPLES_GENE = 3

GSEA_MIN_SIZE = 10
GSEA_MAX_SIZE = 500
GSEA_PERMUTATIONS = 2000
GSEA_SEED = 0

N_CPUS = min(8, os.cpu_count() or 1)
FDR_THRESHOLD = 0.05


# Reported Figure 1c values.
# AUDIT TARGETS ONLY. These values never affect the analysis.
MANUSCRIPT_TARGETS = [
    {
        "cell_type": "Activated Monocyte",
        "pathway": "HALLMARK_INTERFERON_GAMMA_RESPONSE",
        "reported_nes": 2.856,
        "reported_fdr": "<0.001",
    },
    {
        "cell_type": "T cell (CD4+)",
        "pathway": "HALLMARK_INTERFERON_ALPHA_RESPONSE",
        "reported_nes": 2.564,
        "reported_fdr": "<0.001",
    },
    {
        "cell_type": "Classical Monocyte",
        "pathway": "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
        "reported_nes": 2.535,
        "reported_fdr": "<0.001",
    },
    {
        "cell_type": "T cell (CD8+)",
        "pathway": "HALLMARK_MTORC1_SIGNALING",
        "reported_nes": 1.872,
        "reported_fdr": "<0.001",
    },
    {
        "cell_type": "B-cell",
        "pathway": "HALLMARK_MYC_TARGETS_V1",
        "reported_nes": 1.785,
        "reported_fdr": "0.002",
    },
]


# ============================================================
# 2. HELPERS
# ============================================================

def safe_name(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9_-]+", "_", str(value)).strip("_")


def ensure_dir(path: Path) -> Path:
    path.mkdir(parents=True, exist_ok=True)
    return path


def as_csr(matrix):
    if sparse.issparse(matrix):
        return matrix.tocsr()
    return sparse.csr_matrix(matrix)


def assert_integer_like_counts(matrix, tolerance: float = 1e-6) -> None:
    """
    Fail if the requested expression matrix looks normalized/log-transformed
    rather than raw count-like data.
    """
    matrix = as_csr(matrix)

    if matrix.data.size == 0:
        raise ValueError("Expression matrix contains no non-zero values.")

    values = matrix.data[: min(50000, matrix.data.size)]

    if np.nanmin(values) < 0:
        raise ValueError(
            "Negative values were detected. This analysis requires raw RNA counts."
        )

    max_fractional_part = np.max(np.abs(values - np.rint(values)))

    if max_fractional_part > tolerance:
        raise ValueError(
            "The selected expression matrix is not integer-like. "
            "This analysis requires raw RNA counts, not normalized/log-transformed data."
        )


def get_counts_matrix(adata: ad.AnnData, counts_layer: str):
    if counts_layer == "X":
        matrix = adata.X
    else:
        if counts_layer not in adata.layers:
            raise KeyError(
                f"Counts layer '{counts_layer}' was not found. "
                f"Available layers: {list(adata.layers.keys())}"
            )
        matrix = adata.layers[counts_layer]

    matrix = as_csr(matrix)
    assert_integer_like_counts(matrix)
    return matrix


def read_gmt(path: Path) -> dict[str, list[str]]:
    gene_sets: dict[str, list[str]] = {}

    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 3:
                gene_sets[fields[0]] = fields[2:]

    if not gene_sets:
        raise ValueError(f"No gene sets could be read from {path}")

    return gene_sets


def build_hallmark_library(
    hallmark_gmt: Path | None,
) -> tuple[dict[str, list[str]], str]:
    """
    Historical code used gseapy.get_library() and tried:
        MSigDB_Hallmark_2020
        MSigDB_Hallmark_2023
        Hallmark

    A local GMT can optionally be supplied later for frozen/offline
    Code Ocean reproducibility.
    """
    if hallmark_gmt is not None:
        if not hallmark_gmt.exists():
            raise FileNotFoundError(f"Hallmark GMT not found: {hallmark_gmt}")

        return read_gmt(hallmark_gmt), str(hallmark_gmt)

    candidates = [
        "MSigDB_Hallmark_2020",
        "MSigDB_Hallmark_2023",
        "Hallmark",
    ]

    errors = []

    for name in candidates:
        try:
            library = gp.get_library(
                name=name,
                organism="Human",
            )

            library = {
                pathway: list(genes)
                for pathway, genes in library.items()
            }

            return library, name

        except Exception as exc:
            errors.append(
                f"{name}: {type(exc).__name__}: {exc}"
            )

    raise RuntimeError(
        "Could not load a Hallmark library using gseapy.get_library().\n"
        "Attempts:\n"
        + "\n".join(errors)
    )


def standardize_gsea_result(
    result_df: pd.DataFrame,
    cell_type: str,
) -> pd.DataFrame:
    """
    Normalize gseapy result column names across versions.
    """
    out = result_df.copy()

    # Pathway name may be a column or the index.
    if "Term" in out.columns:
        out = out.rename(columns={"Term": "pathway"})
    elif "Name" in out.columns:
        out = out.rename(columns={"Name": "pathway"})
    elif "pathway" not in out.columns:
        out = out.reset_index()
        out = out.rename(columns={out.columns[0]: "pathway"})

    rename_map = {}

    for col in out.columns:
        key = str(col).strip().lower()

        if key == "nes":
            rename_map[col] = "NES"
        elif key in {"fdr q-val", "fdr q-value", "fdr", "fdr_q_val"}:
            rename_map[col] = "FDR"
        elif key in {"nom p-val", "nom p-value", "pval", "p-value"}:
            rename_map[col] = "NOM_P"
        elif key in {"es", "enrichment score"}:
            rename_map[col] = "ES"
        elif key in {
            "lead_genes",
            "lead genes",
            "leadingedge",
            "leading edge",
        }:
            rename_map[col] = "leading_genes"

    out = out.rename(columns=rename_map)

    if "NES" not in out.columns or "FDR" not in out.columns:
        raise ValueError(
            "Could not identify NES/FDR columns in gseapy output. "
            f"Columns: {list(out.columns)}"
        )

    out["NES"] = pd.to_numeric(out["NES"], errors="coerce")
    out["FDR"] = pd.to_numeric(out["FDR"], errors="coerce")

    out["cell_type"] = cell_type
    out["comparison"] = "early_recurrence_vs_recurrence_free"

    return out


def run_prerank(
    ranking: pd.Series,
    gene_sets: dict[str, list[str]],
    cell_type: str,
) -> pd.DataFrame:
    """
    Run Hallmark preranked GSEA using the recovered historical settings.
    """
    ranking = ranking.dropna().copy()

    if ranking.empty:
        raise ValueError(f"No valid ranking statistics for {cell_type}")

    # Historical code guarded against duplicate gene names by retaining
    # the entry with the greatest absolute score.
    if not ranking.index.is_unique:
        tmp = pd.DataFrame(
            {
                "gene": ranking.index.astype(str),
                "score": ranking.to_numpy(),
            }
        )

        tmp["abs_score"] = tmp["score"].abs()

        tmp = (
            tmp.sort_values("abs_score", ascending=False)
            .drop_duplicates("gene", keep="first")
        )

        ranking = tmp.set_index("gene")["score"]

    ranking = ranking.sort_values(ascending=False)

    rank_df = ranking.reset_index()
    rank_df.columns = ["gene", "score"]

    prerank = gp.prerank(
        rnk=rank_df,
        gene_sets=gene_sets,
        processes=N_CPUS,
        permutation_num=GSEA_PERMUTATIONS,
        min_size=GSEA_MIN_SIZE,
        max_size=GSEA_MAX_SIZE,
        seed=GSEA_SEED,
        outdir=None,
        verbose=False,
    )

    return standardize_gsea_result(
        prerank.res2d,
        cell_type=cell_type,
    )


def make_deseq_dataset(
    counts: pd.DataFrame,
    metadata: pd.DataFrame,
):
    """
    Support both the historical and newer PyDESeq2 APIs.
    """
    try:
        return DeseqDataSet(
            counts=counts,
            metadata=metadata,
            design_factors=["analysis_outcome"],
            ref_level=["analysis_outcome", "recurrence_free"],
            n_cpus=N_CPUS,
        )
    except TypeError:
        return DeseqDataSet(
            counts=counts,
            metadata=metadata,
            design="~ analysis_outcome",
            n_cpus=N_CPUS,
        )


def make_deseq_stats(dds):
    """
    Support PyDESeq2 versions with or without n_cpus in DeseqStats().
    """
    kwargs = {
        "dds": dds,
        "contrast": [
            "analysis_outcome",
            "early_recurrence",
            "recurrence_free",
        ],
    }

    try:
        return DeseqStats(
            **kwargs,
            n_cpus=N_CPUS,
        )
    except TypeError:
        return DeseqStats(**kwargs)


# ============================================================
# 3. ARGUMENTS
# ============================================================

def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Pre-BCG sample-level pseudobulk DE and "
            "Hallmark preranked GSEA."
        )
    )

    parser.add_argument(
        "--input-h5ad",
        default="data/raw/RAW_DATA_OBJECT.h5ad",
        help="Annotated AnnData containing raw RNA counts.",
    )

    parser.add_argument(
        "--counts-layer",
        default="X",
        help=(
            "Raw-count source. Use 'X' for adata.X or specify a "
            "layer such as 'counts'."
        ),
    )

    parser.add_argument(
        "--sample-metadata",
        default="metadata/sample_metadata.csv",
    )

    parser.add_argument(
        "--hallmark-gmt",
        default=None,
        help=(
            "Optional local Hallmark GMT. If omitted, the historical "
            "gseapy.get_library() lookup is used."
        ),
    )

    parser.add_argument(
        "--outdir",
        default="results/04_pre_bcg_pathway_enrichment",
    )

    return parser.parse_args()


# ============================================================
# 4. MAIN
# ============================================================

def main():
    args = parse_args()

    input_path = Path(args.input_h5ad)
    sample_metadata_path = Path(args.sample_metadata)
    hallmark_gmt = Path(args.hallmark_gmt) if args.hallmark_gmt else None

    outdir = ensure_dir(Path(args.outdir))
    de_dir = ensure_dir(outdir / "DE")
    gsea_dir = ensure_dir(outdir / "HALLMARK")
    figure_dir = ensure_dir(outdir / "figures")

    if not input_path.exists():
        raise FileNotFoundError(f"Input AnnData not found: {input_path}")

    if not sample_metadata_path.exists():
        raise FileNotFoundError(
            f"Sample metadata not found: {sample_metadata_path}"
        )

    # ========================================================
    # 5. LOAD ANNDATA
    # ========================================================

    print()
    print("=" * 60)
    print("Loading annotated RNA counts")
    print("=" * 60)

    adata = ad.read_h5ad(input_path)

    print(f"Cells: {adata.n_obs:,}")
    print(f"Genes: {adata.n_vars:,}")

    required_obs = ["sample_id", "cell_type"]

    missing_obs = [
        col
        for col in required_obs
        if col not in adata.obs.columns
    ]

    if missing_obs:
        raise ValueError(
            "Input AnnData is missing required obs columns: "
            + ", ".join(missing_obs)
        )

    if not adata.var_names.is_unique:
        raise ValueError(
            "AnnData var_names are not unique. "
            "Resolve duplicate gene names before pseudobulk analysis."
        )

    counts = get_counts_matrix(
        adata,
        args.counts_layer,
    )

    print(f"Raw count source: {args.counts_layer}")

    # ========================================================
    # 6. LOAD CANONICAL SAMPLE METADATA
    # ========================================================

    sample_meta = pd.read_csv(sample_metadata_path)

    required_sample_columns = [
        "sample_number",
        "sample_id",
        "patient_id",
        "timepoint",
        "analysis_outcome",
    ]

    missing_sample_columns = [
        col
        for col in required_sample_columns
        if col not in sample_meta.columns
    ]

    if missing_sample_columns:
        raise ValueError(
            "metadata/sample_metadata.csv is missing required columns: "
            + ", ".join(missing_sample_columns)
        )

    if sample_meta["sample_number"].duplicated().any():
        raise ValueError(
            "Duplicate sample_number values found in sample_metadata.csv."
        )

    if sample_meta.shape[0] != 32:
        raise ValueError(
            f"Expected 32 rows in sample_metadata.csv; found {sample_meta.shape[0]}."
        )

    # ========================================================
    # 7. ATTACH CANONICAL METADATA
    #
    # IMPORTANT:
    # RAW_DATA_OBJECT.h5ad uses numeric obs["sample_id"] values
    # corresponding directly to sample_number (1...32).
    # ========================================================

    obs = adata.obs.copy()

    raw_sample_id = (
        obs["sample_id"]
        .astype(str)
        .str.strip()
    )

    numeric_sample_id = pd.to_numeric(
        raw_sample_id,
        errors="coerce",
    )

    if numeric_sample_id.isna().any():
        bad_ids = sorted(
            raw_sample_id[
                numeric_sample_id.isna()
            ].unique()
        )

        raise ValueError(
            "RAW_DATA_OBJECT.h5ad was expected to contain numeric "
            "sample_id values corresponding to sample_number, but "
            "these values could not be converted to numbers:\n"
            + "\n".join(bad_ids)
        )

    numeric_array = numeric_sample_id.to_numpy(dtype=float)

    if not np.all(
        np.isclose(
            numeric_array,
            np.round(numeric_array),
        )
    ):
        raise ValueError(
            "Some AnnData sample_id values are numeric but not integers."
        )

    obs["sample_number"] = (
        np.round(numeric_array)
        .astype(int)
    )

    sample_map = (
        sample_meta
        .set_index("sample_number")
    )

    missing_sample_numbers = sorted(
        set(obs["sample_number"])
        - set(sample_map.index)
    )

    if missing_sample_numbers:
        raise ValueError(
            "Could not map these AnnData sample numbers to "
            "metadata/sample_metadata.csv:\n"
            + "\n".join(map(str, missing_sample_numbers))
        )

    obs["canonical_sample_id"] = (
        obs["sample_number"]
        .map(sample_map["sample_id"])
    )

    obs["patient_id"] = (
        obs["sample_number"]
        .map(sample_map["patient_id"])
    )

    obs["timepoint"] = (
        obs["sample_number"]
        .map(sample_map["timepoint"])
    )

    obs["analysis_outcome"] = (
        obs["sample_number"]
        .map(sample_map["analysis_outcome"])
    )

    mapping_columns = [
        "canonical_sample_id",
        "patient_id",
        "timepoint",
        "analysis_outcome",
    ]

    if obs[mapping_columns].isna().any().any():
        raise ValueError(
            "Canonical sample metadata mapping produced missing values."
        )

    allowed_outcomes = {
        "early_recurrence",
        "recurrence_free",
    }

    observed_outcomes = set(
        obs["analysis_outcome"].unique()
    )

    if not observed_outcomes.issubset(allowed_outcomes):
        raise ValueError(
            "Unexpected analysis_outcome labels after metadata mapping: "
            + ", ".join(sorted(observed_outcomes))
        )

    sample_mapping_check = (
        obs[
            [
                "sample_number",
                "canonical_sample_id",
                "patient_id",
                "timepoint",
                "analysis_outcome",
            ]
        ]
        .drop_duplicates()
        .sort_values("sample_number")
    )

    if sample_mapping_check.shape[0] != 32:
        raise ValueError(
            "Expected 32 unique mapped AnnData samples; found "
            f"{sample_mapping_check.shape[0]}."
        )

    print()
    print("Sample-number mapping successful:")
    print(
        sample_mapping_check.to_string(
            index=False
        )
    )

    sample_mapping_check.to_csv(
        outdir / "04_sample_mapping_check.tsv",
        sep="\t",
        index=False,
    )

    # ========================================================
    # 8. SUBSET PRE-BCG CELLS
    # ========================================================

    pre_mask = (
        obs["timepoint"]
        .eq("pre")
        .to_numpy()
    )

    obs_pre = obs.loc[pre_mask].copy()
    counts_pre = counts[pre_mask, :]

    print()
    print("=" * 60)
    print("Pre-BCG subset")
    print("=" * 60)

    print(f"Pre-BCG cells: {obs_pre.shape[0]:,}")
    print(
        "Pre-BCG samples: "
        f"{obs_pre['canonical_sample_id'].nunique()}"
    )
    print(
        "Pre-BCG patients: "
        f"{obs_pre['patient_id'].nunique()}"
    )

    pre_sample_meta = (
        obs_pre[
            [
                "canonical_sample_id",
                "patient_id",
                "analysis_outcome",
            ]
        ]
        .drop_duplicates()
        .sort_values("canonical_sample_id")
    )

    print("\nPre-BCG sample outcomes:")
    print(
        pre_sample_meta[
            "analysis_outcome"
        ].value_counts()
    )

    if obs_pre["canonical_sample_id"].nunique() != 16:
        raise ValueError(
            "Expected 16 pre-BCG samples after canonical mapping."
        )

    if obs_pre["patient_id"].nunique() != 16:
        raise ValueError(
            "Expected 16 pre-BCG patients after canonical mapping."
        )

    outcome_sample_counts = (
        pre_sample_meta["analysis_outcome"]
        .value_counts()
        .to_dict()
    )

    if outcome_sample_counts.get("early_recurrence", 0) != 8:
        raise ValueError(
            "Expected 8 early-recurrence pre-BCG samples."
        )

    if outcome_sample_counts.get("recurrence_free", 0) != 8:
        raise ValueError(
            "Expected 8 recurrence-free pre-BCG samples."
        )

    # ========================================================
    # 9. FILTER SAMPLE × CELL-TYPE GROUPS
    #
    # Historical setting:
    #   MIN_CELLS_PER_SAMPLE_CT = 1
    # ========================================================

    group_size = (
        obs_pre
        .groupby(
            [
                "canonical_sample_id",
                "cell_type",
            ],
            observed=True,
        )
        .size()
        .rename("n_cells")
    )

    keep_groups = group_size[
        group_size >= MIN_CELLS_PER_SAMPLE_CT
    ].index

    cell_group_index = pd.MultiIndex.from_frame(
        obs_pre[
            [
                "canonical_sample_id",
                "cell_type",
            ]
        ]
    )

    keep_cells = cell_group_index.isin(keep_groups)

    obs_pre = obs_pre.loc[keep_cells].copy()
    counts_pre = counts_pre[keep_cells, :]

    group_size.reset_index().to_csv(
        outdir / "04_pre_bcg_cells_per_sample_celltype.tsv",
        sep="\t",
        index=False,
    )

    print(
        "\nCells retained after sample×cell-type filter: "
        f"{obs_pre.shape[0]:,}"
    )

    # ========================================================
    # 10. BUILD SAMPLE × CELL-TYPE PSEUDOBULKS
    # ========================================================

    print()
    print("=" * 60)
    print("Constructing pseudobulk profiles")
    print("=" * 60)

    group_keys = pd.MultiIndex.from_frame(
        obs_pre[
            [
                "canonical_sample_id",
                "cell_type",
            ]
        ],
        names=["sample_id", "cell_type"],
    )

    codes, unique_groups = pd.factorize(group_keys)

    order = np.argsort(codes)
    codes_sorted = codes[order]
    counts_sorted = counts_pre[order, :]

    boundaries = np.flatnonzero(
        np.r_[
            True,
            codes_sorted[1:] != codes_sorted[:-1],
            True,
        ]
    )

    pseudobulk_rows = []

    for i in range(len(boundaries) - 1):
        start = boundaries[i]
        stop = boundaries[i + 1]

        pseudobulk_rows.append(
            sparse.csr_matrix(
                counts_sorted[start:stop, :].sum(axis=0)
            )
        )

    pb_matrix = sparse.vstack(
        pseudobulk_rows
    ).tocsr()

    pb_index = pd.MultiIndex.from_tuples(
        list(unique_groups),
        names=["sample_id", "cell_type"],
    )

    pb_counts = pd.DataFrame.sparse.from_spmatrix(
        pb_matrix,
        index=pb_index,
        columns=adata.var_names.astype(str),
    )

    # --------------------------------------------------------
    # Construct pseudobulk metadata directly from the
    # pseudobulk count-matrix index.
    #
    # This guarantees identical row order between pb_counts
    # and pb_meta.
    # --------------------------------------------------------

    pb_meta = (
        pb_counts.index
        .to_frame(index=False)
    )

    # --------------------------------------------------------
    # Sample-level canonical metadata
    # --------------------------------------------------------

    sample_level_meta = (
        obs_pre[
            [
                "canonical_sample_id",
                "patient_id",
                "analysis_outcome",
            ]
        ]
        .drop_duplicates()
        .rename(
            columns={
                "canonical_sample_id": "sample_id"
            }
        )
    )

    # Each pre-BCG sample must map to exactly one
    # patient/outcome record.
    if sample_level_meta["sample_id"].duplicated().any():

        duplicated_samples = (
            sample_level_meta.loc[
                sample_level_meta[
                    "sample_id"
                ].duplicated(
                    keep=False
                ),
                "sample_id",
            ]
            .unique()
            .tolist()
        )

        raise RuntimeError(
            "Some samples map to multiple "
            "patient/outcome records:\n"
            + "\n".join(
                map(
                    str,
                    duplicated_samples,
                )
            )
        )

    # --------------------------------------------------------
    # Attach patient/outcome metadata to each pseudobulk row
    # --------------------------------------------------------

    pb_meta = pb_meta.merge(
        sample_level_meta,
        on="sample_id",
        how="left",
        validate="many_to_one",
    )

    # Restore exactly the same MultiIndex as pb_counts.
    pb_meta = pb_meta.set_index(
        [
            "sample_id",
            "cell_type",
        ]
    )

    # --------------------------------------------------------
    # Validate metadata
    # --------------------------------------------------------

    if pb_meta[
        [
            "patient_id",
            "analysis_outcome",
        ]
    ].isna().any().any():

        missing_rows = pb_meta.loc[
            pb_meta[
                [
                    "patient_id",
                    "analysis_outcome",
                ]
            ]
            .isna()
            .any(axis=1)
        ]

        raise RuntimeError(
            "Missing metadata for one or more "
            "pseudobulk rows:\n"
            + missing_rows.to_string()
        )

    if not pb_counts.index.equals(
        pb_meta.index
    ):
        raise RuntimeError(
            "Pseudobulk count matrix and metadata "
            "are still not aligned."
        )

    print(
        "Pseudobulk count/metadata alignment: OK"
    )

    print(
        f"Pseudobulk profiles: {pb_counts.shape[0]}"
    )

    print(
        f"Genes before filter: {pb_counts.shape[1]}"
    )

    # ========================================================
    # 11. GLOBAL GENE FILTER
    #
    # Historical settings:
    #   total counts >= 10
    #   non-zero in >= 3 pseudobulk profiles
    # ========================================================

    gene_totals = np.asarray(
        pb_matrix.sum(axis=0)
    ).ravel()

    gene_nonzero_samples = np.asarray(
        (pb_matrix > 0).sum(axis=0)
    ).ravel()

    keep_gene = (
        (gene_totals >= MIN_TOTAL_COUNTS_GENE)
        &
        (gene_nonzero_samples >= MIN_SAMPLES_GENE)
    )

    if keep_gene.sum() == 0:
        raise ValueError(
            "Global gene filtering removed every gene."
        )

    pb_counts = pb_counts.loc[
        :,
        pb_counts.columns[keep_gene],
    ]

    print(
        f"Genes after filter: {pb_counts.shape[1]}"
    )

    # Save pseudobulk inputs.
    pb_counts.to_csv(
        outdir / "04_pre_bcg_pseudobulk_counts.tsv.gz",
        sep="\t",
        compression="gzip",
    )

    pb_meta.to_csv(
        outdir / "04_pre_bcg_pseudobulk_metadata.tsv",
        sep="\t",
    )

    # ========================================================
    # 12. LOAD HALLMARK GENE SETS
    # ========================================================

    hallmark_library, hallmark_source = build_hallmark_library(
        hallmark_gmt
    )

    print()
    print("=" * 60)
    print("Hallmark gene sets")
    print("=" * 60)

    print(f"Source: {hallmark_source}")
    print(f"Number of gene sets: {len(hallmark_library)}")

    # ========================================================
    # 13. PSEUDOBULK DE + HALLMARK GSEA BY CELL TYPE
    #
    # Cross-sectional pre-BCG model:
    #   ~ analysis_outcome
    #
    # Contrast:
    #   early_recurrence vs recurrence_free
    # ========================================================

    cell_types = sorted(
        pb_counts.index
        .get_level_values("cell_type")
        .unique()
        .astype(str)
    )

    de_tables = []
    gsea_tables = []
    summary_rows = []

    print()
    print("=" * 60)
    print("Running pre-BCG pseudobulk DE + Hallmark GSEA")
    print("=" * 60)

    for cell_type in cell_types:
        print()
        print("-" * 60)
        print(f"Cell type: {cell_type}")

        counts_ct = pb_counts.xs(
            cell_type,
            level="cell_type",
        )

        metadata_ct = (
            pb_meta.xs(
                cell_type,
                level="cell_type",
            )
            .copy()
            .loc[counts_ct.index]
        )

        group_counts = (
            metadata_ct["analysis_outcome"]
            .value_counts()
        )

        n_early = int(
            group_counts.get(
                "early_recurrence",
                0,
            )
        )

        n_free = int(
            group_counts.get(
                "recurrence_free",
                0,
            )
        )

        print(
            f"  early_recurrence samples: {n_early}"
        )
        print(
            f"  recurrence_free samples: {n_free}"
        )

        if n_early < 2 or n_free < 2:
            print(
                "  SKIP: fewer than two samples in one outcome group."
            )

            summary_rows.append(
                {
                    "cell_type": cell_type,
                    "n_early_recurrence": n_early,
                    "n_recurrence_free": n_free,
                    "n_genes_tested": np.nan,
                    "n_significant_hallmark": np.nan,
                    "status": "SKIP_TOO_FEW_SAMPLES",
                }
            )

            continue

        counts_int = (
            counts_ct.sparse
            .to_dense()
            .round()
            .astype(np.int32)
        )

        model_metadata = pd.DataFrame(
            {
                "analysis_outcome":
                    metadata_ct[
                        "analysis_outcome"
                    ].astype(str)
            },
            index=metadata_ct.index,
        )

        model_metadata["analysis_outcome"] = pd.Categorical(
            model_metadata["analysis_outcome"],
            categories=[
                "recurrence_free",
                "early_recurrence",
            ],
            ordered=True,
        )

        # ----------------------------------------------------
        # PyDESeq2
        # ----------------------------------------------------

        dds = make_deseq_dataset(
            counts=counts_int,
            metadata=model_metadata,
        )

        dds.deseq2()

        stats = make_deseq_stats(dds)
        stats.summary()

        result = stats.results_df.copy()

        result.index = result.index.astype(str)
        result.index.name = "gene"

        result["cell_type"] = cell_type
        result["comparison"] = (
            "early_recurrence_vs_recurrence_free"
        )

        result["direction"] = np.where(
            result["log2FoldChange"] > 0,
            "higher_in_early_recurrence",
            "higher_in_recurrence_free",
        )

        de_tables.append(
            result.reset_index()
        )

        result.to_csv(
            de_dir
            / (
                safe_name(cell_type)
                + "_early_recurrence_vs_recurrence_free.tsv"
            ),
            sep="\t",
        )

        # ----------------------------------------------------
        # GSEA ranking = DESeq2 Wald statistic
        # ----------------------------------------------------

        if "stat" not in result.columns:
            raise ValueError(
                f"DESeq2 result for {cell_type} does not contain 'stat'."
            )

        ranking = result["stat"].copy()
        ranking.index = result.index

        gsea_result = run_prerank(
            ranking=ranking,
            gene_sets=hallmark_library,
            cell_type=cell_type,
        )

        gsea_result["significant"] = (
            gsea_result["FDR"] < FDR_THRESHOLD
        )

        gsea_result["direction"] = np.where(
            gsea_result["NES"] > 0,
            "enriched_in_early_recurrence",
            "enriched_in_recurrence_free",
        )

        gsea_tables.append(gsea_result)

        gsea_result.to_csv(
            gsea_dir
            / (
                safe_name(cell_type)
                + "_HALLMARK_gsea.tsv"
            ),
            sep="\t",
            index=False,
        )

        n_significant = int(
            gsea_result["significant"].sum()
        )

        print(
            "  Significant Hallmark pathways "
            f"(FDR < {FDR_THRESHOLD}): {n_significant}"
        )

        summary_rows.append(
            {
                "cell_type": cell_type,
                "n_early_recurrence": n_early,
                "n_recurrence_free": n_free,
                "n_genes_tested": int(
                    result["stat"].notna().sum()
                ),
                "n_significant_hallmark": n_significant,
                "status": "OK",
            }
        )

    # ========================================================
    # 14. COMBINE RESULTS
    # ========================================================

    if not de_tables:
        raise RuntimeError(
            "No cell type produced a pseudobulk DE result."
        )

    if not gsea_tables:
        raise RuntimeError(
            "No cell type produced a Hallmark GSEA result."
        )

    all_de = pd.concat(
        de_tables,
        ignore_index=True,
    )

    all_gsea = pd.concat(
        gsea_tables,
        ignore_index=True,
    )

    analysis_summary = (
        pd.DataFrame(summary_rows)
        .sort_values("cell_type")
    )

    significant_gsea = (
        all_gsea.loc[
            all_gsea["FDR"] < FDR_THRESHOLD
        ]
        .sort_values(
            [
                "cell_type",
                "FDR",
                "NES",
            ],
            ascending=[
                True,
                True,
                False,
            ],
        )
    )

    all_de.to_csv(
        outdir / "04_pre_bcg_all_pseudobulk_DE.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )

    all_gsea.to_csv(
        outdir / "04_pre_bcg_all_HALLMARK_gsea.tsv",
        sep="\t",
        index=False,
    )

    significant_gsea.to_csv(
        outdir / "04_pre_bcg_significant_HALLMARK_gsea.tsv",
        sep="\t",
        index=False,
    )

    analysis_summary.to_csv(
        outdir / "04_pre_bcg_analysis_summary.tsv",
        sep="\t",
        index=False,
    )


if __name__ == "__main__":
    main()