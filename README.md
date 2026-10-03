# NMIBC-BCG-Multiome

Reproducibility repository for the single-nucleus RNA-seq and ATAC-seq analyses of peripheral blood mononuclear cells from patients with high-risk non-muscle-invasive bladder cancer (NMIBC) treated with intravesical BCG.

This repository contains the reviewer-facing analysis code, canonical sample metadata, and environment specifications used to reproduce the principal computational analyses reported in the manuscript. No custom software was developed; analyses use established open-source R and Python packages.

## Scope

The repository begins from processed, annotated single-cell/multiome objects rather than raw FASTQ files. It reproduces the main pre-BCG analyses used to compare patients who experienced early recurrence with those who remained recurrence-free, including:

| Script | Analysis |
|---|---|
| `03_pre_bcg_de.R` | Cell-level differential gene expression by annotated cell type |
| `04_pre_bcg_pathway_enrichment.py` | Sample-level pseudobulk differential expression and Hallmark GSEA |
| `05_pre_bcg_atac.R` | Differential chromatin accessibility, motif analysis, and chromVAR |
| `06_cellchat_pre_bcg.R` | Cell-cell communication analysis with CellChat |
| `07_milor_pre_activated_monocytes.R` | Activated-monocyte MiloR differential abundance with pathway and AP-1 follow-up |
| `08_ap1_signature_pre_bcg.R` | Six-gene AP-1 monocyte signature, patient-level comparison, and ROC analysis |

The repository is intended to reproduce manuscript-level statistical analyses from the finalized processed objects. Upstream sequencing alignment, initial quality control, and object construction are described in the manuscript Methods and are not rerun here.

## Repository structure

```text
NMIBC-BCG-Multiome/
├── scripts/          # Analysis scripts
├── metadata/         # Canonical sample metadata and sample-ID mappings
├── environment/      # R/Python package requirements and installation files
├── results/          # Generated tables, figures, and intermediate outputs
└── README.md
```

Canonical metadata are stored in:

- `metadata/sample_metadata.csv`
- `metadata/sample_aliases.csv`

## Required data

The primary processed inputs are:

- `merged_obj.rds` — finalized pre-BCG Seurat/Signac object
- `RAW_DATA_OBJECT.h5ad` — annotated single-nucleus RNA object containing raw counts for pseudobulk analysis

The MiloR workflow additionally uses:

- `post_BCG_merged_obj.rds` — finalized post-BCG Seurat/Signac object

Historical shared-ATAC checkpoint matrices (`pre_mat_81657.rds` and `post_mat_81657.rds`) can be used to avoid rebuilding the shared peak matrix. If these checkpoints are not available, the MiloR workflow requires access to the corresponding indexed ATAC fragment files.

Large processed objects are not stored directly in this GitHub repository. Their availability should be interpreted together with the manuscript Data Availability statement.

## Software environment

The analyses were developed using:

- **R 4.4.x**
- **Python 3.10.x**
- Seurat / Signac
- Scanpy
- Harmony / harmonypy
- PyDESeq2
- gseapy
- CellChat
- MiloR
- chromVAR
- JASPAR2020

Python dependencies are listed in `environment/requirements.txt`. R/Bioconductor dependencies are documented in the environment setup files.

For full runs, a machine with approximately **64 GB RAM or more** is recommended. Eight CPU cores are a practical starting point. A GPU is not required.

## Running the analyses

Run commands from the repository root.

```bash
Rscript scripts/03_pre_bcg_de.R
python scripts/04_pre_bcg_pathway_enrichment.py
Rscript scripts/05_pre_bcg_atac.R
Rscript scripts/06_cellchat_pre_bcg.R
Rscript scripts/07_milor_pre_activated_monocytes.R
Rscript scripts/08_ap1_signature_pre_bcg.R
```

Scripts write their outputs beneath `results/`. Individual scripts also save analysis-specific parameter summaries and/or session information where applicable.

## Computational reproducibility

Random seeds, Harmony settings, multimodal integration parameters, differential-analysis thresholds, and other analysis-specific settings are documented in **Supplemental Table 3: Computational parameters and random seeds used for single-nucleus multiomic analyses**.

All Harmony integrations used the **harmonypy v0.0.9 package defaults**; no Harmony integration parameter was manually overridden. Historical weighted-nearest-neighbour settings and UMAP seeds were recovered from the finalized analysis object where available.

Key stochastic settings include:

- joint WNN UMAP seed: **42**
- MiloR neighbourhood sampling seed: **1**
- activated-monocyte neighbourhood k-means seed: **1**
- Hallmark GSEA seed: **0**
- AP-1 `AddModuleScore` seed: **1**

## Output interpretation

For outcome comparisons, positive effect sizes are generally defined in the direction of **early recurrence relative to recurrence-free patients** unless otherwise stated in the script or output table. Each script contains a header describing its comparison, directionality, and principal parameters.

## Code availability

No custom software was developed for this study. This repository provides the analysis scripts, metadata mappings, software requirements, and reproducibility parameters needed to rerun the manuscript analyses from the processed data objects.
