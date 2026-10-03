# ============================================================
# 02_joint_embedding.R
#
# Reproduce the joint RNA-ATAC embedding used for the
# pre-BCG and post-BCG UMAP panels.
#
# IMPORTANT:
# This script computes ONE joint embedding across the full
# 198,652-cell multiome object, then generates separate
# pre-BCG and post-BCG figures by subsetting cells from the
# same embedding.
#
# Input:
#   TFRI_P3_signac.rds
#
# Run from repository root:
#
#   Rscript scripts/02_joint_embedding.R "G:/TFRI_P3_signac.rds"
#
# ============================================================


# ============================================================
# 1. COMMAND-LINE ARGUMENT
# ============================================================

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 1) {
  stop(
    paste0(
      "\nUsage:\n",
      "Rscript scripts/02_joint_embedding.R ",
      "\"path/to/TFRI_P3_signac.rds\"\n"
    )
  )
}

input_rds <- args[1]

if (!file.exists(input_rds)) {
  stop("Input object not found: ", input_rds)
}


# ============================================================
# 2. REQUIRED PACKAGES
# ============================================================

required_packages <- c(
  "Seurat",
  "Signac",
  "harmony",
  "ggplot2"
)

missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    quietly = TRUE,
    FUN.VALUE = logical(1)
  )
]

if (length(missing_packages) > 0) {
  stop(
    "Missing required package(s): ",
    paste(missing_packages, collapse = ", ")
  )
}

suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(harmony)
  library(ggplot2)
})


# ============================================================
# 3. FILE PATHS
# ============================================================

sample_metadata_path <- file.path(
  "metadata",
  "sample_metadata.csv"
)

sample_aliases_path <- file.path(
  "metadata",
  "sample_aliases.csv"
)

results_table_dir <- file.path(
  "results",
  "tables"
)

results_figure_dir <- file.path(
  "results",
  "figures"
)

dir.create(
  results_table_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  results_figure_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ============================================================
# 4. LOAD CANONICAL METADATA
# ============================================================

sample_meta <- read.csv(
  sample_metadata_path,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

sample_aliases <- read.csv(
  sample_aliases_path,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

sample_map <- merge(
  sample_aliases,
  sample_meta,
  by = c(
    "sample_number",
    "sample_name",
    "sample_id"
  ),
  all.x = TRUE,
  all.y = FALSE
)

if (nrow(sample_map) != 32) {
  stop(
    "Canonical sample mapping did not resolve to exactly 32 samples."
  )
}

if (anyDuplicated(sample_map$object_sample_id)) {
  stop(
    "Duplicate object_sample_id values found in sample mapping."
  )
}


# ============================================================
# 5. LOAD FULL MULTIOME OBJECT
# ============================================================

cat("\nLoading full multiome object...\n")

obj <- readRDS(input_rds)

if (!inherits(obj, "Seurat")) {
  stop("Input object is not a Seurat object.")
}

cat(
  "Cells loaded: ",
  ncol(obj),
  "\n",
  sep = ""
)

cat(
  "Features: ",
  nrow(obj),
  "\n",
  sep = ""
)

cat(
  "Assays: ",
  paste(Assays(obj), collapse = ", "),
  "\n",
  sep = ""
)


# ============================================================
# 6. VALIDATE EXPECTED INPUT STATE
# ============================================================

if (ncol(obj) != 198652) {
  stop(
    "Expected 198,652 cells; found ",
    ncol(obj),
    "."
  )
}

required_assays <- c(
  "RNA",
  "ATAC"
)

missing_assays <- setdiff(
  required_assays,
  Assays(obj)
)

if (length(missing_assays) > 0) {
  stop(
    "Missing required assay(s): ",
    paste(missing_assays, collapse = ", ")
  )
}

if (!"lsi" %in% Reductions(obj)) {
  stop(
    "ATAC LSI reduction ('lsi') was not found in the input object."
  )
}

if (!"sample_id" %in% colnames(obj@meta.data)) {
  stop(
    "sample_id metadata column was not found."
  )
}

if (!"cell_type" %in% colnames(obj@meta.data)) {
  stop(
    "cell_type metadata column was not found."
  )
}


# ============================================================
# 7. ATTACH CANONICAL METADATA
# ============================================================

mapping_index <- match(
  as.character(obj$sample_id),
  sample_map$object_sample_id
)

if (any(is.na(mapping_index))) {

  unresolved_ids <- unique(
    as.character(obj$sample_id)[
      is.na(mapping_index)
    ]
  )

  stop(
    paste0(
      "Unresolved object sample IDs:\n",
      paste(unresolved_ids, collapse = "\n")
    )
  )
}

obj$canonical_sample_number <-
  sample_map$sample_number[
    mapping_index
  ]

obj$canonical_sample_name <-
  sample_map$sample_name[
    mapping_index
  ]

obj$canonical_sample_id <-
  sample_map$sample_id[
    mapping_index
  ]

obj$patient_id <-
  sample_map$patient_id[
    mapping_index
  ]

obj$timepoint <-
  sample_map$timepoint[
    mapping_index
  ]

obj$analysis_outcome <-
  sample_map$analysis_outcome[
    mapping_index
  ]

cat("\nCanonical metadata mapping successful.\n")

cat("\nCells by timepoint:\n")
print(table(obj$timepoint))

cat("\nCells by outcome:\n")
print(table(obj$analysis_outcome))

if (sum(obj$timepoint == "pre") != 100523) {
  stop("Expected 100,523 pre-BCG cells.")
}

if (sum(obj$timepoint == "post") != 98129) {
  stop("Expected 98,129 post-BCG cells.")
}

if (length(unique(obj$canonical_sample_id)) != 32) {
  stop("Expected 32 samples.")
}

if (length(unique(obj$patient_id)) != 16) {
  stop("Expected 16 patients.")
}


# ============================================================
# 8. RNA PROCESSING
#
# Reconstruct RNA representation on the full joint object.
# ============================================================

cat("\n============================================\n")
cat("RNA processing\n")
cat("============================================\n")

DefaultAssay(obj) <- "RNA"

set.seed(42)

obj <- NormalizeData(
  obj,
  verbose = FALSE
)

obj <- FindVariableFeatures(
  obj,
  verbose = FALSE
)

obj <- ScaleData(
  obj,
  verbose = FALSE
)

obj <- RunPCA(
  obj,
  npcs = 50,
  reduction.name = "pca",
  verbose = FALSE
)


# ------------------------------------------------------------
# Harmony on RNA PCA
# ------------------------------------------------------------

cat("\nRunning Harmony on RNA PCA...\n")

set.seed(42)

obj <- RunHarmony(
  object = obj,
  group.by.vars = "sample_id",
  reduction = "pca",
  assay.use = "RNA",
  reduction.save = "harmony_rna",
  verbose = FALSE
)


# ============================================================
# 9. ATAC PROCESSING
#
# Use the existing LSI reduction, then run Harmony.
# ============================================================

cat("\n============================================\n")
cat("ATAC processing\n")
cat("============================================\n")

DefaultAssay(obj) <- "ATAC"

cat("\nRunning Harmony on ATAC LSI...\n")

set.seed(42)

obj <- RunHarmony(
  object = obj,
  group.by.vars = "sample_id",
  reduction = "lsi",
  assay.use = "ATAC",
  reduction.save = "harmony_atac",
  project.dim = FALSE,
  verbose = FALSE
)


# ============================================================
# 10. JOINT RNA-ATAC INTEGRATION
#
# This recreates one joint embedding across ALL cells.
# ============================================================

cat("\n============================================\n")
cat("Weighted nearest-neighbor integration\n")
cat("============================================\n")

obj <- FindMultiModalNeighbors(
  object = obj,
  reduction.list = list(
    "harmony_rna",
    "harmony_atac"
  ),
  dims.list = list(
    1:30,
    2:30
  ),
  modality.weight.name = "RNA_ATAC.weight",
  verbose = TRUE
)


# ============================================================
# 11. JOINT UMAP
#
# Compute one UMAP across ALL cells.
# ============================================================

cat("\n============================================\n")
cat("Running joint WNN UMAP\n")
cat("============================================\n")

set.seed(42)

obj <- RunUMAP(
  object = obj,
  nn.name = "weighted.nn",
  reduction.name = "wnn.umap",
  reduction.key = "wnnUMAP_",
  seed.use = 42,
  verbose = TRUE
)


# ============================================================
# 12. EXPORT UMAP COORDINATES
# ============================================================

umap_coordinates <- as.data.frame(
  Embeddings(
    obj,
    reduction = "wnn.umap"
  )
)

if (ncol(umap_coordinates) != 2) {
  stop("WNN UMAP does not contain two dimensions.")
}

colnames(umap_coordinates) <- c(
  "UMAP_1",
  "UMAP_2"
)

umap_coordinates$cell_barcode <- rownames(umap_coordinates)
umap_coordinates$sample_id <- obj$canonical_sample_id
umap_coordinates$patient_id <- obj$patient_id
umap_coordinates$timepoint <- obj$timepoint
umap_coordinates$analysis_outcome <- obj$analysis_outcome
umap_coordinates$cell_type <- obj$cell_type

write.csv(
  umap_coordinates,
  file.path(
    results_table_dir,
    "joint_umap_coordinates.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 13. EXPORT CELL COUNTS
# ============================================================

cell_counts <- as.data.frame(
  table(
    sample_id = obj$canonical_sample_id,
    patient_id = obj$patient_id,
    timepoint = obj$timepoint,
    outcome = obj$analysis_outcome,
    cell_type = obj$cell_type
  )
)

cell_counts <- cell_counts[
  cell_counts$Freq > 0,
]

write.csv(
  cell_counts,
  file.path(
    results_table_dir,
    "cell_counts_by_sample_and_cell_type.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 14. GENERATE PRE-BCG AND POST-BCG FIGURES
#
# IMPORTANT:
# Both plots come from the SAME joint embedding.
# ============================================================

pre_cells <- rownames(obj@meta.data)[obj$timepoint == "pre"]
post_cells <- rownames(obj@meta.data)[obj$timepoint == "post"]


# ------------------------------------------------------------
# Figure 1a - pre-BCG cells only
# ------------------------------------------------------------

fig1a <- DimPlot(
  object = obj,
  reduction = "wnn.umap",
  group.by = "cell_type",
  cells = pre_cells,
  raster = TRUE,
  pt.size = 0.1
) +
  ggtitle("Pre-BCG") +
  theme_classic() +
  theme(
    legend.title = element_blank()
  )

ggsave(
  filename = file.path(
    results_figure_dir,
    "fig1a_pre_bcg_joint_umap.pdf"
  ),
  plot = fig1a,
  width = 8,
  height = 6,
  units = "in"
)


# ------------------------------------------------------------
# Figure 2c - post-BCG cells only
# ------------------------------------------------------------

fig2c <- DimPlot(
  object = obj,
  reduction = "wnn.umap",
  group.by = "cell_type",
  cells = post_cells,
  raster = TRUE,
  pt.size = 0.1
) +
  ggtitle("Pre-sixth BCG") +
  theme_classic() +
  theme(
    legend.title = element_blank()
  )

ggsave(
  filename = file.path(
    results_figure_dir,
    "fig2c_post_bcg_joint_umap.pdf"
  ),
  plot = fig2c,
  width = 8,
  height = 6,
  units = "in"
)


# ============================================================
# 15. SAVE ANALYSIS PARAMETERS
# ============================================================

analysis_parameters <- data.frame(
  parameter = c(
    "input_cells",
    "RNA_assay",
    "RNA_PCs",
    "RNA_Harmony_reduction",
    "ATAC_reduction_input",
    "ATAC_Harmony_reduction",
    "WNN_RNA_dimensions",
    "WNN_ATAC_dimensions",
    "WNN_modality_weight_name",
    "UMAP_seed"
  ),
  value = c(
    ncol(obj),
    "RNA",
    "1:50",
    "harmony_rna",
    "lsi",
    "harmony_atac",
    "1:30",
    "2:30",
    "RNA_ATAC.weight",
    "42"
  ),
  stringsAsFactors = FALSE
)

write.csv(
  analysis_parameters,
  file.path(
    results_table_dir,
    "02_joint_embedding_parameters.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 16. SESSION INFORMATION
# ============================================================

capture.output(
  sessionInfo(),
  file = file.path(
    results_table_dir,
    "session_info_02_joint_embedding.txt"
  )
)


# ============================================================
# 17. FINAL REPORT
# ============================================================

cat("\n")
cat("============================================\n")
cat("Joint embedding analysis complete\n")
cat("============================================\n\n")

cat(
  "Total cells analysed: ",
  ncol(obj),
  "\n",
  sep = ""
)

cat(
  "Pre-BCG cells: ",
  sum(obj$timepoint == "pre"),
  "\n",
  sep = ""
)

cat(
  "Post-BCG cells: ",
  sum(obj$timepoint == "post"),
  "\n",
  sep = ""
)

cat(
  "Patients: ",
  length(unique(obj$patient_id)),
  "\n",
  sep = ""
)

cat(
  "Samples: ",
  length(unique(obj$canonical_sample_id)),
  "\n",
  sep = ""
)

cat("\nOne joint embedding was computed across all cells.\n")
cat("Pre-BCG and post-BCG figures were generated by subsetting\n")
cat("cells from the same joint WNN UMAP.\n")

cat(
  "\nOutputs written to:\n  ",
  results_table_dir,
  "\n  ",
  results_figure_dir,
  "\n",
  sep = ""
)

cat("\nAnalysis complete.\n")