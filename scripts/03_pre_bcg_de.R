# ============================================================
# 03_pre_bcg_de.R
#
# Reproduce the pre-BCG cell-level differential-expression
# analysis used for Figure 1b.
#
#
# Positive avg_log2FC:
#   higher expression in early-recurrence patients
#
# Input:
#   finalized pre-BCG Seurat object (historically merged_obj.rds)
#
# Run from repository root:
#
#   & "C:\Users\andrew\AppData\Local\Programs\R\R-4.4.1\bin\Rscript.exe" `
#       scripts\03_pre_bcg_de.R "G:/merged_obj.rds"
#
# ============================================================


# ============================================================
# 1. COMMAND-LINE ARGUMENT
# ============================================================

args <- commandArgs(
  trailingOnly = TRUE
)

if (length(args) > 1) {
  stop(
    paste0(
      "\nUsage:\n",
      "Rscript scripts/03_pre_bcg_de.R ",
      "\"path/to/merged_obj.rds\"\n"
    )
  )
}

input_rds <- if (
  length(args) == 1
) {
  args[1]
} else {
  "G:/merged_obj.rds"
}

if (!file.exists(input_rds)) {
  stop(
    "Input object not found: ",
    input_rds
  )
}


# ============================================================
# 2. REQUIRED PACKAGES
# ============================================================

required_packages <- c(
  "Seurat",
  "dplyr",
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
    paste(
      missing_packages,
      collapse = ", "
    )
  )
}

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
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

table_dir <- file.path(
  "results",
  "tables",
  "03_pre_bcg_de"
)

figure_dir <- file.path(
  "results",
  "figures"
)

dir.create(
  table_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  figure_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ============================================================
# 4. ANALYSIS PARAMETERS
#
# These reproduce the historical Figure 1b code recovered
# from TFRI P3 1.Rmd.
# ============================================================

MIN_CELLS_PER_GROUP <- 25
MIN_PCT <- 0.10
LOGFC_THRESHOLD <- 0.25
PADJ_THRESHOLD <- 0.05

NORMALIZATION_METHOD <- "LogNormalize"
NORMALIZATION_SCALE_FACTOR <- 10000


# ============================================================
# 5. LOAD CANONICAL SAMPLE METADATA
# ============================================================

if (!file.exists(sample_metadata_path)) {
  stop(
    "Missing metadata file: ",
    sample_metadata_path
  )
}

if (!file.exists(sample_aliases_path)) {
  stop(
    "Missing metadata file: ",
    sample_aliases_path
  )
}

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
    "Canonical sample mapping did not resolve to 32 samples."
  )
}

if (anyDuplicated(
  sample_map$object_sample_id
)) {
  stop(
    "Duplicate object_sample_id values were found in sample mapping."
  )
}


# ============================================================
# 6. LOAD FINALIZED PRE-BCG OBJECT
# ============================================================

cat("\n")
cat("============================================\n")
cat("Loading finalized pre-BCG object\n")
cat("============================================\n")

obj <- readRDS(
  input_rds
)

if (!inherits(
  obj,
  "Seurat"
)) {
  stop(
    "Input file is not a Seurat object."
  )
}

cat(
  "Input: ",
  normalizePath(
    input_rds,
    winslash = "/",
    mustWork = TRUE
  ),
  "\n",
  sep = ""
)

cat(
  "Cells: ",
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

if (ncol(obj) != 100523) {
  warning(
    "Historical pre-BCG object contained 100,523 cells; found ",
    ncol(obj),
    "."
  )
}


# ============================================================
# 7. VALIDATE REQUIRED OBJECT CONTENT
# ============================================================

if (!"RNA" %in% Assays(obj)) {
  stop(
    "RNA assay is missing from the input object."
  )
}

if (!"sample_id" %in% colnames(
  obj@meta.data
)) {
  stop(
    "sample_id is missing from object metadata."
  )
}

if (!"cell_type" %in% colnames(
  obj@meta.data
)) {
  stop(
    "cell_type is missing from object metadata."
  )
}

if (any(
  is.na(obj$cell_type)
)) {
  warning(
    "Some cells have missing cell_type annotations."
  )
}


# ============================================================
# 8. ATTACH CANONICAL METADATA
#
# The original code used the historical Outcome labels
# "recurrence" and "remission".
#
# Here we derive those labels from the canonical metadata so
# the analysis does not depend on inconsistent legacy Outcome
# capitalization in the stored object.
# ============================================================

mapping_index <- match(
  as.character(
    obj$sample_id
  ),
  sample_map$object_sample_id
)

if (any(
  is.na(mapping_index)
)) {

  unresolved <- unique(
    as.character(
      obj$sample_id
    )[
      is.na(mapping_index)
    ]
  )

  stop(
    paste0(
      "Unresolved object sample IDs:\n",
      paste(
        unresolved,
        collapse = "\n"
      )
    )
  )
}


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


# ------------------------------------------------------------
# Confirm this is exclusively the pre-BCG object
# ------------------------------------------------------------

if (!all(
  obj$timepoint == "pre"
)) {
  stop(
    "Input object contains non-pre-BCG cells."
  )
}

if (
  length(
    unique(
      obj$canonical_sample_id
    )
  ) != 16
) {
  stop(
    "Expected 16 pre-BCG samples."
  )
}

if (
  length(
    unique(
      obj$patient_id
    )
  ) != 16
) {
  stop(
    "Expected 16 pre-BCG patients."
  )
}


# ------------------------------------------------------------
# Recreate historical outcome names used in original code
# ------------------------------------------------------------

obj$de_outcome <- ifelse(
  obj$analysis_outcome ==
    "early_recurrence",
  "recurrence",
  ifelse(
    obj$analysis_outcome ==
      "recurrence_free",
    "remission",
    NA_character_
  )
)

if (any(
  is.na(obj$de_outcome)
)) {
  stop(
    "Some cells could not be assigned to recurrence/remission."
  )
}

cat("\nOutcome counts:\n")

print(
  table(
    obj$de_outcome
  )
)


# ============================================================
# 9. PREPARE RNA ASSAY
#
# This reproduces the historical code:
#
#   DefaultAssay(obj) <- "RNA"
#   obj[["RNA"]] <- JoinLayers(obj[["RNA"]])
#   obj <- NormalizeData(obj, verbose = FALSE)
#
# Normalization defaults are made explicit here for strict
# reproducibility.
# ============================================================

cat("\n")
cat("============================================\n")
cat("Preparing RNA expression data\n")
cat("============================================\n")

DefaultAssay(obj) <- "RNA"

cat("\nRNA layers before JoinLayers():\n")
print(
  Layers(
    obj[["RNA"]]
  )
)


# ------------------------------------------------------------
# Join Seurat v5 assay layers exactly as in historical workflow
# ------------------------------------------------------------

obj[["RNA"]] <- JoinLayers(
  obj[["RNA"]]
)

cat("\nRNA layers after JoinLayers():\n")
print(
  Layers(
    obj[["RNA"]]
  )
)


# ------------------------------------------------------------
# Confirm raw counts are available before re-normalization
# ------------------------------------------------------------

if (!"counts" %in% Layers(
  obj[["RNA"]]
)) {
  stop(
    paste0(
      "Joined RNA assay does not contain a counts layer. ",
      "Historical NormalizeData() cannot be reproduced."
    )
  )
}


# ------------------------------------------------------------
# Historical re-normalization
# ------------------------------------------------------------

obj <- NormalizeData(
  object = obj,
  assay = "RNA",
  normalization.method =
    NORMALIZATION_METHOD,
  scale.factor =
    NORMALIZATION_SCALE_FACTOR,
  verbose = FALSE
)

if (!"data" %in% Layers(
  obj[["RNA"]]
)) {
  stop(
    "NormalizeData() did not produce an RNA data layer."
  )
}

cat(
  "\nRNA normalization complete.\n"
)


# ============================================================
# 10. DEFINE CELL TYPES
# ============================================================

cell_types <- sort(
  unique(
    as.character(
      obj$cell_type
    )
  )
)

cell_types <- cell_types[
  !is.na(cell_types) &
    nzchar(cell_types)
]

cat("\nCell types:\n")
print(cell_types)


# ============================================================
# 11. CELL COUNTS BY OUTCOME
# ============================================================

cell_counts <- obj@meta.data %>%
  as.data.frame() %>%
  count(
    cell_type,
    de_outcome,
    name = "n_cells"
  ) %>%
  arrange(
    cell_type,
    de_outcome
  )

write.csv(
  cell_counts,
  file.path(
    table_dir,
    "03_pre_bcg_cell_counts_by_outcome.csv"
  ),
  row.names = FALSE
)

cat(
  "\nCells per cell type and outcome:\n"
)

print(
  as.data.frame(
    cell_counts
  ),
  row.names = FALSE
)


# ============================================================
# 12. CELL-LEVEL WILCOXON DIFFERENTIAL EXPRESSION
#
# Historical settings:
#
#   ident.1         = recurrence
#   ident.2         = remission
#   assay           = RNA
#   slot            = data
#   test.use        = wilcox
#   logfc.threshold = 0.25
#   min.pct         = 0.10
#
# A cell type is tested only when each outcome group contains
# at least 25 cells.
# ============================================================

cat("\n")
cat("============================================\n")
cat("Running pre-BCG differential expression\n")
cat("============================================\n")

de_results <- list()

de_gene_counts <- data.frame(
  cell_type =
    cell_types,

  n_DE_genes =
    NA_integer_,

  recurrence_cells =
    NA_integer_,

  remission_cells =
    NA_integer_,

  genes_tested =
    NA_integer_,

  stringsAsFactors = FALSE
)


for (ct in cell_types) {

  cat(
    "\n--------------------------------------------\n"
  )

  cat(
    "Processing: ",
    ct,
    "\n",
    sep = ""
  )

  cells_ct <- colnames(obj)[
    as.character(
      obj$cell_type
    ) == ct
  ]

  tab <- table(
    obj$de_outcome[
      cells_ct
    ]
  )


  # ----------------------------------------------------------
  # Record cell counts
  # ----------------------------------------------------------

  n_recurrence <- if (
    "recurrence" %in%
      names(tab)
  ) {
    as.integer(
      tab[
        "recurrence"
      ]
    )
  } else {
    0L
  }

  n_remission <- if (
    "remission" %in%
      names(tab)
  ) {
    as.integer(
      tab[
        "remission"
      ]
    )
  } else {
    0L
  }

  de_gene_counts$recurrence_cells[
    de_gene_counts$cell_type ==
      ct
  ] <- n_recurrence

  de_gene_counts$remission_cells[
    de_gene_counts$cell_type ==
      ct
  ] <- n_remission


  cat(
    "  recurrence cells: ",
    n_recurrence,
    "\n",
    sep = ""
  )

  cat(
    "  remission cells: ",
    n_remission,
    "\n",
    sep = ""
  )


  # ----------------------------------------------------------
  # Historical minimum-cell rule
  # ----------------------------------------------------------

  if (
    !all(
      c(
        "remission",
        "recurrence"
      ) %in%
        names(tab)
    ) ||
    min(tab) <
      MIN_CELLS_PER_GROUP
  ) {

    cat(
      "  Skipped: fewer than ",
      MIN_CELLS_PER_GROUP,
      " cells in at least one outcome group.\n",
      sep = ""
    )

    next
  }


  # ----------------------------------------------------------
  # Subset current cell type
  # ----------------------------------------------------------

  obj_ct <- subset(
    obj,
    cells =
      cells_ct
  )

  DefaultAssay(obj_ct) <- "RNA"

  Idents(obj_ct) <-
    obj_ct$de_outcome


  # ----------------------------------------------------------
  # Historical FindMarkers analysis
  # ----------------------------------------------------------

  de <- FindMarkers(
    object =
      obj_ct,

    ident.1 =
      "recurrence",

    ident.2 =
      "remission",

    assay =
      "RNA",

    slot =
      "data",

    test.use =
      "wilcox",

    logfc.threshold =
      LOGFC_THRESHOLD,

    min.pct =
      MIN_PCT,

    verbose =
      FALSE
  )


  # ----------------------------------------------------------
  # Add identifiers
  # ----------------------------------------------------------

  de$gene <-
    rownames(de)

  de$cell_type <-
    ct

  de$comparison <-
    "early_recurrence_vs_recurrence_free"

  de$direction <- ifelse(
    de$avg_log2FC > 0,
    "higher_in_early_recurrence",
    "higher_in_recurrence_free"
  )

  de$significant <-
    de$p_val_adj <
      PADJ_THRESHOLD


  # ----------------------------------------------------------
  # Save result
  # ----------------------------------------------------------

  de_results[[ct]] <-
    de

  safe_ct <- gsub(
    "[^A-Za-z0-9]+",
    "_",
    ct
  )

  write.csv(
    de,
    file.path(
      table_dir,
      paste0(
        "03_pre_bcg_DE_",
        safe_ct,
        ".csv"
      )
    ),
    row.names = FALSE
  )


  # ----------------------------------------------------------
  # Count significant genes exactly as historical code did
  # ----------------------------------------------------------

  n_de <- sum(
    de$p_val_adj <
      PADJ_THRESHOLD,
    na.rm = TRUE
  )

  de_gene_counts$n_DE_genes[
    de_gene_counts$cell_type ==
      ct
  ] <- n_de

  de_gene_counts$genes_tested[
    de_gene_counts$cell_type ==
      ct
  ] <- nrow(de)


  cat(
    "  Genes tested: ",
    nrow(de),
    "\n",
    sep = ""
  )

  cat(
    "  Significant DEGs: ",
    n_de,
    "\n",
    sep = ""
  )
}


# ============================================================
# 13. COMBINE AND EXPORT DE RESULTS
# ============================================================

if (length(
  de_results
) == 0) {
  stop(
    "No cell types passed the differential-expression criteria."
  )
}

all_de <- bind_rows(
  de_results
)

significant_de <- all_de %>%
  filter(
    p_val_adj <
      PADJ_THRESHOLD
  )


write.csv(
  all_de,
  file.path(
    table_dir,
    "03_pre_bcg_all_DE_results.csv"
  ),
  row.names = FALSE
)

write.csv(
  significant_de,
  file.path(
    table_dir,
    "03_pre_bcg_significant_DEGs.csv"
  ),
  row.names = FALSE
)

write.csv(
  de_gene_counts,
  file.path(
    table_dir,
    "03_pre_bcg_DEG_counts.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 14. SUMMARIZE DIRECTION OF DIFFERENTIAL EXPRESSION
# ============================================================

direction_counts <- significant_de %>%
  count(
    cell_type,
    direction,
    name = "n_DE_genes"
  ) %>%
  arrange(
    cell_type,
    direction
  )

write.csv(
  direction_counts,
  file.path(
    table_dir,
    "03_pre_bcg_DEG_counts_by_direction.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 15. MANUSCRIPT REPRODUCIBILITY CHECK
#
# These values are used only as an audit check. They do not
# alter the analysis.
# ============================================================

manuscript_targets <- c(
  "Classical Monocyte" =
    338,

  "Activated Monocyte" =
    220,

  "Non-classical Monocyte" =
    154,

  "T cell (CD4+)" =
    210,

  "T cell (CD8+)" =
    155,

  "NK Cell" =
    147
)

manuscript_total <- 1367


comparison <- data.frame(
  cell_type =
    names(
      manuscript_targets
    ),

  manuscript_DEGs =
    as.integer(
      manuscript_targets
    ),

  stringsAsFactors =
    FALSE
)

comparison$recomputed_DEGs <-
  de_gene_counts$n_DE_genes[
    match(
      comparison$cell_type,
      de_gene_counts$cell_type
    )
  ]

comparison$difference <-
  comparison$recomputed_DEGs -
  comparison$manuscript_DEGs


write.csv(
  comparison,
  file.path(
    table_dir,
    "03_pre_bcg_manuscript_count_check.csv"
  ),
  row.names = FALSE
)


total_recomputed <- sum(
  de_gene_counts$n_DE_genes,
  na.rm = TRUE
)

# ============================================================
# 16. RECREATE FIGURE 1B BAR PLOT
#
# Historical plot:
#
#   x = # DE genes (padj < 0.05)
#   y = cell type, ordered by DEG count
#   fill = DEG count
#   lightgrey -> blue gradient
# ============================================================

deg_df <- de_gene_counts %>%
  filter(
    !is.na(
      n_DE_genes
    )
  ) %>%
  mutate(
    cell_type =
      reorder(
        cell_type,
        n_DE_genes
      )
  )


fig1b <- ggplot(
  deg_df,
  aes(
    x =
      n_DE_genes,
    y =
      cell_type,
    fill =
      n_DE_genes
  )
) +
  geom_col() +
  scale_fill_gradientn(
    colours = c(
      "lightgrey",
      "blue"
    ),
    name =
      "# DE genes"
  ) +
  labs(
    x =
      "# DE genes (padj < 0.05)",
    y =
      NULL,
    title =
      "Differential Expression Burden by Cell Type"
  ) +
  theme_minimal(
    base_size = 12
  ) +
  theme(
    axis.text.y =
      element_text(
        size = 10,
        face = "bold"
      ),

    plot.title =
      element_text(
        size = 14,
        face = "bold"
      ),

    panel.grid.major.y =
      element_blank(),

    panel.grid.minor =
      element_blank()
  )


ggsave(
  filename = file.path(
    figure_dir,
    "fig1b_pre_bcg_DEG_burden.pdf"
  ),
  plot =
    fig1b,
  width =
    6,
  height =
    5,
  units =
    "in"
)

ggsave(
  filename = file.path(
    figure_dir,
    "fig1b_pre_bcg_DEG_burden.png"
  ),
  plot =
    fig1b,
  width =
    6,
  height =
    5,
  units =
    "in",
  dpi =
    600,
  bg =
    "transparent"
)


# ============================================================
# 17. SAVE ANALYSIS PARAMETERS
# ============================================================

analysis_parameters <- data.frame(
  parameter = c(
    "input_object",
    "timepoint",
    "assay",
    "expression_slot",
    "normalization_method",
    "normalization_scale_factor",
    "test",
    "ident_1",
    "ident_2",
    "logfc_threshold",
    "min_pct",
    "minimum_cells_per_group",
    "significance_column",
    "significance_threshold"
  ),

  value = c(
    basename(
      input_rds
    ),
    "pre-BCG",
    "RNA",
    "data",
    NORMALIZATION_METHOD,
    NORMALIZATION_SCALE_FACTOR,
    "wilcox",
    "early_recurrence",
    "recurrence_free",
    LOGFC_THRESHOLD,
    MIN_PCT,
    MIN_CELLS_PER_GROUP,
    "p_val_adj",
    PADJ_THRESHOLD
  ),

  stringsAsFactors =
    FALSE
)


write.csv(
  analysis_parameters,
  file.path(
    table_dir,
    "03_pre_bcg_DE_parameters.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 18. SESSION INFORMATION
# ============================================================

capture.output(
  sessionInfo(),
  file = file.path(
    table_dir,
    "03_pre_bcg_DE_session_info.txt"
  )
)


# ============================================================
# 19. FINAL REPORT
# ============================================================

cat("\n")
cat("============================================\n")
cat("PRE-BCG DIFFERENTIAL EXPRESSION COMPLETE\n")
cat("============================================\n\n")

cat(
  "Cell types assessed: ",
  sum(
    !is.na(
      de_gene_counts$n_DE_genes
    )
  ),
  "\n",
  sep = ""
)

cat(
  "Total significant DEG associations: ",
  total_recomputed,
  "\n",
  sep = ""
)

cat(
  "\nPositive avg_log2FC = higher in early-recurrence cells.\n"
)

cat(
  "Negative avg_log2FC = higher in recurrence-free cells.\n"
)

cat(
  "\nTables written to:\n  ",
  table_dir,
  "\n",
  sep = ""
)

cat(
  "\nFigure written to:\n  ",
  file.path(
    figure_dir,
    "fig1b_pre_bcg_DEG_burden.pdf"
  ),
  "\n",
  sep = ""
)

cat(
  "\nAnalysis complete.\n"
)
