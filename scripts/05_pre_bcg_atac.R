# ============================================================
# 05_pre_bcg_atac.R
#
# Pre-BCG differential accessibility + motif analysis
#
# Intended outputs:
#   - Figure 1d: number of differentially accessible peaks
#                 by cell type
#   - Figure 1e: differential transcription-factor motif
#                 accessibility, highlighting AP-1 motifs
#   - Per-cell-type DA peak tables
#   - Per-cell-type chromVAR differential motif tables
#   - Optional motif enrichment among recurrence-up DA peaks
#
# INPUT:
#   Historical finalized pre-BCG Seurat/Signac object:
#       G:/merged_obj.rds
#
# IMPORTANT PROVENANCE NOTES
# ------------------------------------------------------------
# The peak-level DA settings below were recovered directly
# from the historical analysis:
#
#   DefaultAssay(object) <- "ATAC"
#   Idents(object) <- "Outcome"
#   FindMarkers(
#       ident.1 = "Recurrence",
#       ident.2 = "Remission",
#       test.use = "LR",
#       only.pos = FALSE,
#       min.pct = 0.05,
#       logfc.threshold = 0
#   )
#
# The historical object does not retain a Motif/chromVAR assay,
# so motif annotations and chromVAR deviations are regenerated
# here using JASPAR2020 + hg38.
#
# A historical motif-comparison block used:
#   test.use = "wilcox"
#   min.pct = 0.05
#   logfc.threshold = 0
#
# We apply those settings here to the baseline
# Recurrence-vs-Remission comparison within each cell type.
#
# Positive peak avg_log2FC / positive motif avg_log2FC:
#   greater accessibility in Recurrence
#   (= early recurrence / historical non-responder group)
#
# ============================================================


# ============================================================
# 1. PACKAGES
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(ggrepel)
  library(Matrix)
  library(GenomicRanges)
  library(GenomeInfoDb)
  library(JASPAR2020)
  library(TFBSTools)
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(BiocParallel)
})


# ============================================================
# 2. PATHS
# ============================================================

input_rds <- "G:/merged_obj.rds"

outdir <- file.path(
  "results",
  "05_pre_bcg_atac"
)

da_dir <- file.path(
  outdir,
  "DA_peaks"
)

motif_dir <- file.path(
  outdir,
  "motifs"
)

enrichment_dir <- file.path(
  outdir,
  "motif_enrichment"
)

figure_dir <- file.path(
  outdir,
  "figures"
)

dir.create(
  outdir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  da_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  motif_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  enrichment_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  figure_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ============================================================
# 3. ANALYSIS PARAMETERS
# ============================================================

DA_TEST <- "LR"
DA_MIN_PCT <- 0.05
DA_LOGFC_THRESHOLD <- 0
DA_ONLY_POS <- FALSE
DA_PADJ_THRESHOLD <- 0.05

MOTIF_TEST <- "wilcox"
MOTIF_MIN_PCT <- 0.05
MOTIF_LOGFC_THRESHOLD <- 0
MOTIF_PADJ_THRESHOLD <- 0.05

IDENT_1 <- "Recurrence"
IDENT_2 <- "Remission"

# Historical AP-1 motif family list recovered from the
# original motif work.
AP1_NAMES <- c(
  "FOS",
  "JUN",
  "JUNB",
  "JUND",
  "FOSL1",
  "FOSL2",
  "BATF"
)

# Published Figure 1d values.
# AUDIT TARGETS ONLY: never used to tune the analysis.
MANUSCRIPT_DA_TARGETS <- c(
  "Classical Monocyte" = 5778,
  "Activated Monocyte" = 1470,
  "Non-classical Monocyte" = 416,
  "T cell (CD8+)" = 1033,
  "T cell (CD4+)" = 1090,
  "NK Cell" = 516
)


# ============================================================
# 4. HELPERS
# ============================================================

safe_name <- function(x) {
  x <- gsub("[^A-Za-z0-9_-]+", "_", x)
  x <- gsub("^_+|_+$", "", x)
  x
}


save_csv <- function(x, path, row_names = FALSE) {
  write.csv(
    x,
    file = path,
    row.names = row_names
  )
}


check_required_metadata <- function(object) {

  required <- c(
    "cell_type",
    "Outcome"
  )

  missing <- setdiff(
    required,
    colnames(object@meta.data)
  )

  if (length(missing) > 0) {
    stop(
      paste0(
        "Object is missing required metadata column(s): ",
        paste(missing, collapse = ", ")
      )
    )
  }

  invisible(TRUE)
}


get_motif_name_vector <- function(object) {

  motif_names <- GetMotifData(
    object = object,
    assay = "ATAC",
    slot = "motif.names"
  )

  # Depending on Signac version this may be either a named
  # character vector or a one-column/data-frame-like object.

  if (is.data.frame(motif_names)) {

    if ("motif.name" %in% colnames(motif_names)) {

      out <- motif_names$motif.name
      names(out) <- rownames(motif_names)
      return(out)

    }

    if (ncol(motif_names) == 1) {

      out <- motif_names[[1]]
      names(out) <- rownames(motif_names)
      return(out)

    }

  }

  if (is.list(motif_names) && !is.atomic(motif_names)) {

    out <- unlist(motif_names)
    return(out)

  }

  motif_names
}


collapse_motif_names <- function(x) {

  if (length(x) == 0) {
    return(NA_character_)
  }

  paste(
    unique(as.character(x)),
    collapse = ";"
  )
}


# ============================================================
# 5. LOAD HISTORICAL PRE-BCG OBJECT
# ============================================================

cat("\n")
cat("============================================================\n")
cat("Loading historical pre-BCG Seurat/Signac object\n")
cat("============================================================\n")

if (!file.exists(input_rds)) {
  stop(
    paste0(
      "Input object not found: ",
      input_rds
    )
  )
}

pre_bcg <- readRDS(
  input_rds
)

check_required_metadata(
  pre_bcg
)

if (!"ATAC" %in% Assays(pre_bcg)) {
  stop(
    "ATAC assay not found in merged_obj.rds."
  )
}

cat(
  "Cells: ",
  format(
    ncol(pre_bcg),
    big.mark = ","
  ),
  "\n",
  sep = ""
)

cat(
  "ATAC features: ",
  format(
    nrow(pre_bcg[["ATAC"]]),
    big.mark = ","
  ),
  "\n",
  sep = ""
)

cat("\nOutcome counts:\n")
print(
  table(
    pre_bcg$Outcome,
    useNA = "ifany"
  )
)

cat("\nCell-type counts:\n")
print(
  sort(
    table(pre_bcg$cell_type),
    decreasing = TRUE
  )
)


# ============================================================
# 6. VALIDATE OUTCOME LABELS
# ============================================================

outcome_levels <- unique(
  as.character(
    pre_bcg$Outcome
  )
)

if (!all(
  c(
    IDENT_1,
    IDENT_2
  ) %in% outcome_levels
)) {

  stop(
    paste0(
      "Expected Outcome labels '",
      IDENT_1,
      "' and '",
      IDENT_2,
      "'. Observed: ",
      paste(
        sort(outcome_levels),
        collapse = ", "
      )
    )
  )
}


# ============================================================
# 7. PEAK-LEVEL DIFFERENTIAL ACCESSIBILITY
#
# Historical exact settings:
#
# FindMarkers(
#   object = object_ct,
#   ident.1 = "Recurrence",
#   ident.2 = "Remission",
#   test.use = "LR",
#   only.pos = FALSE,
#   min.pct = 0.05,
#   logfc.threshold = 0
# )
#
# No latent variable was present in the recovered historical
# baseline block, so none is added here.
# ============================================================

cat("\n")
cat("============================================================\n")
cat("Peak-level differential accessibility\n")
cat("============================================================\n")

DefaultAssay(
  pre_bcg
) <- "ATAC"

cell_types <- sort(
  unique(
    as.character(
      pre_bcg$cell_type
    )
  )
)

da_tables <- list()
da_summary <- list()


for (ct in cell_types) {

  cat("\n------------------------------------------------------------\n")
  cat("Cell type: ", ct, "\n", sep = "")

  cells_ct <- WhichCells(
    pre_bcg,
    expression = cell_type == ct
  )

  if (length(cells_ct) == 0) {
    next
  }

  group_n <- table(
    pre_bcg$Outcome[
      match(
        cells_ct,
        colnames(pre_bcg)
      )
    ]
  )

  n_rec <- if (
    IDENT_1 %in% names(group_n)
  ) {
    as.integer(
      group_n[[IDENT_1]]
    )
  } else {
    0L
  }

  n_rem <- if (
    IDENT_2 %in% names(group_n)
  ) {
    as.integer(
      group_n[[IDENT_2]]
    )
  } else {
    0L
  }

  cat(
    "  Recurrence cells: ",
    format(n_rec, big.mark = ","),
    "\n",
    sep = ""
  )

  cat(
    "  Remission cells: ",
    format(n_rem, big.mark = ","),
    "\n",
    sep = ""
  )

  if (
    n_rec == 0L ||
    n_rem == 0L
  ) {

    warning(
      paste0(
        "Skipping ",
        ct,
        ": one outcome group is absent."
      )
    )

    next
  }

  object_ct <- subset(
    pre_bcg,
    cells = cells_ct
  )

  DefaultAssay(
    object_ct
  ) <- "ATAC"

  Idents(
    object_ct
  ) <- "Outcome"

  da <- FindMarkers(
    object = object_ct,
    ident.1 = IDENT_1,
    ident.2 = IDENT_2,
    test.use = DA_TEST,
    only.pos = DA_ONLY_POS,
    min.pct = DA_MIN_PCT,
    logfc.threshold = DA_LOGFC_THRESHOLD
  )

  da$peak <- rownames(
    da
  )

  da$cell_type <- ct

  da$comparison <- paste0(
    IDENT_1,
    "_vs_",
    IDENT_2
  )

  da$direction <- ifelse(
    da$avg_log2FC > 0,
    "higher_in_Recurrence",
    ifelse(
      da$avg_log2FC < 0,
      "higher_in_Remission",
      "no_change"
    )
  )

  da$significant <- (
    !is.na(
      da$p_val_adj
    ) &
      da$p_val_adj < DA_PADJ_THRESHOLD
  )

  da <- da %>%
    arrange(
      p_val_adj,
      desc(
        abs(avg_log2FC)
      )
    )

  da_tables[[ct]] <- da

  n_sig <- sum(
    da$significant,
    na.rm = TRUE
  )

  n_sig_rec <- sum(
    da$significant &
      da$avg_log2FC > 0,
    na.rm = TRUE
  )

  n_sig_rem <- sum(
    da$significant &
      da$avg_log2FC < 0,
    na.rm = TRUE
  )

  da_summary[[ct]] <- data.frame(
    cell_type = ct,
    recurrence_cells = n_rec,
    remission_cells = n_rem,
    peaks_tested = nrow(da),
    significant_DA_peaks = n_sig,
    significant_higher_in_Recurrence = n_sig_rec,
    significant_higher_in_Remission = n_sig_rem,
    stringsAsFactors = FALSE
  )

  cat(
    "  Significant DA peaks: ",
    format(
      n_sig,
      big.mark = ","
    ),
    "\n",
    sep = ""
  )

  save_csv(
    da,
    file.path(
      da_dir,
      paste0(
        safe_name(ct),
        "_Recurrence_vs_Remission_DA_peaks.csv"
      )
    ),
    row_names = FALSE
  )

  rm(
    object_ct
  )

  gc()
}


if (length(da_tables) == 0) {
  stop(
    "No cell type produced a differential-accessibility result."
  )
}


all_da <- bind_rows(
  da_tables
)

da_counts <- bind_rows(
  da_summary
) %>%
  arrange(
    desc(
      significant_DA_peaks
    )
  )


save_csv(
  all_da,
  file.path(
    outdir,
    "05_pre_bcg_all_DA_peaks.csv"
  ),
  row_names = FALSE
)

save_csv(
  da_counts,
  file.path(
    outdir,
    "05_pre_bcg_DA_peak_counts_by_cell_type.csv"
  ),
  row_names = FALSE
)


# ============================================================
# 8. MANUSCRIPT FIGURE 1D AUDIT
#
# These values are audit targets only.
# ============================================================

audit_da <- data.frame(
  cell_type = names(
    MANUSCRIPT_DA_TARGETS
  ),
  reported_n_DAR = as.integer(
    MANUSCRIPT_DA_TARGETS
  ),
  stringsAsFactors = FALSE
)

audit_da <- audit_da %>%
  left_join(
    da_counts %>%
      select(
        cell_type,
        recomputed_n_DAR = significant_DA_peaks
      ),
    by = "cell_type"
  ) %>%
  mutate(
    difference = recomputed_n_DAR - reported_n_DAR
  )


save_csv(
  audit_da,
  file.path(
    outdir,
    "05_pre_bcg_manuscript_DA_count_check.csv"
  ),
  row_names = FALSE
)


cat("\n")
cat("============================================================\n")
cat("Figure 1d manuscript audit\n")
cat("============================================================\n")

print(
  audit_da,
  row.names = FALSE
)


# ============================================================
# 9. FIGURE 1D BAR PLOT
# ============================================================

plot_da <- da_counts %>%
  mutate(
    cell_type = factor(
      cell_type,
      levels = rev(
        cell_type[
          order(
            significant_DA_peaks
          )
        ]
      )
    )
  )


p_da <- ggplot(
  plot_da,
  aes(
    x = cell_type,
    y = significant_DA_peaks
  )
) +
  geom_col(
    width = 0.75
  ) +
  coord_flip() +
  theme_classic(
    base_size = 12
  ) +
  labs(
    x = NULL,
    y = "# DA peaks (adjusted p < 0.05)",
    title = "Differential accessibility burden by cell type"
  )


ggsave(
  filename = file.path(
    figure_dir,
    "fig1d_pre_bcg_DA_peak_counts.pdf"
  ),
  plot = p_da,
  width = 7,
  height = 6,
  units = "in"
)

ggsave(
  filename = file.path(
    figure_dir,
    "fig1d_pre_bcg_DA_peak_counts.png"
  ),
  plot = p_da,
  width = 7,
  height = 6,
  units = "in",
  dpi = 600
)


# ============================================================
# 10. PREPARE A COPY FOR MOTIF ANALYSIS
#
# The historical merged_obj.rds does not contain a retained
# Motif object. Rebuild motif annotations from JASPAR2020.
#
# We keep the DA analysis above untouched.
#
# For motif sequence matching only, restrict the motif-analysis
# copy to standard chromosomes available in the hg38 BSgenome.
# ============================================================

cat("\n")
cat("============================================================\n")
cat("Preparing JASPAR2020 motif annotations\n")
cat("============================================================\n")

motif_obj <- pre_bcg

DefaultAssay(
  motif_obj
) <- "ATAC"

atac_ranges <- granges(
  motif_obj[["ATAC"]]
)

feature_names <- rownames(
  motif_obj[["ATAC"]]
)

names(
  atac_ranges
) <- feature_names

standard_chromosomes <- c(
  paste0(
    "chr",
    1:22
  ),
  "chrX",
  "chrY"
)

keep_features <- (
  !is.na(
    as.character(
      seqnames(
        atac_ranges
      )
    )
  ) &
    as.character(
      seqnames(
        atac_ranges
      )
    ) %in% standard_chromosomes
)

cat(
  "ATAC peaks retained for motif matching: ",
  sum(keep_features),
  " / ",
  length(keep_features),
  "\n",
  sep = ""
)

motif_obj <- subset(
  motif_obj,
  features = feature_names[
    keep_features
  ]
)

DefaultAssay(
  motif_obj
) <- "ATAC"


# ============================================================
# 11. ADD JASPAR2020 MOTIFS
#
# Exact recovered JASPAR options:
#   collection = "CORE"
#   tax_group  = "vertebrates"
# ============================================================

register(
  SerialParam()
)

pfm <- getMatrixSet(
  JASPAR2020,
  opts = list(
    collection = "CORE",
    tax_group = "vertebrates"
  )
)

motif_obj <- AddMotifs(
  object = motif_obj,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  pfm = pfm
)


motif_name_vector <- get_motif_name_vector(
  motif_obj
)

cat(
  "Motifs annotated: ",
  length(
    motif_name_vector
  ),
  "\n",
  sep = ""
)


# ============================================================
# 12. MOTIF ENRICHMENT IN RECURRENCE-UP DA PEAKS
#
# This tests whether motif instances are overrepresented among
# significant DA peaks with greater accessibility in Recurrence.
#
# This is complementary to the chromVAR analysis below.
# ============================================================

cat("\n")
cat("============================================================\n")
cat("Motif enrichment in Recurrence-up DA peaks\n")
cat("============================================================\n")

motif_enrichment_tables <- list()


for (ct in names(da_tables)) {

  cat(
    "\nCell type: ",
    ct,
    "\n",
    sep = ""
  )

  da <- da_tables[[ct]]

  recurrence_up_peaks <- da %>%
    filter(
      p_val_adj < DA_PADJ_THRESHOLD,
      avg_log2FC > 0
    ) %>%
    pull(
      peak
    )

  # Only peaks present in the motif-analysis copy can be tested.
  recurrence_up_peaks <- intersect(
    recurrence_up_peaks,
    rownames(
      motif_obj[["ATAC"]]
    )
  )

  if (
    length(
      recurrence_up_peaks
    ) < 10
  ) {

    cat(
      "  Skipping enrichment: fewer than 10 eligible Recurrence-up peaks.\n"
    )

    next
  }

  cells_ct <- WhichCells(
    motif_obj,
    expression = cell_type == ct
  )

  if (
    length(
      cells_ct
    ) == 0
  ) {
    next
  }

  motif_ct <- subset(
    motif_obj,
    cells = cells_ct
  )

  DefaultAssay(
    motif_ct
  ) <- "ATAC"

  background_peaks <- intersect(
    da$peak,
    rownames(
      motif_ct[["ATAC"]]
    )
  )

  if (
    length(
      background_peaks
    ) == 0
  ) {
    next
  }

  motif_enrichment <- FindMotifs(
    object = motif_ct,
    features = recurrence_up_peaks,
    background = background_peaks
  )

  motif_enrichment$motif_id <- rownames(
    motif_enrichment
  )

  motif_enrichment$cell_type <- ct

  motif_enrichment$comparison <- "Recurrence_up_DA_peaks"

  motif_enrichment$AP1 <- grepl(
    paste(
      AP1_NAMES,
      collapse = "|"
    ),
    motif_enrichment$motif.name,
    ignore.case = TRUE
  )

  motif_enrichment <- motif_enrichment %>%
    arrange(
      p.adjust
    )

  motif_enrichment_tables[[ct]] <- motif_enrichment

  save_csv(
    motif_enrichment,
    file.path(
      enrichment_dir,
      paste0(
        safe_name(ct),
        "_Recurrence_up_DA_peak_motif_enrichment.csv"
      )
    ),
    row_names = FALSE
  )

  cat(
    "  Recurrence-up DA peaks tested: ",
    length(
      recurrence_up_peaks
    ),
    "\n",
    sep = ""
  )

  rm(
    motif_ct
  )

  gc()
}


if (
  length(
    motif_enrichment_tables
  ) > 0
) {

  all_motif_enrichment <- bind_rows(
    motif_enrichment_tables
  )

  save_csv(
    all_motif_enrichment,
    file.path(
      outdir,
      "05_pre_bcg_all_motif_enrichment.csv"
    ),
    row_names = FALSE
  )

}


# ============================================================
# 13. RUN chromVAR
#
# Regenerates motif deviation scores because the historical
# merged_obj.rds does not retain the chromVAR assay.
# ============================================================

cat("\n")
cat("============================================================\n")
cat("Running chromVAR")
cat("============================================================\n")

motif_obj <- RunChromVAR(
  object = motif_obj,
  genome = BSgenome.Hsapiens.UCSC.hg38
)

if (
  !"chromvar" %in% Assays(
    motif_obj
  )
) {
  stop(
    "RunChromVAR completed without creating a chromvar assay."
  )
}


# ============================================================
# 14. DIFFERENTIAL MOTIF ACCESSIBILITY BY CELL TYPE
#
# Recovered historical motif-comparison settings:
#
#   FindMarkers(
#       test.use = "wilcox",
#       min.pct = 0.05,
#       logfc.threshold = 0
#   )
#
# Adapted here to Recurrence vs Remission within each
# pre-BCG cell type.
# ============================================================

cat("\n")
cat("============================================================\n")
cat("Differential chromVAR motif accessibility\n")
cat("============================================================\n")

DefaultAssay(
  motif_obj
) <- "chromvar"

motif_tables <- list()


for (ct in cell_types) {

  cat(
    "\nCell type: ",
    ct,
    "\n",
    sep = ""
  )

  cells_ct <- WhichCells(
    motif_obj,
    expression = cell_type == ct
  )

  if (
    length(
      cells_ct
    ) == 0
  ) {
    next
  }

  motif_ct <- subset(
    motif_obj,
    cells = cells_ct
  )

  DefaultAssay(
    motif_ct
  ) <- "chromvar"

  Idents(
    motif_ct
  ) <- "Outcome"

  groups_present <- unique(
    as.character(
      Idents(
        motif_ct
      )
    )
  )

  if (!all(
    c(
      IDENT_1,
      IDENT_2
    ) %in% groups_present
  )) {

    warning(
      paste0(
        "Skipping motif comparison for ",
        ct,
        ": one outcome group is absent."
      )
    )

    next
  }

  motif_da <- FindMarkers(
    object = motif_ct,
    ident.1 = IDENT_1,
    ident.2 = IDENT_2,
    test.use = MOTIF_TEST,
    min.pct = MOTIF_MIN_PCT,
    logfc.threshold = MOTIF_LOGFC_THRESHOLD
  )

  motif_da$motif_id <- rownames(
    motif_da
  )

  motif_da$tf <- vapply(
    motif_da$motif_id,
    function(id) {

      if (
        id %in% names(
          motif_name_vector
        )
      ) {

        collapse_motif_names(
          motif_name_vector[
            id
          ]
        )

      } else {

        NA_character_

      }

    },
    FUN.VALUE = character(1)
  )

  motif_da$cell_type <- ct

  motif_da$comparison <- paste0(
    IDENT_1,
    "_vs_",
    IDENT_2
  )

  motif_da$AP1 <- grepl(
    paste(
      AP1_NAMES,
      collapse = "|"
    ),
    motif_da$tf,
    ignore.case = TRUE
  )

  motif_da$direction <- ifelse(
    motif_da$avg_log2FC > 0,
    "higher_in_Recurrence",
    ifelse(
      motif_da$avg_log2FC < 0,
      "higher_in_Remission",
      "no_change"
    )
  )

  motif_da$significant <- (
    !is.na(
      motif_da$p_val_adj
    ) &
      motif_da$p_val_adj < MOTIF_PADJ_THRESHOLD
  )

  motif_da <- motif_da %>%
    arrange(
      p_val_adj,
      desc(
        abs(avg_log2FC)
      )
    )

  motif_tables[[ct]] <- motif_da

  save_csv(
    motif_da,
    file.path(
      motif_dir,
      paste0(
        safe_name(ct),
        "_Recurrence_vs_Remission_chromVAR_motifs.csv"
      )
    ),
    row_names = FALSE
  )

  rm(
    motif_ct
  )

  gc()
}


if (
  length(
    motif_tables
  ) == 0
) {

  stop(
    "No cell type produced a chromVAR differential motif result."
  )
}


all_motif_da <- bind_rows(
  motif_tables
)

save_csv(
  all_motif_da,
  file.path(
    outdir,
    "05_pre_bcg_all_chromVAR_motif_results.csv"
  ),
  row_names = FALSE
)


# ============================================================
# 15. AP-1 SUMMARY
# ============================================================

ap1_summary <- all_motif_da %>%
  filter(
    AP1
  ) %>%
  arrange(
    cell_type,
    p_val_adj
  )


save_csv(
  ap1_summary,
  file.path(
    outdir,
    "05_pre_bcg_AP1_chromVAR_results.csv"
  ),
  row_names = FALSE
)


ap1_sig_summary <- all_motif_da %>%
  filter(
    AP1,
    p_val_adj < MOTIF_PADJ_THRESHOLD
  ) %>%
  count(
    cell_type,
    direction,
    name = "n_significant_AP1_motifs"
  )


save_csv(
  ap1_sig_summary,
  file.path(
    outdir,
    "05_pre_bcg_AP1_significant_summary.csv"
  ),
  row_names = FALSE
)


# ============================================================
# 16. FIGURE 1E-LIKE MOTIF VOLCANO
#
# To give one point per motif for the overview plot, calculate
# the mean recurrence-vs-remission chromVAR effect across cell
# types and combine evidence using the minimum adjusted p-value.
#
# IMPORTANT:
# Per-cell-type statistical results are preserved in the output
# tables above. This summary is a visualization only.
# ============================================================

motif_overview <- all_motif_da %>%
  filter(
    !is.na(tf)
  ) %>%
  group_by(
    motif_id,
    tf
  ) %>%
  summarise(
    mean_accessibility_difference = mean(
      avg_log2FC,
      na.rm = TRUE
    ),
    min_adjusted_p = min(
      p_val_adj,
      na.rm = TRUE
    ),
    n_cell_types = n(),
    n_positive_cell_types = sum(
      avg_log2FC > 0,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  mutate(
    AP1 = grepl(
      paste(
        AP1_NAMES,
        collapse = "|"
      ),
      tf,
      ignore.case = TRUE
    ),
    neglog10_FDR = -log10(
      pmax(
        min_adjusted_p,
        .Machine$double.xmin
      )
    )
  )


save_csv(
  motif_overview,
  file.path(
    outdir,
    "05_pre_bcg_motif_overview_for_plot.csv"
  ),
  row_names = FALSE
)


p_motif <- ggplot(
  motif_overview,
  aes(
    x = mean_accessibility_difference,
    y = neglog10_FDR
  )
) +
  geom_point(
    aes(
      shape = AP1
    ),
    alpha = 0.65,
    size = 2
  ) +
  geom_vline(
    xintercept = 0,
    linetype = 2
  ) +
  geom_hline(
    yintercept = -log10(
      MOTIF_PADJ_THRESHOLD
    ),
    linetype = 2
  ) +
  geom_text_repel(
    data = motif_overview %>%
      filter(
        AP1
      ),
    aes(
      label = tf
    ),
    size = 3,
    max.overlaps = Inf
  ) +
  theme_classic(
    base_size = 12
  ) +
  labs(
    x = paste0(
      "Mean motif accessibility difference\n",
      "(Recurrence - Remission)"
    ),
    y = "-log10 minimum adjusted p-value",
    title = "Pre-BCG differential motif accessibility",
    shape = "AP-1 family"
  )


ggsave(
  filename = file.path(
    figure_dir,
    "fig1e_pre_bcg_motif_accessibility_overview.pdf"
  ),
  plot = p_motif,
  width = 7,
  height = 6,
  units = "in"
)

ggsave(
  filename = file.path(
    figure_dir,
    "fig1e_pre_bcg_motif_accessibility_overview.png"
  ),
  plot = p_motif,
  width = 7,
  height = 6,
  units = "in",
  dpi = 600
)


# ============================================================
# 17. SAVE PARAMETERS
# ============================================================

parameter_table <- data.frame(
  parameter = c(
    "input_rds",
    "ident_1",
    "ident_2",
    "DA_test",
    "DA_only_pos",
    "DA_min_pct",
    "DA_logfc_threshold",
    "DA_padj_threshold",
    "motif_database",
    "motif_collection",
    "motif_tax_group",
    "motif_genome",
    "chromVAR",
    "motif_test",
    "motif_min_pct",
    "motif_logfc_threshold",
    "motif_padj_threshold"
  ),
  value = c(
    input_rds,
    IDENT_1,
    IDENT_2,
    DA_TEST,
    DA_ONLY_POS,
    DA_MIN_PCT,
    DA_LOGFC_THRESHOLD,
    DA_PADJ_THRESHOLD,
    "JASPAR2020",
    "CORE",
    "vertebrates",
    "BSgenome.Hsapiens.UCSC.hg38",
    "Signac::RunChromVAR",
    MOTIF_TEST,
    MOTIF_MIN_PCT,
    MOTIF_LOGFC_THRESHOLD,
    MOTIF_PADJ_THRESHOLD
  ),
  stringsAsFactors = FALSE
)


save_csv(
  parameter_table,
  file.path(
    outdir,
    "05_pre_bcg_atac_parameters.csv"
  ),
  row_names = FALSE
)


# ============================================================
# 18. SESSION INFO
# ============================================================

sink(
  file.path(
    outdir,
    "05_pre_bcg_atac_sessionInfo.txt"
  )
)

print(
  sessionInfo()
)

sink()


# ============================================================
# 19. FINAL REPORT
# ============================================================

cat("\n")
cat("============================================================\n")
cat("PRE-BCG ATAC ANALYSIS COMPLETE\n")
cat("============================================================\n")

cat(
  "\nCell types with DA results: ",
  nrow(
    da_counts
  ),
  "\n",
  sep = ""
)

cat(
  "Total significant DA peak calls across cell types: ",
  format(
    sum(
      da_counts$significant_DA_peaks,
      na.rm = TRUE
    ),
    big.mark = ","
  ),
  "\n",
  sep = ""
)

cat(
  "Cell types with chromVAR motif comparisons: ",
  length(
    motif_tables
  ),
  "\n",
  sep = ""
)

cat(
  "Significant AP-1 motif/cell-type associations: ",
  sum(
    ap1_summary$p_val_adj < MOTIF_PADJ_THRESHOLD,
    na.rm = TRUE
  ),
  "\n",
  sep = ""
)

cat(
  "\nPositive DA peak avg_log2FC = greater accessibility in Recurrence.\n"
)

cat(
  "Positive chromVAR avg_log2FC = greater motif accessibility in Recurrence.\n"
)

cat(
  "\nOutputs written to:\n  ",
  normalizePath(
    outdir,
    winslash = "/",
    mustWork = FALSE
  ),
  "\n",
  sep = ""
)
