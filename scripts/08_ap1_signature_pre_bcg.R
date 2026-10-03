# ============================================================
# 08_ap1_signature_pre_bcg.R
#
# Reviewer-facing reproducibility demonstration for the
# pre-BCG AP-1 myeloid signature analysis.
#
# INPUT:
#   G:/merged_obj.rds
#
# ANALYSIS:
#   - Restrict to pre-BCG Classical Monocytes +
#     Activated Monocytes
#   - AP-1 signature:
#       JUN, JUNB, JUND, FOS, FOSB, FOSL2
#   - Seurat AddModuleScore at the single-cell level
#   - Aggregate mean module score per sample
#   - Two-sided Wilcoxon rank-sum test:
#       Recurrence vs Remission
#   - ROC analysis with higher score = recurrence
#   - AUC + 95% CI
#   - Rank-biserial effect size
#   - Per-gene expression/fraction summaries
#
# Historical labels:
#   Remission  = responder
#   Recurrence = non-responder / early recurrence
#
# IMPORTANT:
# The six-gene signature below follows the manuscript Methods:
#   JUN, JUNB, JUND, FOS, FOSB, FOSL2
#
# FOSL1 appeared in some earlier exploratory code / figure
# captions, but is intentionally NOT included here because the
# final Methods define the targeted signature as these six genes.
#
# Run from repository root:
#
# & "C:\Users\andrew\AppData\Local\Programs\R\R-4.4.1\bin\Rscript.exe" `
#     scripts\08_ap1_signature_pre_bcg.R
#
# ============================================================


# ============================================================
# 1. PACKAGES
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(ggpubr)
  library(pROC)
})


# ============================================================
# 2. PATHS
# ============================================================

input_rds <- "G:/merged_obj.rds"

outdir <- file.path(
  "results",
  "08_ap1_signature_pre_bcg"
)

table_dir <- file.path(
  outdir,
  "tables"
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
# 3. PARAMETERS
# ============================================================

AP1_GENES <- c(
  "JUN",
  "JUNB",
  "JUND",
  "FOS",
  "FOSB",
  "FOSL2"
)

MYELOID_CELL_TYPES <- c(
  "Classical Monocyte",
  "Activated Monocyte"
)

RESPONDER_LABEL <- "Remission"
NONRESPONDER_LABEL <- "Recurrence"

MODULE_NAME <- "AP1_signature"

# AddModuleScore samples matched control genes.
# Explicit seed makes this deterministic.
MODULE_SEED <- 1
MODULE_NBIN <- 24
MODULE_CTRL <- 100


# ============================================================
# 4. HELPERS
# ============================================================

save_csv <- function(
  x,
  path
) {

  write.csv(
    x,
    file = path,
    row.names = FALSE
  )

}


# Rank-biserial correlation for two independent groups.
#
# Here:
#   group 1 = Recurrence
#   group 2 = Remission
#
# Positive r_rb means scores tend to be higher in Recurrence.
rank_biserial_from_u <- function(
  x,
  y
) {

  n1 <- length(
    x
  )

  n2 <- length(
    y
  )

  ranks <- rank(
    c(
      x,
      y
    ),
    ties.method = "average"
  )

  rank_sum_x <- sum(
    ranks[
      seq_len(
        n1
      )
    ]
  )

  U1 <- rank_sum_x -
    n1 *
      (
        n1 +
          1
      ) /
      2

  r_rb <- (
    2 *
      U1 /
      (
        n1 *
          n2
      )
  ) -
    1

  list(
    U = U1,
    rank_biserial = r_rb
  )

}


# ============================================================
# 5. LOAD PRE-BCG OBJECT
# ============================================================

cat("\n")
cat("============================================================\n")
cat("Loading pre-BCG object\n")
cat("============================================================\n")

if (
  !file.exists(
    input_rds
  )
) {

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

required_metadata <- c(
  "cell_type",
  "sample_id",
  "Outcome"
)

missing_metadata <- setdiff(
  required_metadata,
  colnames(
    pre_bcg@meta.data
  )
)

if (
  length(
    missing_metadata
  ) >
    0
) {

  stop(
    paste0(
      "Missing required metadata column(s): ",
      paste(
        missing_metadata,
        collapse = ", "
      )
    )
  )

}

if (
  !"RNA" %in%
    Assays(
      pre_bcg
    )
) {

  stop(
    "RNA assay not found."
  )

}

cat(
  "Total pre-BCG cells: ",
  format(
    ncol(
      pre_bcg
    ),
    big.mark = ","
  ),
  "\n",
  sep = ""
)


# ============================================================
# 6. SUBSET MYELOID POPULATIONS
#
# Manuscript Methods:
#   analysis restricted to classical and activated monocytes.
# ============================================================

mono <- subset(
  pre_bcg,
  subset =
    cell_type %in%
      MYELOID_CELL_TYPES
)

cat(
  "Classical + Activated Monocyte cells: ",
  format(
    ncol(
      mono
    ),
    big.mark = ","
  ),
  "\n",
  sep = ""
)

cat("\nCell types:\n")

print(
  table(
    mono$cell_type
  )
)

cat("\nOutcome groups:\n")

print(
  table(
    mono$Outcome
  )
)

cat("\nSamples by outcome:\n")

print(
  mono@meta.data %>%
    distinct(
      sample_id,
      Outcome
    ) %>%
    count(
      Outcome
    )
)


# ============================================================
# 7. PREPARE RNA ASSAY
#
# Join Seurat v5 sample-specific RNA layers and recreate
# normalized expression deterministically.
# ============================================================

DefaultAssay(
  mono
) <- "RNA"

mono <- JoinLayers(
  mono,
  assay = "RNA"
)

mono <- NormalizeData(
  mono,
  assay = "RNA",
  verbose = FALSE
)


# ============================================================
# 8. VALIDATE SIGNATURE GENES
# ============================================================

genes_present <- AP1_GENES[
  AP1_GENES %in%
    rownames(
      mono
    )
]

genes_missing <- setdiff(
  AP1_GENES,
  genes_present
)

cat("\nAP-1 genes present:\n")

print(
  genes_present
)

if (
  length(
    genes_missing
  ) >
    0
) {

  stop(
    paste0(
      "Missing AP-1 signature gene(s): ",
      paste(
        genes_missing,
        collapse = ", "
      )
    )
  )

}


# ============================================================
# 9. SINGLE-CELL AP-1 MODULE SCORE
#
# Final signature:
#   JUN, JUNB, JUND, FOS, FOSB, FOSL2
#
# AddModuleScore requires a list when the six genes are
# intended to represent one expression program.
# ============================================================

cat("\n")
cat("============================================================\n")
cat("Computing AP-1 module score\n")
cat("============================================================\n")

set.seed(
  MODULE_SEED
)

mono <- AddModuleScore(
  object = mono,
  features = list(
    AP1_GENES
  ),
  name = MODULE_NAME,
  nbin = MODULE_NBIN,
  ctrl = MODULE_CTRL,
  seed = MODULE_SEED
)

score_col <- paste0(
  MODULE_NAME,
  "1"
)

if (
  !score_col %in%
    colnames(
      mono@meta.data
    )
) {

  stop(
    paste0(
      "Expected module-score column not found: ",
      score_col
    )
  )

}


# ============================================================
# 10. EXPORT CELL-LEVEL SCORES
# ============================================================

cell_scores <- FetchData(
  mono,
  vars = c(
    score_col,
    "sample_id",
    "Outcome",
    "cell_type"
  )
) %>%
  rownames_to_column(
    "cell"
  )

colnames(
  cell_scores
)[
  colnames(
    cell_scores
  ) ==
    score_col
] <- "AP1_score"

save_csv(
  cell_scores,
  file.path(
    table_dir,
    "AP1_single_cell_scores.csv"
  )
)


# ============================================================
# 11. AGGREGATE TO SAMPLE LEVEL
#
# Manuscript Methods:
#   patient/sample score = mean AP-1 module score across all
#   classical + activated monocytes within that sample.
# ============================================================

patient_scores <- cell_scores %>%
  group_by(
    sample_id,
    Outcome
  ) %>%
  summarise(
    AP1_score = mean(
      AP1_score,
      na.rm = TRUE
    ),
    n_cells = n(),
    .groups = "drop"
  )

patient_scores$Outcome <- factor(
  patient_scores$Outcome,
  levels = c(
    RESPONDER_LABEL,
    NONRESPONDER_LABEL
  )
)

patient_scores <- patient_scores %>%
  arrange(
    Outcome,
    AP1_score
  )

save_csv(
  patient_scores,
  file.path(
    table_dir,
    "AP1_patient_level_scores.csv"
  )
)

cat("\nPatient-level AP-1 scores:\n")

print(
  patient_scores
)


# ============================================================
# 12. TWO-SIDED WILCOXON RANK-SUM TEST
#
# The manuscript specifies a two-sided Wilcoxon rank-sum test.
# We report exact=TRUE when possible and also save the ordinary
# R wilcox.test result for complete provenance.
# ============================================================

rec_scores <- patient_scores$AP1_score[
  patient_scores$Outcome ==
    NONRESPONDER_LABEL
]

rem_scores <- patient_scores$AP1_score[
  patient_scores$Outcome ==
    RESPONDER_LABEL
]

wilcox_exact <- tryCatch(
  wilcox.test(
    rec_scores,
    rem_scores,
    alternative = "two.sided",
    exact = TRUE,
    conf.int = FALSE
  ),
  warning = function(w) {

    message(
      "Exact Wilcoxon test produced warning: ",
      conditionMessage(
        w
      )
    )

    wilcox.test(
      rec_scores,
      rem_scores,
      alternative = "two.sided",
      exact = FALSE,
      conf.int = FALSE
    )

  }
)

wilcox_asymptotic <- wilcox.test(
  rec_scores,
  rem_scores,
  alternative = "two.sided",
  exact = FALSE,
  conf.int = FALSE
)

rb <- rank_biserial_from_u(
  rec_scores,
  rem_scores
)

wilcox_results <- data.frame(
  n_Recurrence = length(
    rec_scores
  ),
  n_Remission = length(
    rem_scores
  ),
  median_Recurrence = median(
    rec_scores
  ),
  median_Remission = median(
    rem_scores
  ),
  mean_Recurrence = mean(
    rec_scores
  ),
  mean_Remission = mean(
    rem_scores
  ),
  W_exact_or_fallback = unname(
    wilcox_exact$statistic
  ),
  p_exact_or_fallback = wilcox_exact$p.value,
  W_asymptotic = unname(
    wilcox_asymptotic$statistic
  ),
  p_asymptotic = wilcox_asymptotic$p.value,
  Mann_Whitney_U_Recurrence = rb$U,
  rank_biserial_Recurrence_vs_Remission =
    rb$rank_biserial,
  stringsAsFactors = FALSE
)

save_csv(
  wilcox_results,
  file.path(
    table_dir,
    "AP1_Wilcoxon_and_effect_size.csv"
  )
)

cat("\n")
cat("============================================================\n")
cat("Wilcoxon test\n")
cat("============================================================\n")

print(
  wilcox_results,
  row.names = FALSE
)


# ============================================================
# 13. ROC ANALYSIS
#
# Historical / manuscript direction:
#   higher AP-1 score = recurrence.
#
# pROC coding:
#   controls = Remission
#   cases    = Recurrence
#   direction = "<"
# ============================================================

cat("\n")
cat("============================================================\n")
cat("ROC analysis\n")
cat("============================================================\n")

roc_obj <- roc(
  response = patient_scores$Outcome,
  predictor = patient_scores$AP1_score,
  levels = c(
    RESPONDER_LABEL,
    NONRESPONDER_LABEL
  ),
  direction = "<",
  quiet = TRUE
)

auc_value <- as.numeric(
  auc(
    roc_obj
  )
)

auc_ci <- as.numeric(
  ci.auc(
    roc_obj
  )
)

# Test AUC against null AUC = 0.5.
auc_test <- tryCatch(
  roc.test(
    roc_obj,
    auc = 0.5,
    method = "delong"
  ),
  error = function(e) {
    NULL
  }
)

auc_p_value <- if (
  is.null(
    auc_test
  )
) {
  NA_real_
} else {
  auc_test$p.value
}

roc_results <- data.frame(
  AUC = auc_value,
  CI95_low = auc_ci[
    1
  ],
  CI95_mid = auc_ci[
    2
  ],
  CI95_high = auc_ci[
    3
  ],
  AUC_vs_0.5_p_value = auc_p_value,
  direction = "higher score = Recurrence",
  stringsAsFactors = FALSE
)

save_csv(
  roc_results,
  file.path(
    table_dir,
    "AP1_ROC_statistics.csv"
  )
)

print(
  roc_results,
  row.names = FALSE
)


# ============================================================
# 14. ROC COORDINATES
# ============================================================

roc_coordinates <- coords(
  roc_obj,
  x = "all",
  ret = c(
    "threshold",
    "specificity",
    "sensitivity"
  ),
  transpose = FALSE
)

roc_coordinates <- as.data.frame(
  roc_coordinates
)

save_csv(
  roc_coordinates,
  file.path(
    table_dir,
    "AP1_ROC_coordinates.csv"
  )
)


# ============================================================
# 15. FIGURE 5A-LIKE PATIENT-LEVEL SCORE PLOT
# ============================================================

plot_scores <- patient_scores

plot_scores$Outcome_display <- ifelse(
  plot_scores$Outcome ==
    RESPONDER_LABEL,
  "Responders",
  "Non-responders"
)

plot_scores$Outcome_display <- factor(
  plot_scores$Outcome_display,
  levels = c(
    "Responders",
    "Non-responders"
  )
)

p_score <- ggplot(
  plot_scores,
  aes(
    x = Outcome_display,
    y = AP1_score,
    fill = Outcome_display
  )
) +
  geom_violin(
    trim = FALSE,
    alpha = 0.4,
    color = NA
  ) +
  geom_boxplot(
    width = 0.2,
    outlier.shape = NA
  ) +
  geom_jitter(
    width = 0.1,
    size = 2.2,
    alpha = 0.8
  ) +
  theme_classic(
    base_size = 14
  ) +
  labs(
    x = NULL,
    y = "AP-1 Monocyte Signature Score",
    title = "Pre-BCG AP-1 Signature by Outcome",
    subtitle = paste0(
      "Two-sided Wilcoxon p = ",
      signif(
        wilcox_exact$p.value,
        4
      )
    )
  ) +
  theme(
    legend.position = "none",
    plot.title = element_text(
      hjust = 0.5,
      face = "bold"
    ),
    plot.subtitle = element_text(
      hjust = 0.5
    )
  )

ggsave(
  filename = file.path(
    figure_dir,
    "fig5a_AP1_patient_scores.pdf"
  ),
  plot = p_score,
  width = 6,
  height = 5
)

ggsave(
  filename = file.path(
    figure_dir,
    "fig5a_AP1_patient_scores.png"
  ),
  plot = p_score,
  width = 6,
  height = 5,
  dpi = 600
)


# ============================================================
# 16. FIGURE 5C-LIKE ROC
# ============================================================

pdf(
  file.path(
    figure_dir,
    "fig5c_AP1_ROC.pdf"
  ),
  width = 6,
  height = 6
)

plot(
  roc_obj,
  lwd = 3,
  legacy.axes = TRUE,
  main = "Pre-BCG AP-1 Signature ROC"
)

abline(
  a = 0,
  b = 1,
  lty = 2
)

legend(
  "bottomright",
  legend = c(
    paste0(
      "AUC = ",
      sprintf(
        "%.3f",
        auc_value
      )
    ),
    paste0(
      "95% CI ",
      sprintf(
        "%.3f",
        auc_ci[
          1
        ]
      ),
      "-",
      sprintf(
        "%.3f",
        auc_ci[
          3
        ]
      )
    )
  ),
  bty = "n"
)

dev.off()


png(
  file.path(
    figure_dir,
    "fig5c_AP1_ROC.png"
  ),
  width = 1800,
  height = 1800,
  res = 300
)

plot(
  roc_obj,
  lwd = 3,
  legacy.axes = TRUE,
  main = "Pre-BCG AP-1 Signature ROC"
)

abline(
  a = 0,
  b = 1,
  lty = 2
)

legend(
  "bottomright",
  legend = c(
    paste0(
      "AUC = ",
      sprintf(
        "%.3f",
        auc_value
      )
    ),
    paste0(
      "95% CI ",
      sprintf(
        "%.3f",
        auc_ci[
          1
        ]
      ),
      "-",
      sprintf(
        "%.3f",
        auc_ci[
          3
        ]
      )
    )
  ),
  bty = "n"
)

dev.off()


# ============================================================
# 17. INDIVIDUAL AP-1 GENE EXPRESSION
#
# Aggregate normalized expression to the sample level.
# ============================================================

gene_expression_tables <- list()

for (
  gene in AP1_GENES
) {

  df_gene <- FetchData(
    mono,
    vars = c(
      gene,
      "sample_id",
      "Outcome"
    )
  )

  colnames(
    df_gene
  )[
    1
  ] <- "expression"

  patient_gene <- df_gene %>%
    group_by(
      sample_id,
      Outcome
    ) %>%
    summarise(
      mean_expression = mean(
        expression,
        na.rm = TRUE
      ),
      .groups = "drop"
    ) %>%
    mutate(
      gene = gene
    )

  gene_expression_tables[
    [
      gene
    ]
  ] <- patient_gene

}

individual_gene_expression <- bind_rows(
  gene_expression_tables
)

save_csv(
  individual_gene_expression,
  file.path(
    table_dir,
    "AP1_individual_gene_patient_expression.csv"
  )
)


# ============================================================
# 18. FRACTION OF CELLS EXPRESSING EACH AP-1 GENE
#
# Figure 5b is described as expression of individual AP-1
# components shown as fraction of expressing cells.
#
# A normalized expression value > 0 is treated as detected.
# ============================================================

gene_fraction_tables <- list()

for (
  gene in AP1_GENES
) {

  df_gene <- FetchData(
    mono,
    vars = c(
      gene,
      "sample_id",
      "Outcome"
    )
  )

  colnames(
    df_gene
  )[
    1
  ] <- "expression"

  patient_fraction <- df_gene %>%
    group_by(
      sample_id,
      Outcome
    ) %>%
    summarise(
      fraction_expressing = mean(
        expression >
          0,
        na.rm = TRUE
      ),
      .groups = "drop"
    ) %>%
    mutate(
      gene = gene
    )

  gene_fraction_tables[
    [
      gene
    ]
  ] <- patient_fraction

}

individual_gene_fraction <- bind_rows(
  gene_fraction_tables
)

save_csv(
  individual_gene_fraction,
  file.path(
    table_dir,
    "AP1_individual_gene_fraction_expressing.csv"
  )
)


# ============================================================
# 19. FIGURE 5B-LIKE FRACTION-EXPRESSING PLOT
# ============================================================

fraction_plot_df <- individual_gene_fraction %>%
  mutate(
    Outcome_display = ifelse(
      Outcome ==
        RESPONDER_LABEL,
      "Responders",
      "Non-responders"
    ),
    Outcome_display = factor(
      Outcome_display,
      levels = c(
        "Responders",
        "Non-responders"
      )
    ),
    gene = factor(
      gene,
      levels = AP1_GENES
    )
  )

p_fraction <- ggplot(
  fraction_plot_df,
  aes(
    x = Outcome_display,
    y = fraction_expressing,
    fill = Outcome_display
  )
) +
  geom_boxplot(
    width = 0.55,
    outlier.shape = NA,
    alpha = 0.5
  ) +
  geom_jitter(
    width = 0.08,
    size = 1.7,
    alpha = 0.75
  ) +
  facet_wrap(
    ~ gene,
    ncol = 3
  ) +
  theme_classic(
    base_size = 11
  ) +
  labs(
    x = NULL,
    y = "Fraction of expressing cells",
    title = "AP-1 Component Expression by Outcome"
  ) +
  theme(
    legend.position = "none",
    strip.text = element_text(
      face = "bold"
    )
  )

ggsave(
  filename = file.path(
    figure_dir,
    "fig5b_AP1_fraction_expressing.pdf"
  ),
  plot = p_fraction,
  width = 8,
  height = 6
)

ggsave(
  filename = file.path(
    figure_dir,
    "fig5b_AP1_fraction_expressing.png"
  ),
  plot = p_fraction,
  width = 8,
  height = 6,
  dpi = 600
)


# ============================================================
# 20. SAVE ANALYSIS OBJECT
# ============================================================

saveRDS(
  mono,
  file.path(
    outdir,
    "pre_bcg_classical_activated_monocytes_AP1_scored.rds"
  )
)


# ============================================================
# 21. PARAMETERS
# ============================================================

parameter_table <- data.frame(
  parameter = c(
    "input_rds",
    "cell_types",
    "AP1_genes",
    "module_scoring",
    "module_seed",
    "module_nbin",
    "module_ctrl",
    "sample_aggregation",
    "outcome_reference",
    "outcome_case",
    "Wilcoxon_alternative",
    "ROC_direction"
  ),
  value = c(
    input_rds,
    paste(
      MYELOID_CELL_TYPES,
      collapse = "; "
    ),
    paste(
      AP1_GENES,
      collapse = "; "
    ),
    "Seurat::AddModuleScore",
    MODULE_SEED,
    MODULE_NBIN,
    MODULE_CTRL,
    "mean module score per sample",
    RESPONDER_LABEL,
    NONRESPONDER_LABEL,
    "two.sided",
    "higher score = Recurrence"
  ),
  stringsAsFactors = FALSE
)

save_csv(
  parameter_table,
  file.path(
    outdir,
    "08_ap1_signature_parameters.csv"
  )
)


# ============================================================
# 22. SESSION INFO
# ============================================================

sink(
  file.path(
    outdir,
    "08_ap1_signature_sessionInfo.txt"
  )
)

print(
  sessionInfo()
)

sink()


# ============================================================
# 23. FINAL REPORT
# ============================================================

cat("\n")
cat("============================================================\n")
cat("PRE-BCG AP-1 SIGNATURE ANALYSIS COMPLETE\n")
cat("============================================================\n")

cat(
  "\nSignature genes: ",
  paste(
    AP1_GENES,
    collapse = ", "
  ),
  "\n",
  sep = ""
)

cat(
  "Cells analysed: ",
  format(
    ncol(
      mono
    ),
    big.mark = ","
  ),
  "\n",
  sep = ""
)

cat(
  "Samples analysed: ",
  nrow(
    patient_scores
  ),
  "\n",
  sep = ""
)

cat(
  "Responder samples: ",
  sum(
    patient_scores$Outcome ==
      RESPONDER_LABEL
  ),
  "\n",
  sep = ""
)

cat(
  "Recurrence samples: ",
  sum(
    patient_scores$Outcome ==
      NONRESPONDER_LABEL
  ),
  "\n",
  sep = ""
)

cat(
  "\nTwo-sided Wilcoxon p-value: ",
  signif(
    wilcox_exact$p.value,
    6
  ),
  "\n",
  sep = ""
)

cat(
  "Rank-biserial correlation: ",
  round(
    rb$rank_biserial,
    4
  ),
  "\n",
  sep = ""
)

cat(
  "ROC AUC: ",
  round(
    auc_value,
    4
  ),
  "\n",
  sep = ""
)

cat(
  "AUC 95% CI: ",
  paste(
    round(
      c(
        auc_ci[
          1
        ],
        auc_ci[
          3
        ]
      ),
      4
    ),
    collapse = " - "
  ),
  "\n",
  sep = ""
)

cat(
  "AUC vs 0.5 p-value: ",
  signif(
    auc_p_value,
    6
  ),
  "\n",
  sep = ""
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
