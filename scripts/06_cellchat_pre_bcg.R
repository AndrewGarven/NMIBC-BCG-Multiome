# ============================================================
# 06_cellchat_pre_bcg.R
#
# Reproducible pre-BCG cell-cell communication analysis
# corresponding to Figure 3a-d.
#
# INPUT:
#   Historical finalized pre-BCG Seurat/Signac object:
#       G:/merged_obj.rds
#
# GROUPING:
#   Major annotated cell types only:
#       group.by = "cell_type"
#
# OUTCOME:
#   Historical labels:
#       Remission  = responder
#       Recurrence = non-responder / early recurrence
#
# HISTORICAL CELLCHAT SETTINGS RECOVERED FROM ORIGINAL CODE:
#   - CellChatDB.human
#   - database restricted to:
#       "Secreted Signaling"
#       "Cell-Cell Contact"
#   - identifyOverExpressedGenes()
#   - identifyOverExpressedInteractions()
#   - computeCommunProb(population.size = TRUE)
#   - filterCommunication(min.cells = 10)
#   - computeCommunProbPathway()
#   - aggregateNet()
#
# FIGURE-LEVEL OUTPUTS:
#   Fig. 3a: communication networks by outcome
#   Fig. 3b: pathway information-flow comparison
#   Fig. 3c: TGFb communication difference heatmap
#   Fig. 3d: CD86 communication difference heatmap
#
# IMPORTANT:
# This script intentionally DOES NOT use the later NRup_DA
# neighbourhood-derived Classical Monocyte split. The manuscript
# Figure 3 analysis is reconstructed using the major cell_type
# annotations only.
#
# Run from repository root with R 4.4.1:
#
# & "C:\Users\andrew\AppData\Local\Programs\R\R-4.4.1\bin\Rscript.exe" `
#     scripts\06_cellchat_pre_bcg.R
#
# ============================================================


# ============================================================
# 1. PACKAGES
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(CellChat)
  library(patchwork)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(future)
})


# ============================================================
# 2. PATHS
# ============================================================

input_rds <- "G:/merged_obj.rds"

outdir <- file.path(
  "results",
  "06_cellchat_pre_bcg"
)

figure_dir <- file.path(
  outdir,
  "figures"
)

table_dir <- file.path(
  outdir,
  "tables"
)

object_dir <- file.path(
  outdir,
  "objects"
)

dir.create(
  outdir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  figure_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  table_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  object_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ============================================================
# 3. PARAMETERS
# ============================================================

GROUP_COLUMN <- "cell_type"
OUTCOME_COLUMN <- "Outcome"

RESPONDER_LABEL <- "Remission"
NONRESPONDER_LABEL <- "Recurrence"

MIN_CELLS <- 10
POPULATION_SIZE <- TRUE

CELLCHAT_DATABASE_CATEGORIES <- c(
  "Secreted Signaling",
  "Cell-Cell Contact"
)

N_WORKERS <- 4

# Small pseudocount used only for log2 ratio summaries.
PSEUDOCOUNT <- 1e-9


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


safe_sum <- function(
  x
) {

  sum(
    x,
    na.rm = TRUE
  )

}


get_network_summary <- function(
  cellchat_obj,
  label
) {

  count_mat <- cellchat_obj@net$count
  weight_mat <- cellchat_obj@net$weight

  data.frame(
    group = label,
    n_cell_groups = length(
      levels(
        cellchat_obj@idents
      )
    ),
    n_inferred_interactions = safe_sum(
      count_mat
    ),
    total_interaction_strength = safe_sum(
      weight_mat
    ),
    stringsAsFactors = FALSE
  )

}


make_pair_comparison <- function(
  df_resp,
  df_nonresp
) {

  full_join(
    df_resp %>%
      group_by(
        source,
        target
      ) %>%
      summarise(
        prob_responder = sum(
          prob,
          na.rm = TRUE
        ),
        .groups = "drop"
      ),
    df_nonresp %>%
      group_by(
        source,
        target
      ) %>%
      summarise(
        prob_nonresponder = sum(
          prob,
          na.rm = TRUE
        ),
        .groups = "drop"
      ),
    by = c(
      "source",
      "target"
    )
  ) %>%
    mutate(
      prob_responder = replace_na(
        prob_responder,
        0
      ),
      prob_nonresponder = replace_na(
        prob_nonresponder,
        0
      ),
      difference_nonresponder_minus_responder =
        prob_nonresponder - prob_responder,
      log2FC_nonresponder_vs_responder =
        log2(
          (
            prob_nonresponder +
              PSEUDOCOUNT
          ) /
            (
              prob_responder +
                PSEUDOCOUNT
            )
        )
    ) %>%
    arrange(
      desc(
        difference_nonresponder_minus_responder
      )
    )

}


make_pathway_comparison <- function(
  df_resp,
  df_nonresp
) {

  full_join(
    df_resp %>%
      group_by(
        pathway_name
      ) %>%
      summarise(
        information_flow_responder = sum(
          prob,
          na.rm = TRUE
        ),
        .groups = "drop"
      ),
    df_nonresp %>%
      group_by(
        pathway_name
      ) %>%
      summarise(
        information_flow_nonresponder = sum(
          prob,
          na.rm = TRUE
        ),
        .groups = "drop"
      ),
    by = "pathway_name"
  ) %>%
    mutate(
      information_flow_responder = replace_na(
        information_flow_responder,
        0
      ),
      information_flow_nonresponder = replace_na(
        information_flow_nonresponder,
        0
      ),
      difference_nonresponder_minus_responder =
        information_flow_nonresponder -
        information_flow_responder,
      log2FC_nonresponder_vs_responder =
        log2(
          (
            information_flow_nonresponder +
              PSEUDOCOUNT
          ) /
            (
              information_flow_responder +
                PSEUDOCOUNT
            )
        )
    ) %>%
    arrange(
      desc(
        difference_nonresponder_minus_responder
      )
    )

}


make_lr_comparison <- function(
  df_resp,
  df_nonresp
) {

  join_columns <- c(
    "source",
    "target",
    "pathway_name",
    "ligand",
    "receptor"
  )

  full_join(
    df_resp %>%
      group_by(
        across(
          all_of(
            join_columns
          )
        )
      ) %>%
      summarise(
        prob_responder = sum(
          prob,
          na.rm = TRUE
        ),
        .groups = "drop"
      ),
    df_nonresp %>%
      group_by(
        across(
          all_of(
            join_columns
          )
        )
      ) %>%
      summarise(
        prob_nonresponder = sum(
          prob,
          na.rm = TRUE
        ),
        .groups = "drop"
      ),
    by = join_columns
  ) %>%
    mutate(
      prob_responder = replace_na(
        prob_responder,
        0
      ),
      prob_nonresponder = replace_na(
        prob_nonresponder,
        0
      ),
      difference_nonresponder_minus_responder =
        prob_nonresponder -
        prob_responder,
      log2FC_nonresponder_vs_responder =
        log2(
          (
            prob_nonresponder +
              PSEUDOCOUNT
          ) /
            (
              prob_responder +
                PSEUDOCOUNT
            )
        )
    ) %>%
    arrange(
      desc(
        difference_nonresponder_minus_responder
      )
    )

}


make_pathway_matrix_df <- function(
  communication_df,
  pathway,
  group_name,
  sources_use,
  targets_use
) {

  observed <- communication_df %>%
    filter(
      pathway_name == pathway,
      source %in% sources_use,
      target %in% targets_use
    ) %>%
    group_by(
      source,
      target
    ) %>%
    summarise(
      prob = sum(
        prob,
        na.rm = TRUE
      ),
      .groups = "drop"
    )

  full_grid <- expand_grid(
    source = sources_use,
    target = targets_use
  )

  full_grid %>%
    left_join(
      observed,
      by = c(
        "source",
        "target"
      )
    ) %>%
    mutate(
      prob = replace_na(
        prob,
        0
      ),
      group = group_name
    )

}


make_pathway_difference_df <- function(
  df_resp,
  df_nonresp,
  pathway,
  shared_cell_types
) {

  responder <- make_pathway_matrix_df(
    communication_df = df_resp,
    pathway = pathway,
    group_name = "Responders",
    sources_use = shared_cell_types,
    targets_use = shared_cell_types
  )

  nonresponder <- make_pathway_matrix_df(
    communication_df = df_nonresp,
    pathway = pathway,
    group_name = "Non-responders",
    sources_use = shared_cell_types,
    targets_use = shared_cell_types
  )

  responder %>%
    select(
      source,
      target,
      prob_responder = prob
    ) %>%
    left_join(
      nonresponder %>%
        select(
          source,
          target,
          prob_nonresponder = prob
        ),
      by = c(
        "source",
        "target"
      )
    ) %>%
    mutate(
      difference_nonresponder_minus_responder =
        prob_nonresponder -
        prob_responder,
      log2FC_nonresponder_vs_responder =
        log2(
          (
            prob_nonresponder +
              PSEUDOCOUNT
          ) /
            (
              prob_responder +
                PSEUDOCOUNT
            )
        )
    )

}


plot_difference_heatmap <- function(
  difference_df,
  pathway,
  output_pdf,
  output_png
) {

  max_abs <- max(
    abs(
      difference_df$difference_nonresponder_minus_responder
    ),
    na.rm = TRUE
  )

  if (
    !is.finite(
      max_abs
    ) ||
      max_abs == 0
  ) {
    max_abs <- 1
  }

  p <- ggplot(
    difference_df,
    aes(
      x = target,
      y = source,
      fill = difference_nonresponder_minus_responder
    )
  ) +
    geom_tile(
      color = "white"
    ) +
    scale_fill_gradient2(
      low = "#2166ac",
      mid = "white",
      high = "#b2182b",
      midpoint = 0,
      limits = c(
        -max_abs,
        max_abs
      )
    ) +
    theme_classic(
      base_size = 11
    ) +
    theme(
      axis.text.x = element_text(
        angle = 45,
        hjust = 1
      )
    ) +
    labs(
      title = paste0(
        pathway,
        " signalling"
      ),
      subtitle = "Non-responder minus responder communication probability",
      x = "Target cell type",
      y = "Source cell type",
      fill = "Difference"
    )

  ggsave(
    filename = output_pdf,
    plot = p,
    width = 9,
    height = 7,
    units = "in"
  )

  ggsave(
    filename = output_png,
    plot = p,
    width = 9,
    height = 7,
    units = "in",
    dpi = 600
  )

  p

}


plot_cellchat_network_pair <- function(
  cellchat_obj,
  title_prefix,
  output_pdf,
  output_png
) {

  group_tab <- table(
    cellchat_obj@idents
  )

  group_size <- as.numeric(
    group_tab
  )

  count_mat <- cellchat_obj@net$count
  weight_mat <- cellchat_obj@net$weight

  pdf(
    output_pdf,
    width = 14,
    height = 7
  )

  old_par <- par(
    no.readonly = TRUE
  )

  par(
    mfrow = c(
      1,
      2
    ),
    xpd = TRUE,
    mar = c(
      1,
      1,
      3,
      1
    )
  )

  netVisual_circle(
    count_mat,
    vertex.weight = sqrt(
      group_size
    ),
    weight.scale = TRUE,
    label.edge = FALSE,
    title.name = paste0(
      title_prefix,
      "\nNumber of interactions"
    ),
    vertex.label.cex = 0.7
  )

  netVisual_circle(
    weight_mat,
    vertex.weight = sqrt(
      group_size
    ),
    weight.scale = TRUE,
    label.edge = FALSE,
    title.name = paste0(
      title_prefix,
      "\nInteraction strength"
    ),
    vertex.label.cex = 0.7
  )

  par(
    old_par
  )

  dev.off()

  png(
    output_png,
    width = 4200,
    height = 2100,
    res = 300
  )

  old_par <- par(
    no.readonly = TRUE
  )

  par(
    mfrow = c(
      1,
      2
    ),
    xpd = TRUE,
    mar = c(
      1,
      1,
      3,
      1
    )
  )

  netVisual_circle(
    count_mat,
    vertex.weight = sqrt(
      group_size
    ),
    weight.scale = TRUE,
    label.edge = FALSE,
    title.name = paste0(
      title_prefix,
      "\nNumber of interactions"
    ),
    vertex.label.cex = 0.7
  )

  netVisual_circle(
    weight_mat,
    vertex.weight = sqrt(
      group_size
    ),
    weight.scale = TRUE,
    label.edge = FALSE,
    title.name = paste0(
      title_prefix,
      "\nInteraction strength"
    ),
    vertex.label.cex = 0.7
  )

  par(
    old_par
  )

  dev.off()

}


# ============================================================
# 5. LOAD PRE-BCG OBJECT
# ============================================================

cat("\n")
cat("============================================================\n")
cat("Loading pre-BCG Seurat object\n")
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
  GROUP_COLUMN,
  OUTCOME_COLUMN
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
  ) > 0
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
  !"RNA" %in% Assays(
    pre_bcg
  )
) {
  stop(
    "RNA assay not found."
  )
}

cat(
  "Cells: ",
  format(
    ncol(
      pre_bcg
    ),
    big.mark = ","
  ),
  "\n",
  sep = ""
)

cat("\nOutcome counts:\n")

print(
  table(
    pre_bcg[[OUTCOME_COLUMN]][, 1],
    useNA = "ifany"
  )
)

cat("\nCell-type counts:\n")

print(
  sort(
    table(
      pre_bcg[[GROUP_COLUMN]][, 1]
    ),
    decreasing = TRUE
  )
)


# ============================================================
# 6. PREPARE RNA ASSAY
#
# CellChat expects normalized expression.
#
# The original CellChat code used the RNA assay directly.
# To make the reviewer-facing script robust to Seurat v5
# split layers in merged_obj.rds, we join RNA layers and run
# standard NormalizeData() before creating CellChat objects.
#
# This does not use the ATAC assay.
# ============================================================

cat("\n")
cat("============================================================\n")
cat("Preparing RNA expression for CellChat\n")
cat("============================================================\n")

DefaultAssay(
  pre_bcg
) <- "RNA"

pre_bcg <- JoinLayers(
  pre_bcg,
  assay = "RNA"
)

pre_bcg <- NormalizeData(
  pre_bcg,
  assay = "RNA",
  verbose = FALSE
)

cat("RNA layers:\n")

print(
  Layers(
    pre_bcg[["RNA"]]
  )
)


# ============================================================
# 7. SPLIT BY PRE-BCG OUTCOME
# ============================================================

outcome_values <- as.character(
  pre_bcg[[OUTCOME_COLUMN]][, 1]
)

if (
  !all(
    c(
      RESPONDER_LABEL,
      NONRESPONDER_LABEL
    ) %in% outcome_values
  )
) {
  stop(
    paste0(
      "Expected Outcome labels '",
      RESPONDER_LABEL,
      "' and '",
      NONRESPONDER_LABEL,
      "'."
    )
  )
}

responder_cells <- colnames(
  pre_bcg
)[
  outcome_values ==
    RESPONDER_LABEL
]

nonresponder_cells <- colnames(
  pre_bcg
)[
  outcome_values ==
    NONRESPONDER_LABEL
]

responders <- subset(
  pre_bcg,
  cells = responder_cells
)

nonresponders <- subset(
  pre_bcg,
  cells = nonresponder_cells
)


cat(
  "\nResponder cells: ",
  format(
    ncol(
      responders
    ),
    big.mark = ","
  ),
  "\n",
  sep = ""
)

cat(
  "Non-responder cells: ",
  format(
    ncol(
      nonresponders
    ),
    big.mark = ","
  ),
  "\n",
  sep = ""
)


# ============================================================
# 8. CREATE CELLCHAT OBJECTS
#
# Manuscript analysis:
#   all major annotated cell types
#   group.by = "cell_type"
# ============================================================

cat("\n")
cat("============================================================\n")
cat("Creating CellChat objects\n")
cat("============================================================\n")

cellchat_responder <- createCellChat(
  object = responders,
  group.by = GROUP_COLUMN,
  assay = "RNA"
)

cellchat_nonresponder <- createCellChat(
  object = nonresponders,
  group.by = GROUP_COLUMN,
  assay = "RNA"
)


# ============================================================
# 9. CELLCHAT DATABASE
#
# Exact recovered historical categories:
#   Secreted Signaling
#   Cell-Cell Contact
# ============================================================

CellChatDB <- CellChatDB.human

CellChatDB.use <- subsetDB(
  CellChatDB,
  search = CELLCHAT_DATABASE_CATEGORIES,
  key = "annotation"
)

cellchat_responder@DB <- CellChatDB.use
cellchat_nonresponder@DB <- CellChatDB.use


# ============================================================
# 10. SUBSET SIGNALING DATA
# ============================================================

cellchat_responder <- subsetData(
  cellchat_responder
)

cellchat_nonresponder <- subsetData(
  cellchat_nonresponder
)


# ============================================================
# 11. OVEREXPRESSED GENES / INTERACTIONS
#
# Recovered historical parallel setting:
#   future.globals.maxSize = 4 GB
#   multisession workers = 4
# ============================================================

cat("\n")
cat("============================================================\n")
cat("Identifying overexpressed genes and interactions\n")
cat("============================================================\n")

options(
  future.globals.maxSize = 8 * 1024^3
)

plan(
  sequential
)

cellchat_responder <- identifyOverExpressedGenes(
  cellchat_responder
)

cellchat_responder <- identifyOverExpressedInteractions(
  cellchat_responder
)

cellchat_nonresponder <- identifyOverExpressedGenes(
  cellchat_nonresponder
)

cellchat_nonresponder <- identifyOverExpressedInteractions(
  cellchat_nonresponder
)


# ============================================================
# 12. COMMUNICATION PROBABILITIES
#
# Historical exact parameter:
#   population.size = TRUE
# ============================================================

cat("\n")
cat("============================================================\n")
cat("Computing communication probabilities\n")
cat("============================================================\n")

cellchat_responder <- computeCommunProb(
  cellchat_responder,
  population.size = POPULATION_SIZE
)

cellchat_nonresponder <- computeCommunProb(
  cellchat_nonresponder,
  population.size = POPULATION_SIZE
)


# ============================================================
# 13. MINIMUM CELL FILTER
#
# Historical exact parameter:
#   min.cells = 10
# ============================================================

cellchat_responder <- filterCommunication(
  cellchat_responder,
  min.cells = MIN_CELLS
)

cellchat_nonresponder <- filterCommunication(
  cellchat_nonresponder,
  min.cells = MIN_CELLS
)


# ============================================================
# 14. PATHWAY PROBABILITIES + NETWORK AGGREGATION
# ============================================================

cellchat_responder <- computeCommunProbPathway(
  cellchat_responder
)

cellchat_nonresponder <- computeCommunProbPathway(
  cellchat_nonresponder
)

cellchat_responder <- aggregateNet(
  cellchat_responder
)

cellchat_nonresponder <- aggregateNet(
  cellchat_nonresponder
)

plan(
  sequential
)


# ============================================================
# 15. SAVE CELLCHAT OBJECTS
# ============================================================

saveRDS(
  cellchat_responder,
  file = file.path(
    object_dir,
    "cellchat_pre_bcg_responders.rds"
  )
)

saveRDS(
  cellchat_nonresponder,
  file = file.path(
    object_dir,
    "cellchat_pre_bcg_nonresponders.rds"
  )
)


# ============================================================
# 16. EXTRACT ALL INFERRED COMMUNICATIONS
# ============================================================

df_responder <- subsetCommunication(
  cellchat_responder
)

df_nonresponder <- subsetCommunication(
  cellchat_nonresponder
)

df_responder$outcome_group <- "Responder"

df_nonresponder$outcome_group <- "Non-responder"

save_csv(
  df_responder,
  file.path(
    table_dir,
    "all_interactions_responders.csv"
  )
)

save_csv(
  df_nonresponder,
  file.path(
    table_dir,
    "all_interactions_nonresponders.csv"
  )
)

save_csv(
  bind_rows(
    df_responder,
    df_nonresponder
  ),
  file.path(
    table_dir,
    "all_interactions_combined.csv"
  )
)


# ============================================================
# 17. GLOBAL NETWORK SUMMARY
#
# Useful audit for manuscript statement that pre-BCG
# non-responders show a greater number of inferred
# ligand-receptor interactions.
# ============================================================

network_summary <- bind_rows(
  get_network_summary(
    cellchat_responder,
    "Responder"
  ),
  get_network_summary(
    cellchat_nonresponder,
    "Non-responder"
  )
)

save_csv(
  network_summary,
  file.path(
    table_dir,
    "global_network_summary.csv"
  )
)

cat("\n")
cat("============================================================\n")
cat("Global network summary\n")
cat("============================================================\n")

print(
  network_summary,
  row.names = FALSE
)


# ============================================================
# 18. PAIRWISE SOURCE -> TARGET COMPARISON
# ============================================================

pair_comparison <- make_pair_comparison(
  df_resp = df_responder,
  df_nonresp = df_nonresponder
)

save_csv(
  pair_comparison,
  file.path(
    table_dir,
    "pairwise_source_target_comparison.csv"
  )
)


# ============================================================
# 19. PATHWAY INFORMATION FLOW COMPARISON
# ============================================================

pathway_comparison <- make_pathway_comparison(
  df_resp = df_responder,
  df_nonresp = df_nonresponder
)

save_csv(
  pathway_comparison,
  file.path(
    table_dir,
    "pathway_information_flow_comparison.csv"
  )
)


# ============================================================
# 20. LIGAND-RECEPTOR COMPARISON
# ============================================================

lr_comparison <- make_lr_comparison(
  df_resp = df_responder,
  df_nonresp = df_nonresponder
)

save_csv(
  lr_comparison,
  file.path(
    table_dir,
    "ligand_receptor_comparison.csv"
  )
)


# ============================================================
# 21. COMMON CELL TYPES
#
# Use exactly the cell groups retained in BOTH CellChat
# networks for differential heatmaps.
# ============================================================

responder_cell_types <- levels(
  cellchat_responder@idents
)

nonresponder_cell_types <- levels(
  cellchat_nonresponder@idents
)

shared_cell_types <- intersect(
  responder_cell_types,
  nonresponder_cell_types
)

cat("\nShared CellChat cell types:\n")

print(
  shared_cell_types
)


# ============================================================
# 22. FIGURE 3A NETWORK PLOTS
# ============================================================

plot_cellchat_network_pair(
  cellchat_obj = cellchat_responder,
  title_prefix = "Pre-BCG Responders",
  output_pdf = file.path(
    figure_dir,
    "fig3a_pre_bcg_responders_network.pdf"
  ),
  output_png = file.path(
    figure_dir,
    "fig3a_pre_bcg_responders_network.png"
  )
)

plot_cellchat_network_pair(
  cellchat_obj = cellchat_nonresponder,
  title_prefix = "Pre-BCG Non-responders",
  output_pdf = file.path(
    figure_dir,
    "fig3a_pre_bcg_nonresponders_network.pdf"
  ),
  output_png = file.path(
    figure_dir,
    "fig3a_pre_bcg_nonresponders_network.png"
  )
)


# ============================================================
# 23. MERGED CELLCHAT OBJECT
# ============================================================

object_list <- list(
  Responders = cellchat_responder,
  Non_responders = cellchat_nonresponder
)

cellchat_merged <- mergeCellChat(
  object_list,
  add.names = names(
    object_list
  )
)

saveRDS(
  cellchat_merged,
  file = file.path(
    object_dir,
    "cellchat_pre_bcg_merged.rds"
  )
)


# ============================================================
# 24. FIGURE 3A COMPARISON BAR PLOTS
#
# Number and strength of interactions.
# ============================================================

pdf(
  file.path(
    figure_dir,
    "fig3a_compare_interactions.pdf"
  ),
  width = 10,
  height = 5
)

par(
  mfrow = c(
    1,
    2
  )
)

compareInteractions(
  cellchat_merged,
  show.legend = FALSE,
  group = c(
    1,
    2
  )
)

compareInteractions(
  cellchat_merged,
  show.legend = FALSE,
  group = c(
    1,
    2
  ),
  measure = "weight"
)

dev.off()


png(
  file.path(
    figure_dir,
    "fig3a_compare_interactions.png"
  ),
  width = 3000,
  height = 1500,
  res = 300
)

par(
  mfrow = c(
    1,
    2
  )
)

compareInteractions(
  cellchat_merged,
  show.legend = FALSE,
  group = c(
    1,
    2
  )
)

compareInteractions(
  cellchat_merged,
  show.legend = FALSE,
  group = c(
    1,
    2
  ),
  measure = "weight"
)

dev.off()


# ============================================================
# 25. DIFFERENTIAL NETWORK PLOTS
# ============================================================

pdf(
  file.path(
    figure_dir,
    "pre_bcg_differential_network.pdf"
  ),
  width = 12,
  height = 6
)

par(
  mfrow = c(
    1,
    2
  )
)

netVisual_diffInteraction(
  cellchat_merged,
  weight.scale = TRUE
)

netVisual_diffInteraction(
  cellchat_merged,
  weight.scale = TRUE,
  measure = "weight"
)

dev.off()


# ============================================================
# 26. FIGURE 3B PATHWAY INFORMATION FLOW
#
# Historical downstream function:
#   rankNet(
#       cellchat_merged,
#       mode = "comparison",
#       stacked = TRUE,
#       do.stat = TRUE
#   )
# ============================================================

pdf(
  file.path(
    figure_dir,
    "fig3b_pre_bcg_pathway_information_flow.pdf"
  ),
  width = 10,
  height = 8
)

rankNet(
  cellchat_merged,
  mode = "comparison",
  stacked = TRUE,
  do.stat = TRUE
)

dev.off()


png(
  file.path(
    figure_dir,
    "fig3b_pre_bcg_pathway_information_flow.png"
  ),
  width = 3000,
  height = 2400,
  res = 300
)

rankNet(
  cellchat_merged,
  mode = "comparison",
  stacked = TRUE,
  do.stat = TRUE
)

dev.off()


# ============================================================
# 27. PATHWAYS REPORTED IN PRE-BCG FIGURE
#
# Manuscript Figure 3 emphasizes:
#   TGFb
#   CD86
#
# Also export MHC-II because it is discussed in the pre-BCG
# Results section.
# ============================================================

pathways_of_interest <- c(
  "TGFb",
  "CD86",
  "MHC-II"
)

available_pathways <- union(
  unique(
    df_responder$pathway_name
  ),
  unique(
    df_nonresponder$pathway_name
  )
)

cat("\nAvailable pathways of interest:\n")

print(
  pathways_of_interest[
    pathways_of_interest %in%
      available_pathways
  ]
)


for (
  pathway in pathways_of_interest
) {

  pathway_interactions <- bind_rows(
    df_responder %>%
      filter(
        pathway_name == pathway
      ),
    df_nonresponder %>%
      filter(
        pathway_name == pathway
      )
  )

  save_csv(
    pathway_interactions,
    file.path(
      table_dir,
      paste0(
        gsub(
          "[^A-Za-z0-9_-]+",
          "_",
          pathway
        ),
        "_interactions.csv"
      )
    )
  )

}


# ============================================================
# 28. FIGURE 3C - TGFb DIFFERENTIAL HEATMAP
# ============================================================

if (
  "TGFb" %in%
    available_pathways
) {

  tgfb_difference <- make_pathway_difference_df(
    df_resp = df_responder,
    df_nonresp = df_nonresponder,
    pathway = "TGFb",
    shared_cell_types = shared_cell_types
  )

  save_csv(
    tgfb_difference,
    file.path(
      table_dir,
      "TGFb_source_target_difference.csv"
    )
  )

  plot_difference_heatmap(
    difference_df = tgfb_difference,
    pathway = "TGFb",
    output_pdf = file.path(
      figure_dir,
      "fig3c_pre_bcg_TGFb_difference_heatmap.pdf"
    ),
    output_png = file.path(
      figure_dir,
      "fig3c_pre_bcg_TGFb_difference_heatmap.png"
    )
  )

} else {

  warning(
    "TGFb pathway was not detected in either CellChat object."
  )

}


# ============================================================
# 29. FIGURE 3D - CD86 DIFFERENTIAL HEATMAP
# ============================================================

if (
  "CD86" %in%
    available_pathways
) {

  cd86_difference <- make_pathway_difference_df(
    df_resp = df_responder,
    df_nonresp = df_nonresponder,
    pathway = "CD86",
    shared_cell_types = shared_cell_types
  )

  save_csv(
    cd86_difference,
    file.path(
      table_dir,
      "CD86_source_target_difference.csv"
    )
  )

  plot_difference_heatmap(
    difference_df = cd86_difference,
    pathway = "CD86",
    output_pdf = file.path(
      figure_dir,
      "fig3d_pre_bcg_CD86_difference_heatmap.pdf"
    ),
    output_png = file.path(
      figure_dir,
      "fig3d_pre_bcg_CD86_difference_heatmap.png"
    )
  )

} else {

  warning(
    "CD86 pathway was not detected in either CellChat object."
  )

}


# ============================================================
# 30. OPTIONAL MHC-II DIFFERENTIAL HEATMAP
# ============================================================

if (
  "MHC-II" %in%
    available_pathways
) {

  mhcii_difference <- make_pathway_difference_df(
    df_resp = df_responder,
    df_nonresp = df_nonresponder,
    pathway = "MHC-II",
    shared_cell_types = shared_cell_types
  )

  save_csv(
    mhcii_difference,
    file.path(
      table_dir,
      "MHC-II_source_target_difference.csv"
    )
  )

  plot_difference_heatmap(
    difference_df = mhcii_difference,
    pathway = "MHC-II",
    output_pdf = file.path(
      figure_dir,
      "pre_bcg_MHC-II_difference_heatmap.pdf"
    ),
    output_png = file.path(
      figure_dir,
      "pre_bcg_MHC-II_difference_heatmap.png"
    )
  )

}


# ============================================================
# 31. MYELOID -> LYMPHOID CD86 SUMMARY
#
# Useful explicit audit of manuscript statement that CD86
# signalling is concentrated in monocyte-to-CD4/Treg axes.
# ============================================================

myeloid_cell_types <- intersect(
  c(
    "Classical Monocyte",
    "Activated Monocyte",
    "Non-classical Monocyte",
    "Dendritic Cell",
    "Plasmacytoid Dendritic Cell"
  ),
  shared_cell_types
)

lymphoid_targets <- intersect(
  c(
    "T cell (CD4+)",
    "T-reg",
    "T cell (CD8+)",
    "NK Cell",
    "NKT cell",
    "B-cell",
    "B-cell / Plasma Cell"
  ),
  shared_cell_types
)

cd86_myeloid_lymphoid <- bind_rows(
  df_responder %>%
    filter(
      pathway_name == "CD86",
      source %in% myeloid_cell_types,
      target %in% lymphoid_targets
    ) %>%
    mutate(
      outcome_group = "Responder"
    ),
  df_nonresponder %>%
    filter(
      pathway_name == "CD86",
      source %in% myeloid_cell_types,
      target %in% lymphoid_targets
    ) %>%
    mutate(
      outcome_group = "Non-responder"
    )
)

save_csv(
  cd86_myeloid_lymphoid,
  file.path(
    table_dir,
    "CD86_myeloid_to_lymphoid_interactions.csv"
  )
)


# ============================================================
# 32. PARAMETERS
# ============================================================

parameter_table <- data.frame(
  parameter = c(
    "input_rds",
    "group_by",
    "responder_label",
    "nonresponder_label",
    "database",
    "database_categories",
    "population_size",
    "min_cells",
    "future_workers",
    "RNA_preparation"
  ),
  value = c(
    input_rds,
    GROUP_COLUMN,
    RESPONDER_LABEL,
    NONRESPONDER_LABEL,
    "CellChatDB.human",
    paste(
      CELLCHAT_DATABASE_CATEGORIES,
      collapse = "; "
    ),
    POPULATION_SIZE,
    MIN_CELLS,
    N_WORKERS,
    "JoinLayers(RNA) + NormalizeData(RNA)"
  ),
  stringsAsFactors = FALSE
)

save_csv(
  parameter_table,
  file.path(
    outdir,
    "06_cellchat_pre_bcg_parameters.csv"
  )
)


# ============================================================
# 33. SESSION INFO
# ============================================================

sink(
  file.path(
    outdir,
    "06_cellchat_pre_bcg_sessionInfo.txt"
  )
)

print(
  sessionInfo()
)

sink()


# ============================================================
# 34. FINAL REPORT
# ============================================================

cat("\n")
cat("============================================================\n")
cat("PRE-BCG CELLCHAT ANALYSIS COMPLETE\n")
cat("============================================================\n")

cat(
  "\nResponders: ",
  format(
    ncol(
      responders
    ),
    big.mark = ","
  ),
  " cells\n",
  sep = ""
)

cat(
  "Non-responders: ",
  format(
    ncol(
      nonresponders
    ),
    big.mark = ","
  ),
  " cells\n",
  sep = ""
)

cat(
  "\nResponder inferred interactions: ",
  network_summary$n_inferred_interactions[
    network_summary$group ==
      "Responder"
  ],
  "\n",
  sep = ""
)

cat(
  "Non-responder inferred interactions: ",
  network_summary$n_inferred_interactions[
    network_summary$group ==
      "Non-responder"
  ],
  "\n",
  sep = ""
)

cat(
  "\nResponder total interaction strength: ",
  network_summary$total_interaction_strength[
    network_summary$group ==
      "Responder"
  ],
  "\n",
  sep = ""
)

cat(
  "Non-responder total interaction strength: ",
  network_summary$total_interaction_strength[
    network_summary$group ==
      "Non-responder"
  ],
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
