#!/usr/bin/env Rscript

# ============================================================
# 07_milor_pre_activated_monocytes.R
#
# Activated Monocyte pre-BCG MiloR workflow reconstructed from
# the original TFRI P3 3.Rmd and TFRI P3 4.Rmd notebooks.
#
# Key legacy features preserved:
#   - shared pre/post ATAC peak space
#   - 0.5% global shared-peak detection filter
#   - Activated Monocyte embedding built from pre + post cells
#   - RNA PCA 1:30 + ATAC LSI 2:30, column-wise z-scored
#   - Milo graph built only after subsetting the SCE to pre-BCG
#   - k = 30, prop = 0.05, refined = TRUE
#   - NonResponder reference; responseResponder contrast
#   - full nhoods(milo) membership used to map neighbourhoods
#   - five Hallmark module scores averaged by NR-up neighbourhood
#   - k-means (k = 3, seed = 1) in Hallmark-score space
#   - legacy cluster 2 ("pop2") used for chromVAR / ATAC follow-up
#
# Default local inputs:
#   G:/merged_obj.rds
#   G:/post_BCG_merged_obj.rds
#
# Paths may be overridden with:
#   PRE_BCG_RDS
#   POST_BCG_RDS
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(SingleCellExperiment)
  library(miloR)
  library(Matrix)
  library(S4Vectors)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(ggrepel)
  library(msigdbr)
  library(fgsea)
  library(GenomicRanges)
  library(GenomeInfoDb)
  library(JASPAR2020)
  library(TFBSTools)
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(BiocParallel)
  library(uwot)
  library(EnsDb.Hsapiens.v86)
})

# ============================================================
# 1. PATHS
# ============================================================

pre_rds <- Sys.getenv("PRE_BCG_RDS", unset = "G:/merged_obj.rds")
post_rds <- Sys.getenv("POST_BCG_RDS", unset = "G:/post_BCG_merged_obj.rds")
fragment_root <- Sys.getenv("FRAGMENT_ROOT", unset = "G:/TFRI P3 ATAC files")

outdir <- file.path("results", "07_milor_pre_activated_monocytes")
milo_dir <- file.path(outdir, "milo")
pathway_dir <- file.path(outdir, "pathway")
motif_dir <- file.path(outdir, "motif")
atac_dir <- file.path(outdir, "ATAC")
figure_dir <- file.path(outdir, "figures")
intermediate_dir <- file.path(outdir, "intermediate")

for (d in c(outdir, milo_dir, pathway_dir, motif_dir, atac_dir, figure_dir, intermediate_dir)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

if (!file.exists(pre_rds)) stop("Pre-BCG object not found: ", pre_rds)
if (!file.exists(post_rds)) stop("Post-BCG object not found: ", post_rds)
if (!dir.exists(fragment_root)) stop("Fragment root not found: ", fragment_root)

# ============================================================
# 2. PARAMETERS
# ============================================================

CELL_TYPE <- "Activated Monocyte"
CELL_TYPE_COL <- "cell_type"

SEED <- 1

SHARED_TOP_PER_TIMEPOINT <- 80000
GLOBAL_MIN_DETECTION_FRAC <- 0.005

RNA_NFEATURES <- 5000
PCA_NPCS <- 40
PCA_DIMS_USE <- 1:30

ATAC_TOP <- 20000
LSI_NPCS <- 40
LSI_DIMS_USE <- 2:30
MIN_PEAKS <- 2000

K_PARAM <- 30
PROP <- 0.05
REFINED <- TRUE
REDUCED_DIM <- "joint"

MIN_CELLS_TOTAL <- 1500
MIN_CELLS_PER_SAMPLE <- 25
MIN_SAMPLES_PER_GROUP <- 3
SPATIAL_FDR_CUTOFF <- 0.05

PATHWAY_NAMES <- c(
  "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  "HALLMARK_INTERFERON_ALPHA_RESPONSE",
  "HALLMARK_INFLAMMATORY_RESPONSE",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "HALLMARK_COMPLEMENT"
)

KMEANS_CENTERS <- 3
LEGACY_TF_CLUSTER <- 2

MOTIF_MIN_PCT <- 0.05
MOTIF_LOGFC_THRESHOLD <- 0
MOTIF_PADJ_CUTOFF <- 0.05
PEAK_PADJ_CUTOFF <- 0.05

AP1_CHROMVAR_PATTERN <- "JUN|FOS|FOSL|ATF|BATF|MAF"
AP1_PEAK_PATTERN <- "JUN|FOS|FOSL|BATF|ATF3|AP-1|AP1"

GSEA_NPERM <- 10000
GSEA_MIN_SIZE <- 10
GSEA_MAX_SIZE <- 500

# ============================================================
# 3. HELPERS
# ============================================================

save_csv <- function(x, path, row_names = FALSE) {
  write.csv(x, path, row.names = row_names)
}

save_fgsea_csv <- function(x, path) {
  out <- as.data.frame(x)

  if ("leadingEdge" %in% colnames(out)) {
    out$leadingEdge <- vapply(
      out$leadingEdge,
      paste,
      collapse = ";",
      FUN.VALUE = character(1)
    )
  }

  list_cols <- vapply(out, is.list, logical(1))
  if (any(list_cols)) {
    out[list_cols] <- lapply(
      out[list_cols],
      function(z) vapply(z, paste, collapse = ";", FUN.VALUE = character(1))
    )
  }

  save_csv(out, path)
}

normalize_response <- function(x) {
  x <- trimws(as.character(x))

  x[x %in% c(
    "Non-responder", "NonResponder", "Non_responder", "Non responder",
    "Recurrence", "early_recurrence", "Early recurrence"
  )] <- "NonResponder"

  x[x %in% c(
    "Responder", "Responder ", "Remission",
    "recurrence_free", "Recurrence free"
  )] <- "Responder"

  x
}

extract_response <- function(object) {
  md <- object@meta.data

  for (nm in c("response", "Outcome", "analysis_outcome")) {
    if (nm %in% colnames(md)) {
      out <- normalize_response(md[[nm]])
      if (all(out %in% c("Responder", "NonResponder"))) return(out)
    }
  }

  stop("Could not identify a clean response/outcome column.")
}

extract_sample_base <- function(object) {
  md <- object@meta.data

  for (nm in c("patient_id", "sample_id", "orig.ident", "sample")) {
    if (nm %in% colnames(md)) {
      out <- trimws(as.character(md[[nm]]))
      if (!anyNA(out) && all(nzchar(out))) return(out)
    }
  }

  stop("Could not identify patient_id, sample_id, orig.ident, or sample.")
}

relink_fragment_paths <- function(object, fragment_root, assay = "ATAC") {
  frags <- Fragments(object[[assay]])

  if (length(frags) == 0) {
    stop("No Fragment objects found in assay ", assay, ".")
  }

  new_paths <- character(length(frags))

  for (i in seq_along(frags)) {
    old_path <- GetFragmentData(frags[[i]], slot = "path")
    old_path_slash <- chartr("\\", "/", old_path)

    parent <- basename(dirname(old_path_slash))
    sample_name <- if (identical(parent, "outs")) {
      basename(dirname(dirname(old_path_slash)))
    } else {
      parent
    }

    new_path <- file.path(fragment_root, sample_name, "atac_fragments.tsv.gz")
    new_index <- paste0(new_path, ".tbi")

    if (!file.exists(new_path)) {
      stop(
        "Relinked fragment file not found for ", sample_name, ": ", new_path
      )
    }

    if (!file.exists(new_index)) {
      stop(
        "Tabix index not found for ", sample_name, ": ", new_index
      )
    }

    message("  Relinking ", sample_name, " -> ", new_path)

    # UpdatePath preserves the existing Seurat-cell -> fragment-barcode
    # mapping and verifies that the relocated fragment file/index match
    # the files used to create the original Fragment object.
    frags[[i]] <- UpdatePath(
      frags[[i]],
      new.path = new_path,
      verbose = FALSE
    )

    new_paths[i] <- new_path
  }

  Fragments(object[[assay]]) <- NULL
  Fragments(object[[assay]]) <- frags

  invisible(list(object = object, paths = new_paths))
}

collapse_counts_layers <- function(object, assay = "RNA") {
  layer_names <- Layers(object[[assay]])
  count_layers <- layer_names[grepl("^counts", layer_names)]

  if (length(count_layers) == 0) {
    stop("No count layers found in assay ", assay)
  }

  mats <- lapply(
    count_layers,
    function(layer_name) GetAssayData(object, assay = assay, layer = layer_name)
  )

  common_features <- Reduce(intersect, lapply(mats, rownames))
  mats <- lapply(mats, function(m) m[common_features, , drop = FALSE])

  out <- if (length(mats) == 1) mats[[1]] else Reduce(Matrix::cbind2, mats)

  missing_cells <- setdiff(colnames(object), colnames(out))
  if (length(missing_cells) > 0) {
    stop("RNA count layers are missing ", length(missing_cells), " cells.")
  }

  out[, colnames(object), drop = FALSE]
}

get_top_atac_peaks <- function(object, n_top = 80000) {
  if (!"ATAC" %in% names(object@assays)) {
    stop("ATAC assay is missing.")
  }

  DefaultAssay(object) <- "ATAC"
  object <- RunTFIDF(object)
  object <- FindTopFeatures(object, min.cutoff = "q0")

  head(VariableFeatures(object), n_top)
}

build_shared_counts_and_detection <- function(
  object,
  peaks_gr,
  activated_cells,
  sample_col = "orig.ident"
) {
  md <- object@meta.data

  if (!sample_col %in% colnames(md)) {
    stop("Required sample column is missing: ", sample_col)
  }

  DefaultAssay(object) <- "ATAC"
  fragments <- Fragments(object)

  samples <- unique(as.character(md[[sample_col]]))
  detected <- integer(length(peaks_gr))
  names(detected) <- NULL
  activated_parts <- list()

  for (s in samples) {
    message("  FeatureMatrix: ", s)

    cells_s <- rownames(md)[as.character(md[[sample_col]]) == s]

    m_s <- FeatureMatrix(
      fragments = fragments,
      features = peaks_gr,
      cells = cells_s
    )

    detected <- detected + Matrix::rowSums(m_s > 0)

    act_s <- intersect(activated_cells, colnames(m_s))
    if (length(act_s) > 0) {
      activated_parts[[length(activated_parts) + 1L]] <- m_s[, act_s, drop = FALSE]
    }

    rm(m_s)
    gc(verbose = FALSE)
  }

  if (length(activated_parts) == 0) {
    stop("No Activated Monocyte ATAC counts were generated.")
  }

  activated_counts <- if (length(activated_parts) == 1) {
    activated_parts[[1]]
  } else {
    Reduce(Matrix::cbind2, activated_parts)
  }

  activated_counts <- activated_counts[, activated_cells, drop = FALSE]

  list(
    counts = activated_counts,
    detected = detected,
    n_cells = ncol(object)
  )
}

get_motif_names <- function(object) {
  x <- GetMotifData(object = object, assay = "ATAC", slot = "motif.names")

  if (is.data.frame(x)) {
    if ("motif.name" %in% colnames(x)) {
      out <- x$motif.name
      names(out) <- rownames(x)
      return(out)
    }

    if (ncol(x) == 1) {
      out <- x[[1]]
      names(out) <- rownames(x)
      return(out)
    }
  }

  x
}

nhood_membership_long <- function(milo) {
  nh <- nhoods(milo)

  if (is.null(rownames(nh))) {
    rownames(nh) <- colnames(milo)
  }

  if (is.null(colnames(nh))) {
    colnames(nh) <- as.character(seq_len(ncol(nh)))
  }

  if (!inherits(nh, "dgCMatrix")) {
    nh <- as(nh, "dgCMatrix")
  }

  sm <- Matrix::summary(nh)

  data.frame(
    cell_prefixed = rownames(nh)[sm$i],
    Nhood = suppressWarnings(as.integer(colnames(nh)[sm$j])),
    stringsAsFactors = FALSE
  )
}

# ============================================================
# 4. LOAD PRE + POST OBJECTS
# ============================================================

cat("\n============================================================\n")
cat("Loading pre- and post-BCG objects\n")
cat("============================================================\n")

pre_bcg <- readRDS(pre_rds)
post_bcg <- readRDS(post_rds)

cat("\nRelinking fragment files to indexed local copies\n")
pre_relinked <- relink_fragment_paths(pre_bcg, fragment_root, assay = "ATAC")
pre_bcg <- pre_relinked$object

post_relinked <- relink_fragment_paths(post_bcg, fragment_root, assay = "ATAC")
post_bcg <- post_relinked$object

fragment_manifest <- data.frame(
  timepoint = c(rep("pre", length(pre_relinked$paths)), rep("post", length(post_relinked$paths))),
  fragment_path = c(pre_relinked$paths, post_relinked$paths),
  index_path = paste0(c(pre_relinked$paths, post_relinked$paths), ".tbi"),
  stringsAsFactors = FALSE
)

save_csv(
  fragment_manifest,
  file.path(intermediate_dir, "fragment_path_manifest.csv")
)

rm(pre_relinked, post_relinked)
gc(verbose = FALSE)

for (nm in c("RNA", "ATAC")) {
  if (!nm %in% names(pre_bcg@assays)) stop("Pre object missing assay: ", nm)
  if (!nm %in% names(post_bcg@assays)) stop("Post object missing assay: ", nm)
}

if (!CELL_TYPE_COL %in% colnames(pre_bcg@meta.data)) {
  stop("Pre object missing metadata column: ", CELL_TYPE_COL)
}
if (!CELL_TYPE_COL %in% colnames(post_bcg@meta.data)) {
  stop("Post object missing metadata column: ", CELL_TYPE_COL)
}

pre_act_cells <- colnames(pre_bcg)[pre_bcg@meta.data[[CELL_TYPE_COL]] == CELL_TYPE]
post_act_cells <- colnames(post_bcg)[post_bcg@meta.data[[CELL_TYPE_COL]] == CELL_TYPE]

cat("Pre Activated Monocytes: ", format(length(pre_act_cells), big.mark = ","), "\n", sep = "")
cat("Post Activated Monocytes: ", format(length(post_act_cells), big.mark = ","), "\n", sep = "")

if (length(pre_act_cells) < MIN_CELLS_TOTAL) {
  stop("Too few pre-BCG Activated Monocytes.")
}
if (length(post_act_cells) == 0) {
  stop("No post-BCG Activated Monocytes found.")
}

# Keep the original pre Activated Monocyte object for the legacy
# downstream RNA / ATAC analyses.
act <- subset(pre_bcg, cells = pre_act_cells)

# ============================================================
# 5. DEFINE THE LEGACY SHARED PRE/POST ATAC PEAK SPACE
# ============================================================

cat("\n============================================================\n")
cat("Building shared pre/post ATAC peak space\n")
cat("============================================================\n")

top_pre <- get_top_atac_peaks(pre_bcg, SHARED_TOP_PER_TIMEPOINT)
gc(verbose = FALSE)

top_post <- get_top_atac_peaks(post_bcg, SHARED_TOP_PER_TIMEPOINT)
gc(verbose = FALSE)

gr_pre <- StringToGRanges(top_pre, sep = c("-", "-"))
gr_post <- StringToGRanges(top_post, sep = c("-", "-"))

combined_peaks <- reduce(c(gr_pre, gr_post))
combined_peaks <- keepStandardChromosomes(combined_peaks, pruning.mode = "coarse")

peak_width <- width(combined_peaks)
combined_peaks <- combined_peaks[peak_width > 20 & peak_width < 10000]

cat("Shared candidate peaks: ", format(length(combined_peaks), big.mark = ","), "\n", sep = "")

saveRDS(
  combined_peaks,
  file.path(intermediate_dir, "combined_shared_peaks_before_detection_filter.rds")
)

# ============================================================
# 6. COUNT SHARED PEAKS
#
# Legacy workflow counted this shared peak set across all cells,
# then retained peaks detected in >=0.5% of all pre+post cells.
# To reduce peak-memory use, this implementation accumulates the
# all-cell detection counts sample-by-sample while retaining the
# full matrix only for Activated Monocytes.
# ============================================================

cat("\n============================================================\n")
cat("Counting shared peaks in pre-BCG cells\n")
cat("============================================================\n")

pre_shared <- build_shared_counts_and_detection(
  pre_bcg,
  combined_peaks,
  activated_cells = pre_act_cells,
  sample_col = "orig.ident"
)

cat("\n============================================================\n")
cat("Counting shared peaks in post-BCG cells\n")
cat("============================================================\n")

post_shared <- build_shared_counts_and_detection(
  post_bcg,
  combined_peaks,
  activated_cells = post_act_cells,
  sample_col = "orig.ident"
)

if (!identical(rownames(pre_shared$counts), rownames(post_shared$counts))) {
  stop("Pre/post shared ATAC feature orders do not match.")
}

detected_total <- pre_shared$detected + post_shared$detected
n_cells_total <- pre_shared$n_cells + post_shared$n_cells
min_cells_global <- ceiling(n_cells_total * GLOBAL_MIN_DETECTION_FRAC)

keep_peak_idx <- detected_total >= min_cells_global

pre_atac <- pre_shared$counts[keep_peak_idx, , drop = FALSE]
post_atac <- post_shared$counts[keep_peak_idx, , drop = FALSE]

cat(
  "Shared peaks after 0.5% global detection filter: ",
  format(nrow(pre_atac), big.mark = ","),
  "\n",
  sep = ""
)

save_csv(
  data.frame(
    peak = rownames(pre_shared$counts),
    detected_cells = detected_total,
    keep = keep_peak_idx
  ),
  file.path(intermediate_dir, "shared_peak_detection_filter.csv")
)

rm(pre_shared, post_shared, detected_total, keep_peak_idx)
gc(verbose = FALSE)

# ============================================================
# 7. BUILD THE PRE+POST ACTIVATED-MONOCYTE OBJECT
# ============================================================

cat("\n============================================================\n")
cat("Building combined pre/post Activated Monocyte object\n")
cat("============================================================\n")

pre_act_for_rna <- subset(pre_bcg, cells = pre_act_cells)
post_act_for_rna <- subset(post_bcg, cells = post_act_cells)

pre_rna <- collapse_counts_layers(pre_act_for_rna, assay = "RNA")
post_rna <- collapse_counts_layers(post_act_for_rna, assay = "RNA")

common_genes <- intersect(rownames(pre_rna), rownames(post_rna))
pre_rna <- pre_rna[common_genes, , drop = FALSE]
post_rna <- post_rna[common_genes, , drop = FALSE]

meta_pre <- pre_act_for_rna@meta.data[colnames(pre_rna), , drop = FALSE]
meta_post <- post_act_for_rna@meta.data[colnames(post_rna), , drop = FALSE]

meta_pre$response <- extract_response(pre_act_for_rna)
meta_post$response <- extract_response(post_act_for_rna)

pre_sample_base <- extract_sample_base(pre_act_for_rna)
post_sample_base <- extract_sample_base(post_act_for_rna)

meta_pre$timepoint <- "pre"
meta_post$timepoint <- "post"
meta_pre$sample <- paste0(pre_sample_base, "_pre")
meta_post$sample <- paste0(post_sample_base, "_post")

pre_prefixed <- paste0("pre_", colnames(pre_rna))
post_prefixed <- paste0("post_", colnames(post_rna))

colnames(pre_rna) <- pre_prefixed
colnames(post_rna) <- post_prefixed
colnames(pre_atac) <- pre_prefixed
colnames(post_atac) <- post_prefixed

rownames(meta_pre) <- pre_prefixed
rownames(meta_post) <- post_prefixed

rna_all <- Matrix::cbind2(pre_rna, post_rna)
atac_all <- Matrix::cbind2(pre_atac, post_atac)
meta_all <- rbind(meta_pre, meta_post)

stopifnot(identical(colnames(rna_all), colnames(atac_all)))
meta_all <- meta_all[colnames(rna_all), , drop = FALSE]

obj <- CreateSeuratObject(
  counts = rna_all,
  assay = "RNA",
  meta.data = meta_all
)

obj[["ATAC_shared"]] <- CreateChromatinAssay(counts = atac_all)

rm(
  post_bcg, pre_act_for_rna, post_act_for_rna,
  pre_rna, post_rna, pre_atac, post_atac,
  rna_all, atac_all, meta_pre, meta_post, meta_all
)
gc(verbose = FALSE)

# ============================================================
# 8. LEGACY ACTIVATED-MONOCYTE JOINT EMBEDDING
#
# The original Milo workflow re-embedded each cell type using
# BOTH pre- and post-treatment cells, then subset the SCE to
# pre-BCG before constructing the Milo graph.
# ============================================================

cat("\n============================================================\n")
cat("Building pre+post Activated Monocyte joint embedding\n")
cat("============================================================\n")

set.seed(SEED)

DefaultAssay(obj) <- "RNA"
obj <- NormalizeData(obj, verbose = FALSE)
obj <- FindVariableFeatures(obj, nfeatures = RNA_NFEATURES, verbose = FALSE)
obj <- ScaleData(obj, features = VariableFeatures(obj), verbose = FALSE)
obj <- RunPCA(
  obj,
  npcs = PCA_NPCS,
  features = VariableFeatures(obj),
  verbose = FALSE
)

DefaultAssay(obj) <- "ATAC_shared"
obj <- RunTFIDF(obj)
obj <- FindTopFeatures(obj, min.cutoff = "q0")

top_peaks <- VariableFeatures(obj)
if (length(top_peaks) > ATAC_TOP) {
  top_peaks <- top_peaks[seq_len(ATAC_TOP)]
}

if (length(top_peaks) < MIN_PEAKS) {
  stop("Too few ATAC_shared peaks after feature selection: ", length(top_peaks))
}

VariableFeatures(obj) <- top_peaks

obj <- RunSVD(
  obj,
  n = LSI_NPCS,
  reduction.name = "lsi",
  reduction.key = "LSI_",
  features = top_peaks
)

pca_mat <- Embeddings(obj, "pca")[, PCA_DIMS_USE, drop = FALSE]
lsi_mat <- Embeddings(obj, "lsi")[, LSI_DIMS_USE, drop = FALSE]

pca_z <- scale(pca_mat)
lsi_z <- scale(lsi_mat)
joint <- cbind(pca_z, lsi_z)

cat(
  "Joint embedding: ",
  nrow(joint), " cells x ", ncol(joint), " dimensions\n",
  sep = ""
)

sce <- as.SingleCellExperiment(obj)
reducedDim(sce, "PCA") <- pca_mat
reducedDim(sce, "LSI") <- lsi_mat
reducedDim(sce, REDUCED_DIM) <- joint

md <- obj@meta.data[colnames(obj), , drop = FALSE]
colData(sce) <- S4Vectors::DataFrame(md)

saveRDS(
  sce,
  file.path(intermediate_dir, "sce_ActivatedMonocyte_pre_post_joint.rds")
)

# ============================================================
# 9. SUBSET TO PRE-BCG AND FILTER LOW-CELL SAMPLES
# ============================================================

md_all <- as.data.frame(colData(sce))
md_all$response <- normalize_response(md_all$response)

keep_pre <- md_all$timepoint == "pre"
sce_pre <- sce[, keep_pre]
md_pre <- as.data.frame(colData(sce_pre))
md_pre$response <- normalize_response(md_pre$response)

sample_tab <- table(md_pre$sample, md_pre$response)
per_sample_total <- rowSums(sample_tab)
keep_samples <- names(per_sample_total)[per_sample_total >= MIN_CELLS_PER_SAMPLE]

sce_pre <- sce_pre[, md_pre$sample %in% keep_samples]
md_pre <- as.data.frame(colData(sce_pre))
md_pre$response <- normalize_response(md_pre$response)

n_responder_samples <- length(unique(md_pre$sample[md_pre$response == "Responder"]))
n_nonresponder_samples <- length(unique(md_pre$sample[md_pre$response == "NonResponder"]))

if (ncol(sce_pre) < MIN_CELLS_TOTAL) {
  stop("Too few pre-BCG Activated Monocytes after filtering.")
}

if (
  n_responder_samples < MIN_SAMPLES_PER_GROUP ||
  n_nonresponder_samples < MIN_SAMPLES_PER_GROUP
) {
  stop(
    "Too few samples per group after filtering: Responder=",
    n_responder_samples,
    "; NonResponder=",
    n_nonresponder_samples
  )
}

cat("Pre-BCG cells entering Milo: ", format(ncol(sce_pre), big.mark = ","), "\n", sep = "")
cat("Responder samples: ", n_responder_samples, "\n", sep = "")
cat("NonResponder samples: ", n_nonresponder_samples, "\n", sep = "")

# ============================================================
# 10. MILOR DIFFERENTIAL ABUNDANCE
# ============================================================

cat("\n============================================================\n")
cat("Running MiloR\n")
cat("============================================================\n")

set.seed(SEED)
milo <- Milo(sce_pre)

d_use <- ncol(reducedDim(milo, REDUCED_DIM))

milo <- buildGraph(
  milo,
  k = K_PARAM,
  d = d_use,
  reduced.dim = REDUCED_DIM
)

set.seed(SEED)
milo <- makeNhoods(
  milo,
  prop = PROP,
  k = K_PARAM,
  d = d_use,
  refined = REFINED,
  reduced_dims = REDUCED_DIM
)

md_milo <- as.data.frame(colData(milo))
md_milo$response <- normalize_response(md_milo$response)
md_milo <- md_milo[colnames(milo), , drop = FALSE]

milo <- countCells(
  milo,
  meta.data = md_milo,
  sample = "sample"
)

Y <- nhoodCounts(milo)
samples <- colnames(Y)

sample_design <- unique(md_milo[, c("sample", "response"), drop = FALSE])
sample_design <- sample_design[!duplicated(sample_design$sample), , drop = FALSE]
sample_design <- sample_design[match(samples, sample_design$sample), , drop = FALSE]

if (anyNA(sample_design$sample) || anyNA(sample_design$response)) {
  stop("Milo sample design could not be aligned to nhoodCounts.")
}

design_df <- data.frame(
  response = factor(
    normalize_response(sample_design$response),
    levels = c("NonResponder", "Responder")
  ),
  row.names = sample_design$sample
)

da <- testNhoods(
  milo,
  design = ~ response,
  design.df = design_df,
  model.contrasts = "responseResponder",
  reduced.dim = REDUCED_DIM
)

da <- as.data.frame(da)

if (!"Nhood" %in% colnames(da)) {
  da$Nhood <- rownames(da)
}

da$Nhood <- suppressWarnings(as.integer(as.character(da$Nhood)))

if (anyNA(da$Nhood)) {
  da$Nhood <- seq_len(nrow(da))
}

if (!"SpatialFDR" %in% colnames(da)) {
  if ("FDR" %in% colnames(da)) {
    da$SpatialFDR <- da$FDR
  } else {
    stop("Milo output contains neither SpatialFDR nor FDR.")
  }
}

da$direction <- ifelse(
  da$logFC < 0,
  "NonResponder_enriched",
  ifelse(da$logFC > 0, "Responder_enriched", "No_direction")
)

da$significant <- !is.na(da$SpatialFDR) & da$SpatialFDR < SPATIAL_FDR_CUTOFF

save_csv(da, file.path(milo_dir, "ActivatedMonocyte_pre_Milo_DA.csv"))
saveRDS(da, file.path(milo_dir, "ActivatedMonocyte_pre_Milo_DA.rds"))
saveRDS(milo, file.path(milo_dir, "ActivatedMonocyte_pre_Milo_object.rds"))

n_nr_nhoods <- sum(
  da$significant & da$logFC < 0,
  na.rm = TRUE
)
n_r_nhoods <- sum(
  da$significant & da$logFC > 0,
  na.rm = TRUE
)

cat("Milo neighbourhoods: ", nrow(da), "\n", sep = "")
cat("Significant NonResponder-enriched neighbourhoods: ", n_nr_nhoods, "\n", sep = "")
cat("Significant Responder-enriched neighbourhoods: ", n_r_nhoods, "\n", sep = "")

# ============================================================
# 11. FULL NEIGHBOURHOOD MEMBERSHIP
#
# IMPORTANT: nhoods(milo), not nhoodIndex(milo).
# ============================================================

nh_map <- nhood_membership_long(milo)

if (anyNA(nh_map$Nhood)) {
  stop("Neighbourhood IDs could not be converted to integers.")
}

nh_map$cell <- sub("^pre_", "", nh_map$cell_prefixed)

nh_map_act <- nh_map %>%
  filter(cell %in% colnames(act)) %>%
  select(Nhood, cell, cell_prefixed)

save_csv(
  nh_map_act,
  file.path(milo_dir, "ActivatedMonocyte_pre_neighbourhood_cells.csv")
)

nr_up_nhoods <- da %>%
  filter(
    SpatialFDR < SPATIAL_FDR_CUTOFF,
    logFC < 0
  ) %>%
  pull(Nhood) %>%
  unique()

nr_up_cells <- nh_map_act %>%
  filter(Nhood %in% nr_up_nhoods) %>%
  pull(cell) %>%
  unique()

cat(
  "Unique cells in any NR-up neighbourhood: ",
  format(length(nr_up_cells), big.mark = ","),
  "\n",
  sep = ""
)

# ============================================================
# 12. MILO OVERLAY
# ============================================================

try({
  set.seed(SEED)

  um <- uwot::umap(
    reducedDim(milo, REDUCED_DIM),
    n_neighbors = 30,
    min_dist = 0.3,
    metric = "cosine",
    verbose = FALSE
  )

  rownames(um) <- colnames(milo)
  reducedDim(milo, "UMAP") <- um
  milo <- buildNhoodGraph(milo)

  p_milo <- plotNhoodGraphDA(
    milo,
    milo_res = da,
    layout = "UMAP"
  ) +
    ggtitle("Activated Monocyte | pre-BCG | Responder vs NonResponder")

  ggsave(
    file.path(figure_dir, "fig4a_ActivatedMonocyte_pre_Milo_DA.pdf"),
    p_milo,
    width = 7,
    height = 6
  )
}, silent = FALSE)

# ============================================================
# 13. FIVE HALLMARK MODULE SCORES
# ============================================================

cat("\n============================================================\n")
cat("Scoring NR-up neighbourhoods by Hallmark activity\n")
cat("============================================================\n")

hallmark <- msigdbr(
  species = "Homo sapiens",
  collection = "H"
)

hallmark_list <- split(hallmark$gene_symbol, hallmark$gs_name)
hallmark_list <- lapply(hallmark_list, unique)

missing_pathways <- setdiff(PATHWAY_NAMES, names(hallmark_list))
if (length(missing_pathways) > 0) {
  stop("Missing Hallmark pathways: ", paste(missing_pathways, collapse = ", "))
}

pathways_use <- hallmark_list[PATHWAY_NAMES]

DefaultAssay(act) <- "RNA"
act <- JoinLayers(act, assay = "RNA")
act <- NormalizeData(act, verbose = FALSE)
act <- ScaleData(act, verbose = FALSE)

set.seed(SEED)
act <- AddModuleScore(
  act,
  features = pathways_use,
  name = "HALLMARK_",
  seed = SEED
)

score_cols_chr <- paste0("HALLMARK_", seq_along(pathways_use))

meta_scores <- act@meta.data %>%
  rownames_to_column("cell") %>%
  select(cell, all_of(score_cols_chr))

meta_nh_scores <- nh_map_act %>%
  inner_join(meta_scores, by = "cell") %>%
  filter(Nhood %in% nr_up_nhoods)

nhood_scores_NR <- meta_nh_scores %>%
  group_by(Nhood) %>%
  summarise(
    nCells = n(),
    across(all_of(score_cols_chr), ~ mean(.x, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  left_join(da, by = "Nhood")

rename_map <- setNames(score_cols_chr, PATHWAY_NAMES)
nhood_scores_NR <- nhood_scores_NR %>% rename(!!!rename_map)

save_csv(
  nhood_scores_NR,
  file.path(pathway_dir, "ActivatedMonocyte_NRup_neighbourhood_Hallmark_scores.csv")
)

# ============================================================
# 14. LEGACY K-MEANS SUBSTRUCTURE OF NR-UP NEIGHBOURHOODS
#
# The legacy notebook calls kmeans(scale(X), centers = 3) and
# subsequently uses cluster 2 for TF analysis. X was not defined
# in the saved notebook; here it is explicitly reconstructed as
# the five-neighbourhood Hallmark score matrix used immediately
# around that code block.
# ============================================================

plot_df <- nhood_scores_NR

X <- as.matrix(plot_df[, PATHWAY_NAMES, drop = FALSE])

if (any(!is.finite(X))) {
  stop("Non-finite Hallmark values detected before k-means.")
}

set.seed(SEED)
km <- kmeans(scale(X), centers = KMEANS_CENTERS)
plot_df$cluster <- factor(km$cluster)

cluster_summary <- plot_df %>%
  group_by(cluster) %>%
  summarise(
    n_neighbourhoods = n(),
    across(all_of(PATHWAY_NAMES), ~ mean(.x, na.rm = TRUE)),
    .groups = "drop"
  )

save_csv(
  plot_df,
  file.path(pathway_dir, "ActivatedMonocyte_NRup_neighbourhood_clusters.csv")
)

save_csv(
  cluster_summary,
  file.path(pathway_dir, "ActivatedMonocyte_NRup_cluster_Hallmark_means.csv")
)

inflam_cluster <- cluster_summary %>%
  arrange(desc(HALLMARK_TNFA_SIGNALING_VIA_NFKB)) %>%
  slice(1) %>%
  pull(cluster)

cat("Highest-TNF Hallmark cluster: ", as.character(inflam_cluster), "\n", sep = "")
cat("Legacy TF-analysis cluster: ", LEGACY_TF_CLUSTER, "\n", sep = "")

pop2_nhoods <- plot_df %>%
  filter(as.integer(as.character(cluster)) == LEGACY_TF_CLUSTER) %>%
  pull(Nhood) %>%
  unique()

pop2_cells <- nh_map_act %>%
  filter(Nhood %in% pop2_nhoods) %>%
  pull(cell) %>%
  unique()

act$pop2 <- ifelse(colnames(act) %in% pop2_cells, "pop2", "other")

save_csv(
  data.frame(
    cell = colnames(act),
    pop2 = act$pop2,
    stringsAsFactors = FALSE
  ),
  file.path(pathway_dir, "ActivatedMonocyte_pop2_cell_labels.csv")
)

save_csv(
  data.frame(Nhood = pop2_nhoods),
  file.path(pathway_dir, "ActivatedMonocyte_pop2_neighbourhoods.csv")
)

cat("Cluster-2 neighbourhoods: ", length(pop2_nhoods), "\n", sep = "")
cat("Unique cluster-2 cells: ", format(length(pop2_cells), big.mark = ","), "\n", sep = "")

# Historical pre-BCG UMAP centroids are used only for display.
if ("umap" %in% names(act@reductions)) {
  umap_df <- as.data.frame(Embeddings(act, "umap"))
  umap_df$cell <- rownames(umap_df)
  colnames(umap_df)[1:2] <- c("UMAP1", "UMAP2")

  nh_centroids <- nh_map_act %>%
    inner_join(umap_df[, c("cell", "UMAP1", "UMAP2")], by = "cell") %>%
    group_by(Nhood) %>%
    summarise(
      UMAP1 = median(UMAP1, na.rm = TRUE),
      UMAP2 = median(UMAP2, na.rm = TRUE),
      .groups = "drop"
    )

  plot_df_umap <- plot_df %>%
    left_join(nh_centroids, by = "Nhood")

  p_cluster <- ggplot(
    plot_df_umap,
    aes(x = UMAP1, y = UMAP2)
  ) +
    geom_point(aes(size = nCells, shape = cluster), alpha = 0.9) +
    theme_classic() +
    labs(
      title = "NR-up Activated Monocyte neighbourhood clusters",
      size = "Cells",
      shape = "Cluster"
    )

  ggsave(
    file.path(figure_dir, "fig4_NRup_Hallmark_clusters.pdf"),
    p_cluster,
    width = 7,
    height = 6
  )
}

cluster_long <- cluster_summary %>%
  pivot_longer(
    cols = all_of(PATHWAY_NAMES),
    names_to = "Pathway",
    values_to = "MeanScore"
  )

p_heat <- ggplot(
  cluster_long,
  aes(x = Pathway, y = cluster, fill = MeanScore)
) +
  geom_tile() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    title = "Mean Hallmark activity per NR-up neighbourhood cluster",
    x = NULL,
    y = "Cluster"
  )

ggsave(
  file.path(figure_dir, "fig4_NRup_cluster_Hallmark_heatmap.pdf"),
  p_heat,
  width = 9,
  height = 4
)

# ============================================================
# 15. PREPARE ORIGINAL PRE-BCG ATAC ASSAY FOR MOTIF ANALYSIS
# ============================================================

cat("\n============================================================\n")
cat("Preparing JASPAR2020 motifs and chromVAR for cluster 2\n")
cat("============================================================\n")

DefaultAssay(act) <- "ATAC"

atac <- act[["ATAC"]]
gr <- granges(atac)
names(gr) <- rownames(atac)

std_chrs <- c(paste0("chr", 1:22), "chrX", "chrY", "chrM")
chr_vec <- as.character(seqnames(gr))
keep_idx <- !is.na(chr_vec) & chr_vec %in% std_chrs

gr_keep <- gr[keep_idx]

counts_keep <- GetAssayData(
  act,
  assay = "ATAC",
  layer = "counts"
)[keep_idx, , drop = FALSE]

frags <- Fragments(atac)
annotation_original <- tryCatch(Annotation(atac), error = function(e) NULL)

atac_new <- CreateChromatinAssay(
  counts = counts_keep,
  ranges = gr_keep,
  fragments = frags,
  genome = "hg38"
)

if (!is.null(annotation_original) && length(annotation_original) > 0) {
  Annotation(atac_new) <- annotation_original
}

act[["ATAC"]] <- atac_new
DefaultAssay(act) <- "ATAC"

register(SerialParam())

# Rebuilt assays do not reliably retain GC-bias metadata.
act <- RegionStats(
  object = act,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  assay = "ATAC",
  verbose = FALSE
)

pfm <- getMatrixSet(
  JASPAR2020,
  opts = list(
    species = 9606,
    all_versions = FALSE
  )
)

act <- AddMotifs(
  act,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  pfm = pfm
)

act <- RunChromVAR(
  act,
  genome = BSgenome.Hsapiens.UCSC.hg38
)

# ============================================================
# 16. DIFFERENTIAL CHROMVAR ACTIVITY: POP2 VS OTHER
# ============================================================

DefaultAssay(act) <- "chromvar"
Idents(act) <- "pop2"

motif_da <- FindMarkers(
  act,
  ident.1 = "pop2",
  ident.2 = "other",
  test.use = "wilcox",
  min.pct = MOTIF_MIN_PCT,
  logfc.threshold = MOTIF_LOGFC_THRESHOLD
)

motif_names <- get_motif_names(act)

motif_da$motif_id <- rownames(motif_da)
motif_da$tf <- vapply(
  motif_da$motif_id,
  function(id) {
    if (!is.null(names(motif_names)) && id %in% names(motif_names)) {
      paste(unique(as.character(motif_names[id])), collapse = ";")
    } else {
      NA_character_
    }
  },
  FUN.VALUE = character(1)
)

motif_da$AP1 <- grepl(
  AP1_CHROMVAR_PATTERN,
  motif_da$tf,
  ignore.case = TRUE
)

motif_da <- motif_da %>% arrange(p_val_adj)

save_csv(
  motif_da,
  file.path(motif_dir, "ActivatedMonocyte_pop2_vs_other_chromVAR.csv")
)

ap1_pop2_up <- motif_da %>%
  filter(
    AP1,
    avg_log2FC > 0,
    p_val_adj < MOTIF_PADJ_CUTOFF
  ) %>%
  arrange(p_val_adj)

save_csv(
  ap1_pop2_up,
  file.path(motif_dir, "ActivatedMonocyte_pop2_significant_AP1_chromVAR.csv")
)

motif_plot <- motif_da %>%
  mutate(
    neglog10padj = -log10(pmax(p_val_adj, .Machine$double.xmin)),
    label = ifelse(AP1, tf, NA_character_)
  )

p_motif <- ggplot(
  motif_plot,
  aes(x = avg_log2FC, y = neglog10padj)
) +
  geom_point(alpha = 0.35, size = 1.2) +
  geom_point(
    data = motif_plot %>% filter(AP1),
    size = 2
  ) +
  geom_text_repel(
    data = motif_plot %>% filter(AP1),
    aes(label = tf),
    size = 3,
    max.overlaps = Inf
  ) +
  geom_vline(xintercept = 0, linetype = 2) +
  theme_classic(base_size = 12) +
  labs(
    x = "chromVAR deviation (pop2 - other)",
    y = "-log10(adjusted p-value)",
    title = "Activated Monocyte TF motif activity"
  )

ggsave(
  file.path(figure_dir, "fig4_ActivatedMonocyte_pop2_chromVAR_AP1.pdf"),
  p_motif,
  width = 7,
  height = 6
)

# ============================================================
# 17. PEAK-LEVEL DIFFERENTIAL ACCESSIBILITY: POP2 VS OTHER
# ============================================================

cat("\n============================================================\n")
cat("Peak-level differential accessibility: pop2 vs other\n")
cat("============================================================\n")

DefaultAssay(act) <- "ATAC"
Idents(act) <- "pop2"

if (!"nCount_ATAC" %in% colnames(act@meta.data)) {
  stop("nCount_ATAC is required for logistic-regression peak DA.")
}

da_peaks <- FindMarkers(
  act,
  ident.1 = "pop2",
  ident.2 = "other",
  test.use = "LR",
  latent.vars = "nCount_ATAC",
  only.pos = TRUE
)

save_csv(
  da_peaks,
  file.path(atac_dir, "ActivatedMonocyte_pop2_vs_other_DA_peaks.csv"),
  row_names = TRUE
)

# ============================================================
# 18. NEAREST-GENE ANNOTATION
# ============================================================

peak_gr <- StringToGRanges(rownames(da_peaks), sep = c("-", "-"))
names(peak_gr) <- rownames(da_peaks)

annotations <- suppressWarnings(
  GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
)

seqlevelsStyle(peak_gr) <- "UCSC"
seqlevelsStyle(annotations) <- "UCSC"

common_seq <- intersect(seqlevels(peak_gr), seqlevels(annotations))
peak_gr2 <- keepSeqlevels(peak_gr, common_seq, pruning.mode = "coarse")
ann2 <- keepSeqlevels(annotations, common_seq, pruning.mode = "coarse")

hits <- distanceToNearest(peak_gr2, ann2, ignore.strand = TRUE)
q <- queryHits(hits)
s <- subjectHits(hits)

peak_to_gene <- data.frame(
  peak = names(peak_gr2)[q],
  gene = ann2$gene_name[s],
  gene_id = ann2$gene_id[s],
  distance = mcols(hits)$distance,
  stringsAsFactors = FALSE
)

peak_to_gene <- peak_to_gene[order(peak_to_gene$peak, peak_to_gene$distance), ]
peak_to_gene <- peak_to_gene[!duplicated(peak_to_gene$peak), ]

da_peaks_annot <- cbind(
  da_peaks,
  peak_to_gene[
    match(rownames(da_peaks), peak_to_gene$peak),
    c("gene", "gene_id", "distance"),
    drop = FALSE
  ]
)

save_csv(
  da_peaks_annot,
  file.path(atac_dir, "ActivatedMonocyte_pop2_vs_other_DA_peaks_annotated.csv"),
  row_names = TRUE
)

# ============================================================
# 19. AP-1 MOTIF CONTENT OF DA PEAKS
# ============================================================

motif_matrix <- GetMotifData(
  object = act,
  assay = "ATAC",
  slot = "data"
)

ap1_ids <- names(motif_names)[
  grepl(
    AP1_PEAK_PATTERN,
    motif_names,
    ignore.case = TRUE
  )
]

ap1_ids <- intersect(ap1_ids, colnames(motif_matrix))

if (length(ap1_ids) > 0) {
  ap1_peak_flag <- rowSums(
    motif_matrix[, ap1_ids, drop = FALSE]
  ) > 0

  da_peaks$AP1 <- ap1_peak_flag[rownames(da_peaks)]

  ap1_da_peaks <- da_peaks[
    da_peaks$AP1 %in% TRUE &
      da_peaks$p_val_adj < PEAK_PADJ_CUTOFF,
    ,
    drop = FALSE
  ]

  save_csv(
    ap1_da_peaks,
    file.path(atac_dir, "ActivatedMonocyte_pop2_significant_AP1_DA_peaks.csv"),
    row_names = TRUE
  )

  top_n <- min(1000L, nrow(da_peaks))
  top_peak_names <- rownames(
    da_peaks[order(-abs(da_peaks$avg_log2FC)), , drop = FALSE]
  )[seq_len(top_n)]

  fisher_tab <- table(
    top = rownames(da_peaks) %in% top_peak_names,
    AP1 = da_peaks$AP1
  )

  if (all(dim(fisher_tab) == c(2, 2))) {
    fisher_ap1 <- fisher.test(fisher_tab)

    fisher_out <- data.frame(
      top_n = top_n,
      prop_AP1_top = mean(da_peaks$AP1[rownames(da_peaks) %in% top_peak_names], na.rm = TRUE),
      prop_AP1_all = mean(da_peaks$AP1, na.rm = TRUE),
      odds_ratio = unname(fisher_ap1$estimate),
      conf_low = fisher_ap1$conf.int[1],
      conf_high = fisher_ap1$conf.int[2],
      p_value = fisher_ap1$p.value
    )

    save_csv(
      fisher_out,
      file.path(atac_dir, "ActivatedMonocyte_AP1_top1000_peak_Fisher_test.csv")
    )
  }
}

# ============================================================
# 20. RNA DE + HALLMARK GSEA: POP2 VS OTHER
# ============================================================

cat("\n============================================================\n")
cat("RNA differential expression and Hallmark GSEA: pop2 vs other\n")
cat("============================================================\n")

DefaultAssay(act) <- "RNA"
Idents(act) <- "pop2"

de_genes <- FindMarkers(
  act,
  ident.1 = "pop2",
  ident.2 = "other",
  logfc.threshold = 0,
  min.pct = 0.1,
  test.use = "wilcox"
)

de_genes$gene <- rownames(de_genes)

save_csv(
  de_genes,
  file.path(pathway_dir, "ActivatedMonocyte_pop2_vs_other_RNA_DE.csv")
)

ranks <- de_genes$avg_log2FC
names(ranks) <- de_genes$gene
ranks <- ranks[!is.na(ranks)]
ranks <- ranks[!duplicated(names(ranks))]
ranks <- sort(ranks, decreasing = TRUE)

hallmark_all <- split(hallmark$gene_symbol, hallmark$gs_name)
hallmark_all <- lapply(hallmark_all, unique)

set.seed(SEED)
fg <- fgsea(
  pathways = hallmark_all,
  stats = ranks,
  nperm = GSEA_NPERM,
  minSize = GSEA_MIN_SIZE,
  maxSize = GSEA_MAX_SIZE
)

fg <- fg[order(fg$padj), ]

save_fgsea_csv(
  fg,
  file.path(pathway_dir, "ActivatedMonocyte_pop2_vs_other_Hallmark_GSEA.csv")
)

# ============================================================
# 21. SAVE FINAL OBJECTS / PARAMETERS / SESSION INFO
# ============================================================

saveRDS(
  act,
  file.path(outdir, "ActivatedMonocyte_pre_with_pop2_chromVAR.rds")
)

parameter_table <- data.frame(
  parameter = c(
    "pre_rds",
    "post_rds",
    "fragment_root",
    "cell_type",
    "seed",
    "shared_top_peaks_per_timepoint",
    "global_peak_detection_fraction",
    "RNA_nfeatures",
    "PCA_npcs",
    "PCA_dims_used",
    "ATAC_top_features",
    "LSI_npcs",
    "LSI_dims_used",
    "joint_embedding",
    "Milo_k",
    "Milo_prop",
    "Milo_refined",
    "min_cells_per_sample",
    "min_samples_per_group",
    "SpatialFDR_cutoff",
    "Milo_reference",
    "Milo_contrast",
    "Hallmark_neighbourhood_pathways",
    "Hallmark_kmeans_centers",
    "legacy_TF_cluster",
    "motif_database",
    "peak_DA_test",
    "peak_DA_latent_var"
  ),
  value = c(
    pre_rds,
    post_rds,
    fragment_root,
    CELL_TYPE,
    SEED,
    SHARED_TOP_PER_TIMEPOINT,
    GLOBAL_MIN_DETECTION_FRAC,
    RNA_NFEATURES,
    PCA_NPCS,
    paste(PCA_DIMS_USE, collapse = ","),
    ATAC_TOP,
    LSI_NPCS,
    paste(LSI_DIMS_USE, collapse = ","),
    "z(PCA[1:30]) + z(LSI[2:30])",
    K_PARAM,
    PROP,
    REFINED,
    MIN_CELLS_PER_SAMPLE,
    MIN_SAMPLES_PER_GROUP,
    SPATIAL_FDR_CUTOFF,
    "NonResponder",
    "responseResponder",
    paste(PATHWAY_NAMES, collapse = ";"),
    KMEANS_CENTERS,
    LEGACY_TF_CLUSTER,
    "JASPAR2020",
    "LR",
    "nCount_ATAC"
  ),
  stringsAsFactors = FALSE
)

save_csv(
  parameter_table,
  file.path(outdir, "07_milor_pre_activated_monocytes_parameters.csv")
)

sink(file.path(outdir, "07_milor_pre_activated_monocytes_sessionInfo.txt"))
print(sessionInfo())
sink()

cat("\n============================================================\n")
cat("PRE-BCG ACTIVATED MONOCYTE MILOR ANALYSIS COMPLETE\n")
cat("============================================================\n")
cat("Milo neighbourhoods: ", nrow(da), "\n", sep = "")
cat("NR-enriched neighbourhoods: ", n_nr_nhoods, "\n", sep = "")
cat("Responder-enriched neighbourhoods: ", n_r_nhoods, "\n", sep = "")
cat("Cells in any NR-up neighbourhood: ", length(nr_up_cells), "\n", sep = "")
cat("Legacy cluster-2 neighbourhoods: ", length(pop2_nhoods), "\n", sep = "")
cat("Legacy cluster-2 cells: ", length(pop2_cells), "\n", sep = "")
cat("Significant AP-1 chromVAR motifs higher in pop2: ", nrow(ap1_pop2_up), "\n", sep = "")
cat("Outputs: ", normalizePath(outdir, winslash = "/", mustWork = FALSE), "\n", sep = "")
