# =============================================================================
# Mouse extrahepatic biliary fibroblasts — single-cell RNA-seq analysis
#
# Study    : Identification of distinct Meflin+ fibroblast subsets in the
#            biliary stroma reveals a therapeutic strategy for extrahepatic
#            cholangiocarcinoma
# Dataset  : GSE163777 (normal mouse bilio-vascular tree, 10x Genomics)
# Figures  : Figure 3, Supplementary Figure S1
#
# -----------------------------------------------------------------------------
# WHAT THIS SCRIPT DOES
#
#   0  configuration              1  setup and helper functions
#   2  download and build the Seurat object
#   3  quality control
#   4  doublet removal (scDblFinder)
#   5  normalisation, dimensionality reduction, clustering of all cells
#   6  cell-type annotation            <-- MANUAL DECISION 1
#   7  reference expression tables (mural / contractile / stellate genes)
#   8  in silico extraction of fibroblasts, with quantified purity
#   9  re-clustering of fibroblasts    <-- MANUAL DECISION 2
#  10  depth-controlled marker genes
#  11  transcriptional programme scores
#  12  differential expression between peribiliary and periportal fibroblasts
#  13  Gene Ontology over-representation analysis
#  14  gene set enrichment analysis
#  15  cross-species signature scoring   (not part of the published study)
#  16  human annotation audit            (not part of the published study)
#  17  provenance
#
# -----------------------------------------------------------------------------
# RUNNING PART OF THIS SCRIPT
#
# Every section is preceded by a machine-readable marker, "#@ STEP n". The
# helper file run.R reads those markers so that sections can be run by number
# instead of by hunting for line numbers:
#
#   source("run.R")
#   steps()          # list the sections
#   run_steps(9)     # CONFIG and helpers, then section 9
#   run_steps(9:16)
#
# CONFIG (section 0) and the helper functions (section 1) are always run first,
# so a change to either is picked up without being looked for. Objects a section
# needs are restored from the saved checkpoints when they are not already in the
# session. run.R is a development convenience and is not part of the analysis.
#
# -----------------------------------------------------------------------------
# THE TWO PLACES THAT NEED A HUMAN
#
# The script assigns what the evidence assigns unambiguously and STOPS wherever
# a judgement is required, having first written out the numbers needed to make
# it. Nothing ambiguous is decided automatically. Search for ">>>" to find them.
#
#   >>> MANUAL CALLS  (section 6)
#       Clusters whose panel score does not win by a clear margin, or that carry
#       markers of several lineages at once. Fill in CELLTYPE_MANUAL and
#       CELLTYPE_FINGERPRINT. Evidence: printed read-out plus
#       S1_01f_assignment_decisions, S1_01g_markers_of_flagged_clusters,
#       S1_01c_diagnostic_genes, S1_01d_multilineage_clusters, S1_01e_qc_by_cluster.
#
#   >>> FIBROBLAST CLUSTERS TO DROP  (section 9)
#       Fibroblast clusters that are ambient contamination rather than a
#       fibroblast state. Fill in FIBROBLAST_DROP. Evidence:
#       S1_05e_fibroblast_cluster_qc, S1_05f_fibroblast_cluster_markers,
#       S1_05d_cluster_separability, S1_05_resolution_scan.
#
# CELLTYPE_FINGERPRINT records the cluster sizes the manual calls were written
# against. If anything upstream changes, the fingerprint changes and the calls
# are rejected instead of being applied to differently numbered clusters. This
# is deliberate: cluster numbers are not stable, and a stale manual call is a
# silent error.
#
# -----------------------------------------------------------------------------
# REPRODUCIBILITY
#
# Seurat's functions each carry their own seed argument and ignore the global
# RNG, so RunPCA (seed.use), FindClusters (random.seed) and RunUMAP (seed.use)
# are already deterministic. scDblFinder is not: it must be given a seed
# immediately before the call and a seeded parallel back-end, otherwise a
# different set of cells is removed on every run and everything downstream
# moves. Its calls are therefore cached in 00_doublet_calls.rds -- DO NOT DELETE
# THAT FILE if you want the published numbers to reproduce.
#
# Section 6 rewrites Idents(sc) at its end, so it resets them to seurat_clusters
# at its start. Re-running it from part-way through would otherwise key the
# annotation on cell-type labels instead of cluster numbers.
#
# -----------------------------------------------------------------------------
# ANALYTICAL DECISIONS AND THEIR REASONS
#
#   Fibroblast definition.  Fibroblasts are defined by high Col1a1 and Pdgfra
#   together with absence of the mural markers Rgs5 and Notch3. Pdgfrb is NOT
#   used: in this tissue it is expressed more highly by pericytes (4.81) than by
#   fibroblasts (1.35-3.33), so it does not separate the two. Kcnj8 and Abcc9 are
#   also not used: both are expressed at low absolute levels in every cluster
#   here and do not discriminate.
#
#   No pseudotime.  The normal extrahepatic biliary stroma contains no
#   myofibroblasts (section 6 quantifies this), so no terminal state exists in
#   this dataset. Myofibroblastic conversion is assessed experimentally by
#   lineage tracing after bile duct ligation and in orthotopic tumours.
#
#   Depth control.  Clusters differ in transcriptome complexity, and genes
#   detected at trace levels appear enriched in the cluster with the greater
#   sequencing depth. Marker genes are therefore re-tested after downsampling
#   every cell to a common depth (section 9), and only depth-robust genes are
#   used for interpretation.
#
#   Ribosomal genes.  Ribosomal protein and mitochondrially encoded transcripts
#   differ systematically between clusters and dominate rank-based enrichment.
#   They are excluded before over-representation and gene set enrichment
#   analysis (sections 12-13).
#
# -----------------------------------------------------------------------------
# HOW TO RUN
#
#   Set the parameters in section 0 and run top to bottom. All outputs are
#   written under OUT_DIR; nothing is written to the working directory.
#   Cluster numbering is not stable across runs, so cluster identities are
#   assigned by marker-panel scoring rather than by index, and every assignment
#   is checked with stopifnot().
#
# Software (versions used)
#   R 4.2.0 | Seurat 4.3.0 | scDblFinder 1.12.0 | patchwork 1.1.2 |
#   clusterProfiler 4.6.0 | org.Mm.eg.db 3.16.0 | rstatix 0.7.2
# =============================================================================


#@ STEP 0
# =============================================================================
# 0. CONFIGURATION
# =============================================================================

CONFIG <- list(
  seed            = 1234,

  data_dir        = "data/GSE163777",
  out_dir         = file.path("output", "mouse_biliary_fibroblast"),
  gsm             = "GSM4987107",

  # quality control
  min_cells       = 3,
  min_features    = 200,
  max_percent_mt  = 10,
  max_features    = NA,        # NA = no upper bound; set a value to apply one

  # clustering, all cells
  n_hvg           = 2000,
  n_pcs_all       = 30,
  resolution_all  = 0.5,

  # clustering, fibroblasts
  n_pcs_fib       = 10,
  # Two clustering passes. The first is deliberately finer than the final one:
  # its job is to expose contamination as its own small cluster so it can be
  # removed, which a coarse partition would bury inside a subset. The second
  # defines the subsets that are reported.
  resolution_fib_qc = 0.3,
  resolution_fib    = 0.2,
  # resolutions compared in the separability test in section 9. Set these to the
  # values that bracket the decision: the one that gives the number of clusters
  # you expect, and the next one up.
  resolution_candidates = c(0.1, 0.2, 0.3),

  # depth control
  # NA = downsample to the median library size of the shallowest subset. Do NOT
  # use the global minimum: a single low-depth cell then sets the target for
  # every cell, detection collapses (pct.1 fell to 23 % of its native value when
  # this was tried here), and genuinely restricted markers such as Notum drop
  # below min.pct and are scored as depth artefacts. The aim is to remove the
  # depth DIFFERENCE between the groups being compared, not to reduce every
  # library to the worst one in the dataset.
  downsample_umi  = NA,
  # ADDED 2026-10-01. The downsampled object is cached in
  # objects/03b_fibroblasts_downsampled.rds and reused, because the marker
  # counts, the post-downsampling median UMI and the whole GSEA are computed
  # from it and are already quoted in the manuscript. TRUE forces a fresh draw;
  # the draw is seeded, so it is reproducible either way, but a fresh draw is a
  # DIFFERENT draw and those numbers have to be re-read from it.
  redo_downsample = FALSE,

  # ADDED 2026-10-05, same pattern and same reason as redo_downsample. The
  # one-versus-rest GSEA is now run with pvalueCutoff = 1 so that EVERY tested
  # set keeps its NES -- a panel that shows a set as "not enriched" needs the
  # non-significant value, and the saved table only ever held the significant
  # ones. Benjamini-Hochberg is applied across all tested sets before the
  # cutoff is used, so the adjusted P values and therefore the significant
  # table are unchanged by this; what changes is only how much is kept.
  # The full result is cached, and the significant table is FILTERED from the
  # cache rather than produced by a second run, so the two can never disagree.
  redo_gsea       = FALSE,

  # enrichment
  exclude_pattern = "^Rp[sl]|^Mrp[sl]|^mt-",   # ribosomal and mitochondrial
  min_cells_gene  = 10,        # gene must be detected in >= this many cells
                               # to enter the enrichment universe

  # MSigDB gene sets for section 14. msigdb_dir is the directory holding the
  # downloaded GMT files; NULL skips the MSigDB analysis entirely. Use the
  # *.symbols.gmt files -- the ranked list is keyed on gene symbols, and an
  # entrez GMT silently matches nothing.
  #
  # Collections are PRE-SPECIFIED here rather than chosen after seeing the
  # results. Each collection is its own testing family (Benjamini-Hochberg is
  # applied within, not across), so adding collections inflates the chance of a
  # significant hit somewhere. Report every collection that was run, including
  # the empty ones.
  msigdb_dir         = NULL,
  msigdb_collections = c(
    hallmark     = "mh.all.v2024.1.Mm.symbols.gmt",
    reactome     = "m2.cp.reactome.v2024.1.Mm.symbols.gmt",
    wikipathways = "m2.cp.wikipathways.v2024.1.Mm.symbols.gmt"),
    # Also available in the 2024.1.Mm download, deliberately not run by default:
    #   m5.go.bp.*        duplicates gseGO above with a different annotation
    #                     snapshot, so the two disagree slightly and invite
    #                     reporting whichever looks better
    #   m2.cgp.*          perturbation signatures; thousands of sets, mostly
    #                     from unrelated systems
    #   m3.gtrd.*         transcription-factor target sets; useful if the
    #                     question is which regulators differ, but pre-specify
    #   m8.all.*          cell-type signatures; informative for annotation,
    #                     not for a within-fibroblast contrast
    #   msigdb.*          everything at once; do not use as a fishing net

  # retained for compatibility with the working copy; the optional
  # cross-species sections are not part of the published analysis
  human_rds       = NULL,
  human_group_by  = NULL,
  human_caf_rds      = NULL,
  human_caf_group_by = NULL
)

# Machine-specific paths live in config_local.R, which is not part of the
# released code. Keep it out of the repository (add it to .gitignore) so that
# the script published alongside the paper carries only relative paths.
#
#   config_local.R
#     CONFIG$data_dir  <- "/path/to/GSE163777"
#     CONFIG$out_dir   <- file.path("output", "run_YYYYMMDD")
#
if (file.exists("config_local.R")) {
  source("config_local.R")
  message("Local configuration applied from config_local.R")
}

set.seed(CONFIG$seed)


#@ STEP 1
# =============================================================================
# 1. SETUP
# =============================================================================

library(Seurat)
library(Matrix)
library(tidyverse)
library(patchwork)
library(SingleCellExperiment)
library(scDblFinder)
library(clusterProfiler)
library(org.Mm.eg.db)
library(rstatix)

OUT_DIR <- CONFIG$out_dir
for (sub in c("figures", "tables", "objects"))
  dir.create(file.path(OUT_DIR, sub), recursive = TRUE, showWarnings = FALSE)

out <- function(...) file.path(OUT_DIR, ...)
fig <- function(x)   out("figures", x)
tab <- function(x)   out("tables",  x)
obj <- function(x)   out("objects", x)

# Diagnostic plots from this script go to figures/exploratory/, NOT to figures/,
# and every one of them is named explor_*. RENAMED 2026-09-29: eight of them
# still carried manuscript names (Fig5A_umap_fibroblasts, Fig5B_featureplots,
# Fig5C_violin, Fig5D_marker_heatmap, Fig5E_programme_scores,
# FigS1A/B/C_*). The directory split already stopped them overwriting the real
# panels, but the panel letters have since been reassigned -- Fig. 5D is now
# the marker dot plot, 5E the GO panel -- so a file called Fig5D_marker_heatmap
# sitting in the output was an invitation to put the wrong plot in the paper.
#
# Both this script and figures.R used to write into figures/, under fourteen
# identical file names (FigS1A_umap_allcells, Fig5B_featureplots, FigS1D_GO and
# so on). Whichever ran last won. Running run_steps() after figures.R therefore
# replaced every publication figure with this script's unstyled version --
# no jittered points, no rasterised layers, default theme, default font -- and
# nothing said so, because overwriting a file is not an error.
#
# The two scripts have different jobs. This one checks the analysis as it goes
# and its plots are for the person running it. figures.R produces what goes to
# the journal. Separating the directories is what makes that division real
# rather than a convention that a re-run can silently break.
save_fig <- function(p, name, width = 6, height = 4) {
  d <- out("figures", "exploratory")
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(file.path(d, paste0(name, ".pdf")), p,
                  width = width, height = height)
  invisible(p)
}
save_tab <- function(x, name, row.names = FALSE) {
  utils::write.csv(x, tab(paste0(name, ".csv")), row.names = row.names)
  invisible(x)
}

# Seurat 4 uses slot=, Seurat 5 uses layer=.
get_counts <- function(x)
  tryCatch(GetAssayData(x, layer = "counts"),
           error = function(e) GetAssayData(x, slot = "counts"))

# Genes absent from the matrix must be dropped before indexing.
present <- function(genes, x) intersect(genes, rownames(x))

# Seurat 5 prefixes numeric identity names with "g"; strip it so that column
# names match the values held in the metadata.
#
# group.by is passed only when it is set. Seurat's default is "ident", and
# passing group.by = NULL explicitly resolves to object[[NULL]], which fails
# with "undefined columns selected".
avg_expr <- function(x, features, group.by = NULL) {
  features <- present(features, x)
  if (!length(features))
    stop("none of the requested features are present in the object")
  args <- list(object = x, features = features, assays = "RNA")
  if (!is.null(group.by)) args$group.by <- group.by
  m <- do.call(AverageExpression, args)$RNA
  m <- as.matrix(m)
  colnames(m) <- sub("^g", "", colnames(m))
  m
}

# Mean z-scored expression of each marker panel in each cluster. Used to assign
# identities by evidence rather than by cluster index, which is not stable.
# Genes absent from the object are dropped and reported; a panel left with no
# genes scores NA rather than silently scoring zero.
# Mean z-scored expression of each marker panel in each cluster. Identities are
# assigned by evidence rather than by cluster index, which is not stable.
#
# Two filters keep a panel from producing a spurious winner:
#
#   min_expr   a gene not expressed anywhere still gets z-scored and contributes
#              noise to whichever panel contains it
#   min_genes  a panel reduced to a single gene always yields a z-score near 1
#              for whichever cluster is highest, however low the absolute
#              expression. A cell type is called only when at least this many of
#              its markers are expressed somewhere in the dataset.
#
# A panel that falls below min_genes scores NA and cannot win. The genes it lost
# are themselves a result: if none of the canonical markers of a lineage clear
# the floor, that lineage is not present in the tissue.
score_panels <- function(x, panels, group.by = NULL, min_expr = 1.0,
                         min_genes = 2, verbose = TRUE) {

  wanted  <- unique(unlist(panels))
  missing <- setdiff(wanted, rownames(x))

  m <- avg_expr(x, wanted, group.by = group.by)
  silent <- rownames(m)[apply(m, 1, max) < min_expr]
  usable <- setdiff(rownames(m), silent)
  m <- m[usable, , drop = FALSE]

  status <- tibble::tibble(
    panel      = names(panels),
    n_total    = lengths(panels),
    n_usable   = vapply(panels, function(g) length(intersect(g, usable)), integer(1)),
    called     = n_usable >= min_genes,
    not_in_object = vapply(panels, function(g)
                      paste(intersect(g, missing), collapse = ", "), character(1)),
    below_floor   = vapply(panels, function(g)
                      paste(intersect(g, silent), collapse = ", "), character(1)))

  if (verbose) {
    print(status, n = Inf)
    if (any(!status$called))
      message("Lineages not called, fewer than ", min_genes,
              " markers expressed above ", min_expr, ": ",
              paste(status$panel[!status$called], collapse = ", "))
  }

  z <- t(scale(t(log1p(m))))          # genes with no variance become NaN
  z[!is.finite(z)] <- 0

  res <- sapply(panels, function(g) {
    g <- intersect(g, rownames(z))
    if (length(g) < min_genes) return(rep(NA_real_, ncol(z)))
    colMeans(z[g, , drop = FALSE])
  })
  res <- as.matrix(res)
  rownames(res) <- colnames(m)
  attr(res, "panel_status") <- status
  res
}

## --- composition-driven structure in the embedding ---------------------------
# Ribosomal protein transcripts are among the most abundant in the library, are
# tightly co-regulated as a block of roughly 80 genes, and differ systematically
# in fraction between cell states. If they enter the variable-feature set they
# can define a principal component on their own, and any PC inside the range
# used for the neighbour graph then contributes to the clustering. A cluster
# boundary drawn on that axis separates cells by library composition rather than
# by identity.
#
# This is not hypothetical here: ribosomal transcripts are the single largest
# block of genes separating the peribiliary from the periportal subset, and the
# difference persists after UMI downsampling. The following two functions
# measure the problem and test whether the clustering depends on it.

# Association between each principal component and the technical covariates.
# Spearman is used because the relationship need not be linear.
pc_covariate_cor <- function(x, dims, label) {
  emb  <- Embeddings(x, "pca")[, dims, drop = FALSE]
  covs <- data.frame(percent_ribo = x$percent.ribo,
                     percent_mt   = x$percent.mt,
                     log_nCount   = log10(x$nCount_RNA))
  out <- expand.grid(PC = dims, covariate = names(covs),
                     stringsAsFactors = FALSE) %>%
    mutate(rho = round(mapply(function(p, cv)
             suppressWarnings(cor(emb[, paste0("PC_", p)], covs[[cv]],
                                  method = "spearman")),
             PC, covariate), 3)) %>%
    arrange(desc(abs(rho)))
  save_tab(out, paste0("QC_pc_covariate_correlation_", label))
  flagged <- out %>% filter(covariate == "percent_ribo", abs(rho) > 0.4)
  if (nrow(flagged))
    warning(label, ": PC ", paste(flagged$PC, collapse = ", "),
            " correlate with ribosomal content (rho up to ",
            round(max(abs(flagged$rho)), 2),
            "); clustering may be partly compositional")
  print(head(out, 10))
  out
}

# Adjusted Rand index, written out rather than taken from mclust to keep the
# dependency list short. 1 = identical partitions, 0 = chance agreement.
adj_rand <- function(a, b) {
  tab <- table(a, b); n <- sum(tab)
  cs  <- function(v) sum(choose(v, 2))
  i   <- cs(rowSums(tab)); j <- cs(colSums(tab))
  idx <- cs(as.vector(tab)); expected <- i * j / choose(n, 2)
  (idx - expected) / ((i + j) / 2 - expected)
}

# Re-run variable features, PCA and clustering with ribosomal and mitochondrial
# transcripts removed from the variable-feature set only (they stay in the
# object and in every expression value reported). Returns the alternative
# cluster labels, aligned to the cells of x.
cluster_without_housekeeping <- function(x, n_hvg, dims, resolution, seed) {
  y <- FindVariableFeatures(x, selection.method = "vst",
                            nfeatures = n_hvg + 500, verbose = FALSE)
  vf <- setdiff(VariableFeatures(y), grep(CONFIG$exclude_pattern,
                                          rownames(y), value = TRUE))
  vf <- head(vf, n_hvg)
  VariableFeatures(y) <- vf
  y <- RunPCA(y, features = vf, npcs = max(dims), seed.use = seed,
              verbose = FALSE)
  y <- FindNeighbors(y, dims = dims, verbose = FALSE)
  y <- FindClusters(y, resolution = resolution, random.seed = seed,
                    verbose = FALSE)
  Idents(y)
}
## --- lineage sentinels --------------------------------------------------------
# Markers of mutually exclusive lineages. A cluster carrying several of them at
# once is residual ambient RNA, unresolved doublets or low-complexity cells
# rather than a cell type; it is flagged for adjudication, never auto-assigned.
#
# These live in section 1, with the helpers, rather than next to the annotation
# logic in section 9 that first uses them, because section 9 uses them AS WELL.
# run.R always executes sections 0 and 1 and only the sections you name, so a
# constant defined in section 6 and read in section 9 makes run_steps(9) fail in
# a fresh session with "object 'LINEAGE_SENTINELS' not found". Anything used by
# more than one section belongs here.
LINEAGE_SENTINELS <- list(hepatocyte  = "Alb",
                          endothelium = c("Cdh5", "Pecam1"),
                          mural       = "Rgs5",
                          glia        = "Plp1",
                          fibroblast  = "Col1a1",
                          epithelium  = "Epcam",
                          immune      = "Ptprc")
SENTINEL_FLOOR <- c(Alb = 3, Cdh5 = 1, Pecam1 = 0.5, Rgs5 = 1, Plp1 = 0.5,
                    Col1a1 = 5, Epcam = 0.5, Ptprc = 0.5)
MAX_LINEAGES   <- 3     # more than this and the cluster must be adjudicated

## -----------------------------------------------------------------------------

message("Output directory: ", normalizePath(OUT_DIR, mustWork = FALSE))


#@ STEP 2
# =============================================================================
# 2. DOWNLOAD AND OBJECT CONSTRUCTION
# =============================================================================
# GSE163777 comprises a single 10x sample.

if (!dir.exists(CONFIG$data_dir)) {
  # BiocManager::install("GEOquery")
  GEOquery::getGEOSuppFiles("GSE163777", baseDir = dirname(CONFIG$data_dir),
                            makeDirectory = TRUE)
}

mat      <- readMM(file.path(CONFIG$data_dir,
                             paste0(CONFIG$gsm, "_matrix.mtx.gz")))
features <- read.delim(file.path(CONFIG$data_dir,
                                 paste0(CONFIG$gsm, "_features.tsv.gz")),
                       header = FALSE, stringsAsFactors = FALSE)
barcodes <- read.delim(file.path(CONFIG$data_dir,
                                 paste0(CONFIG$gsm, "_barcodes.tsv.gz")),
                       header = FALSE, stringsAsFactors = FALSE)

colnames(mat) <- barcodes$V1
rownames(mat) <- features$V1
keep_symbol   <- !duplicated(features$V2)       # gene symbols where unique
rownames(mat)[keep_symbol] <- features$V2[keep_symbol]

sc <- CreateSeuratObject(counts = mat, project = "GSE163777",
                         min.cells    = CONFIG$min_cells,
                         min.features = CONFIG$min_features)
sc[["percent.mt"]] <- PercentageFeatureSet(sc, pattern = "^mt-")

# Ribosomal fraction. Recorded as a covariate, not as a filter: it is used below
# to test whether the clustering is driven by library composition. No cell is
# excluded on this metric.
sc[["percent.ribo"]] <- PercentageFeatureSet(sc, pattern = "^Rp[sl]")

n_raw <- ncol(sc)


#@ STEP 3
# =============================================================================
# 3. QUALITY CONTROL
# =============================================================================

save_fig(VlnPlot(sc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                 ncol = 3, pt.size = 0), "QC_01_prefilter", 9, 3)

sc <- subset(sc, subset = nFeature_RNA > CONFIG$min_features &
                          percent.mt   < CONFIG$max_percent_mt)
if (!is.na(CONFIG$max_features))
  sc <- subset(sc, subset = nFeature_RNA < CONFIG$max_features)

n_qc <- ncol(sc)
message(sprintf("Cells: %d raw -> %d after QC", n_raw, n_qc))


#@ STEP 4
# =============================================================================
# 4. DOUBLET REMOVAL
# -----------------------------------------------------------------------------
# Run on all cells before any subsetting, so that doublets formed between any
# pair of cell types can be detected. The same tool and version are used for the
# human dataset, so the two analyses are treated identically.
# =============================================================================

# scDblFinder samples randomly and is NOT reproducible unless the seed is set
# immediately before the call and the parallel back-end is given its own seed.
# Setting the seed once at the top of the script is not enough: any intervening
# call that touches the RNG shifts the state, and a different set of cells is
# then removed. That propagates into every downstream result, including how many
# clusters a given resolution produces.
#
# The calls are cached. On a re-run the file is reused, so the analysis is
# reproducible even if the package version changes. Delete the file to recompute.
dbl_file <- obj("00_doublet_calls.rds")
if (file.exists(dbl_file)) {
  dbl <- readRDS(dbl_file)
  stopifnot(setequal(names(dbl$class), colnames(sc)))
  message("Reusing cached doublet calls from ", dbl_file)
} else {
  set.seed(CONFIG$seed)
  sce_all <- as.SingleCellExperiment(sc)
  sce_all <- scDblFinder(
    sce_all,
    BPPARAM = BiocParallel::SerialParam(RNGseed = CONFIG$seed))
  stopifnot(identical(colnames(sce_all), colnames(sc)))
  dbl <- list(class = setNames(as.character(sce_all$scDblFinder.class),
                               colnames(sce_all)),
              score = setNames(as.numeric(sce_all$scDblFinder.score),
                               colnames(sce_all)))
  saveRDS(dbl, dbl_file)
  rm(sce_all)
}
sce_all <- NULL

sc$doublet_class <- unname(dbl$class[colnames(sc)])
sc$doublet_score <- unname(dbl$score[colnames(sc)])

save_tab(as.data.frame(table(class = sc$doublet_class)), "QC_02_doublet_counts")
message(sprintf("Doublets called: %d / %d (%.1f%%)",
                sum(sc$doublet_class == "doublet"), ncol(sc),
                100 * mean(sc$doublet_class == "doublet")))

sc <- subset(sc, subset = doublet_class == "singlet")
n_singlet <- ncol(sc)


#@ STEP 5
# =============================================================================
# 5. NORMALISATION, DIMENSIONALITY REDUCTION, CLUSTERING (ALL CELLS)
# =============================================================================

sc <- NormalizeData(sc, normalization.method = "LogNormalize",
                    scale.factor = 1e4)
sc <- FindVariableFeatures(sc, selection.method = "vst",
                           nfeatures = CONFIG$n_hvg)
sc <- ScaleData(sc, features = rownames(sc))
sc <- RunPCA(sc, features = VariableFeatures(sc), npcs = 50)

save_fig(ElbowPlot(sc, ndims = 50), "QC_03_elbow_allcells", 4, 3)

dims_all <- seq_len(CONFIG$n_pcs_all)
sc <- FindNeighbors(sc, dims = dims_all)
sc <- FindClusters(sc, resolution = CONFIG$resolution_all)
sc <- RunUMAP(sc, dims = dims_all)

message("Clusters (all cells): ", nlevels(Idents(sc)))

# Composition check on the all-cell embedding. Cell-type signals here are much
# stronger than any compositional axis, so this is expected to pass; it is run
# because the fibroblast extraction depends on it.
message(sprintf("Ribosomal/mitochondrial genes among the %d variable features: %d",
                length(VariableFeatures(sc)),
                sum(grepl(CONFIG$exclude_pattern, VariableFeatures(sc)))))
pc_cor_all <- pc_covariate_cor(sc, dims_all, "allcells")

alt_all <- cluster_without_housekeeping(sc, CONFIG$n_hvg, dims_all,
                                        CONFIG$resolution_all, CONFIG$seed)
ari_all <- adj_rand(Idents(sc), alt_all)
message(sprintf("All-cell clustering without ribosomal/mitochondrial HVGs: %d clusters, ARI = %.3f",
                nlevels(alt_all), ari_all))
save_tab(as.data.frame.matrix(table(original = Idents(sc), no_housekeeping = alt_all)),
         "QC_04b_clustering_without_housekeeping_allcells", row.names = TRUE)

saveRDS(sc, obj("01_allcells_clustered.rds"))


#@ STEP 6
# =============================================================================
# 6. CELL-TYPE ANNOTATION
# -----------------------------------------------------------------------------
# Identities are assigned by marker-panel scoring, then checked. Assignment is
# made on mean expression per cluster rather than on detection rate: ambient RNA
# inflates detection rate for abundant transcripts (for example Alb is detected
# in 33-56% of cells in every cluster while its mean expression is two orders of
# magnitude higher in hepatocytes than elsewhere).
#
# This section must be re-runnable. The end of it replaces Idents(sc) with the
# cell-type labels, so running it a second time in the same session would key the
# panel scores, the annotation table and the manual calls on labels instead of
# cluster numbers -- and the manual calls would then silently address the wrong
# things. Resetting the identities here makes the section idempotent.
# =============================================================================

if (is.null(sc$seurat_clusters))
  stop("sc$seurat_clusters is missing; run section 5 first")
Idents(sc) <- sc$seurat_clusters
message("Identities reset to seurat_clusters (", nlevels(Idents(sc)),
        " clusters) for annotation")

# Panels must not share genes between lineages that need to be told apart.
# Mesothelium expresses Krt19 and Krt7, so those two are kept out of the
# cholangiocyte panel, which relies on Epcam, Sox9, Cftr and Hnf1b instead.
# Hepatic stellate cells are Col1a1-positive, so the fibroblast and stellate
# panels are separated by Lrat, Reln, Lama1 and Ecm1 rather than by collagen.

PANELS <- list(
  Fibroblast    = c("Col1a1", "Col1a2", "Pdgfra", "Dcn", "Lum", "Dpt"),
  HSC           = c("Lrat", "Reln", "Lama1", "Ecm1"),
  Pericyte      = c("Rgs5", "Notch3", "Cspg4", "Higd1b", "Vtn"),
  VSM           = c("Myh11", "Cnn1", "Pln", "Des", "Acta2"),
  sEC           = c("Stab2", "Clec4g", "Gata4", "Oit3"),
  vEC           = c("Vwf", "Tek", "Pecam1", "Cdh5", "Nr2f2"),
  aEC           = c("Sema3g", "Ly6c1", "Gja5", "Gja4", "Efnb2", "Bmx"),
  lyEC          = c("Lyve1", "Prox1", "Pdpn", "Ccl21a"),
  Mesothelium   = c("Msln", "Upk3b", "Lrrn4", "Lrrn4cl", "Wt1"),
  Hepatocyte    = c("Alb", "Fabp1", "Apoa5", "Rgn", "Ttr"),
  Cholangiocyte = c("Epcam", "Sox9", "Cftr", "Hnf1b", "Onecut1", "Spp1"),
  Immune        = c("Ptprc", "Cd3e", "Cd68", "Cd14", "Lyz2", "C1qa", "C1qb"),
  Schwann       = c("Sox10", "Plp1", "Mpz", "S100b", "Ncmap", "Mbp")
)

# Genes that decide the ambiguous calls, printed next to the panel scores so
# that a marginal assignment can be adjudicated directly rather than by score.
DIAGNOSTIC_GENES <- c("Col1a1", "Pdgfra", "Pdgfrb", "Islr",
                      "Lrat", "Reln",
                      "Rgs5", "Notch3", "Myh11", "Cnn1",
                      "Pecam1", "Cdh5",
                      "Msln", "Upk3b", "Krt19", "Epcam", "Sox9",
                      "Sox10", "Plp1",
                      "Alb", "Ptprc")

panel_scores <- score_panels(sc, PANELS, min_expr = 1.0, min_genes = 2)
save_tab(attr(panel_scores, "panel_status"), "S1_01a_panel_status")
save_tab(as.data.frame(round(panel_scores, 3)),
         "S1_01_celltype_panel_scores", row.names = TRUE)
print(round(panel_scores, 2))

# A cluster with no clear winner is low quality: it carries a little of every
# lineage and nothing specific.
# A panel scored NA (no genes present) must not win, and must not break the
# ranking. sort() drops NA, so the margin falls back to NA when only one panel
# could be scored.
best      <- apply(panel_scores, 1, function(v)
                   colnames(panel_scores)[which.max(replace(v, is.na(v), -Inf))])
best_val  <- apply(panel_scores, 1, max, na.rm = TRUE)
second    <- apply(panel_scores, 1, function(v) {
                   v <- sort(v, decreasing = TRUE)
                   if (length(v) >= 2) v[2] else NA_real_ })
margin    <- best_val - second

annotation <- tibble(cluster = rownames(panel_scores),
                     best = best, score = round(best_val, 3),
                     margin = round(margin, 3)) %>%
  mutate(flag = ifelse(margin < 0.5 | is.na(margin), "REVIEW", ""))
save_tab(annotation, "S1_01b_celltype_assignment")
print(annotation, n = Inf)

diagnostics <- avg_expr(sc, DIAGNOSTIC_GENES) %>% round(2) %>% as.data.frame()
save_tab(diagnostics, "S1_01c_diagnostic_genes", row.names = TRUE)
print(diagnostics)

## =============================================================================
## DECISION POINT -- automatic where the evidence is unambiguous, manual otherwise
## -----------------------------------------------------------------------------
## Cell-type assignment is the step every later result depends on, and it is the
## step where a silent error is least likely to be noticed: a mislabelled cluster
## does not throw an error, it just changes the answer. Cluster numbers are also
## not stable -- they have changed twice in this analysis -- so anything keyed on
## a number goes stale without saying so.
##
## The design is therefore:
##   * a cluster is assigned automatically ONLY if its panel score wins by a wide
##     margin and it shows no sign of being a mixture;
##   * every other cluster is left unassigned, its evidence is printed and
##     written to disk, and the script STOPS;
##   * the user enters those calls in CELLTYPE_MANUAL below;
##   * on the next run the script verifies that the manual block was written for
##     this clustering (via a fingerprint) and that it covers exactly the
##     clusters that need it, and stops again if either has changed.
##
## Nothing about an ambiguous cluster is decided by the script.
## =============================================================================

AUTO_MARGIN    <- 0.5   # panel score must beat the runner-up by this much
AUTO_MIN_SCORE <- 0.5   # and must itself be at least this high

## --- evidence 1: markers of several mutually exclusive lineages at once --------
# LINEAGE_SENTINELS, SENTINEL_FLOOR and MAX_LINEAGES are defined in section 1,
# because section 9 uses them too and the runner only guarantees sections 0 and
# 1. A constant defined here and used there breaks every partial re-run.

# A cluster can trip the multi-lineage test purely through ambient RNA: a small,
# transcriptionally enormous population such as hepatocytes carries a little of
# every other lineage, and every other cluster carries a little Alb. What
# separates that from a genuine mixture is DOMINANCE -- how far the winning panel
# is ahead, and how far the leading sentinel is ahead of the next one. A cluster
# that is far ahead on both is not ambiguous whatever the sentinel count says, so
# the flag can be waived. Both thresholds are deliberately strict.
MULTILINEAGE_WAIVER <- list(min_margin = 2.0, min_dominance = 10)

# Cluster numbers must come from a source that later code cannot overwrite.
# Idents(sc) is reassigned to the cell-type labels further down, so re-running
# this section in the same session would otherwise silently key everything on
# labels instead of numbers -- which is how the cluster fingerprint came out with
# 12 entries for a 16-cluster run.
cluster_vec <- if (!is.null(sc$seurat_clusters))
  as.character(sc$seurat_clusters) else as.character(Idents(sc))
if (!setequal(unique(cluster_vec), annotation$cluster))
  stop("the annotation table is keyed on ", paste(head(annotation$cluster, 3),
       collapse = ", "), " ... but the clusters are ",
       paste(head(sort(unique(cluster_vec)), 3), collapse = ", "),
       " ... . Section 6 was run against cell-type labels rather than cluster ",
       "numbers, which happens when it is re-run from part-way through. ",
       "Re-run section 6 from its start.")

sentinel_expr <- avg_expr(sc, unique(unlist(LINEAGE_SENTINELS)),
                          group.by = "seurat_clusters")
lineages_present <- function(k)
  names(which(vapply(LINEAGE_SENTINELS, function(g) {
    g <- intersect(g, rownames(sentinel_expr))
    length(g) > 0 && any(sentinel_expr[g, k] >= SENTINEL_FLOOR[g])
  }, logical(1))))

multilineage <- tibble(
  cluster    = colnames(sentinel_expr),
  n_lineages = vapply(colnames(sentinel_expr),
                      function(k) length(lineages_present(k)), integer(1)),
  lineages   = vapply(colnames(sentinel_expr),
                      function(k) paste(lineages_present(k), collapse = "+"),
                      character(1)),
  # ratio of the leading sentinel to the next one, on the floor-normalised scale
  # so that genes with different floors are comparable
  dominance  = vapply(colnames(sentinel_expr), function(k) {
                 r <- sort(sentinel_expr[, k] /
                           SENTINEL_FLOOR[rownames(sentinel_expr)],
                           decreasing = TRUE)
                 if (length(r) < 2 || r[2] <= 0) Inf else round(r[1] / r[2], 1)
               }, numeric(1)),
  values     = vapply(colnames(sentinel_expr), function(k)
                 paste(sprintf("%s %.2f", rownames(sentinel_expr),
                               sentinel_expr[, k]), collapse = "; "),
                 character(1)))
save_tab(multilineage, "S1_01d_multilineage_clusters")
print(multilineage %>% select(-values), n = Inf)

## --- evidence 2: QC profile ---------------------------------------------------
if (is.null(sc$percent.ribo))
  sc[["percent.ribo"]] <- PercentageFeatureSet(sc, pattern = "^Rp[sl]")
if (is.null(sc$doublet_score)) sc$doublet_score <- NA_real_

qc_by_cluster <- tibble(cluster = cluster_vec,
                        nFeature = sc$nFeature_RNA, nCount = sc$nCount_RNA,
                        percent_mt = sc$percent.mt,
                        percent_ribo = sc$percent.ribo,
                        doublet_score = sc$doublet_score) %>%
  group_by(cluster) %>%
  summarise(n_cells = n(),
            nFeature_median = median(nFeature),
            nCount_median   = median(nCount),
            percent_mt_median   = round(median(percent_mt), 2),
            percent_ribo_median = round(median(percent_ribo), 2),
            doublet_score_median = round(median(doublet_score), 3),
            .groups = "drop")
save_tab(qc_by_cluster, "S1_01e_qc_by_cluster")

## --- which clusters can be assigned automatically ------------------------------
ml <- multilineage[match(annotation$cluster, multilineage$cluster), ]

# The multi-lineage flag is waived where the assignment is far ahead on both the
# panel margin and the sentinel dominance, i.e. where the extra lineages are
# clearly ambient rather than a mixture.
waived <- ml$n_lineages > MAX_LINEAGES &
          annotation$margin >= MULTILINEAGE_WAIVER$min_margin &
          ml$dominance      >= MULTILINEAGE_WAIVER$min_dominance
waived[is.na(waived)] <- FALSE

auto_ok <- annotation$margin >= AUTO_MARGIN &
           annotation$score  >= AUTO_MIN_SCORE &
           !is.na(annotation$margin) &
           (ml$n_lineages <= MAX_LINEAGES | waived)
if (any(waived))
  message("Multi-lineage flag waived (ambient, not a mixture) for cluster(s): ",
          paste(annotation$cluster[waived], collapse = ", "))

CELLTYPE <- setNames(rep(NA_character_, nrow(annotation)), annotation$cluster)
CELLTYPE[auto_ok] <- annotation$best[auto_ok]
needs_decision <- annotation$cluster[!auto_ok]

decision_table <- annotation %>%
  left_join(multilineage %>% select(cluster, n_lineages, lineages, dominance),
            by = "cluster") %>%
  left_join(qc_by_cluster, by = "cluster") %>%
  mutate(assignment = ifelse(cluster %in% needs_decision,
                             "MANUAL REQUIRED", paste("auto:", best)),
         reason = case_when(
           !(cluster %in% needs_decision)      ~ "",
           n_lineages > MAX_LINEAGES           ~ sprintf("%d lineages (%s), dominance %.1f",
                                                         n_lineages, lineages, dominance),
           is.na(margin)                       ~ "only one panel scorable",
           margin < AUTO_MARGIN & score < AUTO_MIN_SCORE
                                               ~ sprintf("margin %.3f and score %.3f both low",
                                                         margin, score),
           margin < AUTO_MARGIN                ~ sprintf("margin %.3f below %.2f",
                                                         margin, AUTO_MARGIN),
           TRUE                                ~ sprintf("score %.3f below %.2f",
                                                         score, AUTO_MIN_SCORE))) %>%
  select(cluster, n_cells, assignment, reason, best, score, margin,
         n_lineages, lineages, dominance, nFeature_median, nCount_median,
         percent_mt_median, percent_ribo_median, doublet_score_median)
save_tab(decision_table, "S1_01f_assignment_decisions")
print(decision_table %>% select(cluster, n_cells, assignment, reason), n = Inf)

## --- evidence 3: full read-out for each cluster that needs a decision ----------
# Everything needed to make the call, in one place: the ranked panel scores, the
# diagnostic genes, the sentinel values, and the genes that actually distinguish
# the cluster from the rest of the dataset.
if (length(needs_decision)) {
  flagged_markers <- map_dfr(needs_decision, function(k) {
    m <- tryCatch(
      FindMarkers(sc, ident.1 = k, only.pos = TRUE, min.pct = 0.25,
                  logfc.threshold = 0.5, verbose = FALSE),
      error = function(e) { warning("FindMarkers failed for cluster ", k,
                                    ": ", conditionMessage(e)); NULL })
    none <- tibble(cluster = k, rank = NA_integer_, gene = NA_character_,
                   avg_log2FC = NA_real_, pct.1 = NA_real_, pct.2 = NA_real_)
    if (is.null(m) || !nrow(m)) return(none)
    # p_val_adj can be NA; a bare comparison then yields NA and breaks if()
    keep <- which(!is.na(m$p_val_adj) & m$p_val_adj < 0.05)
    if (!length(keep)) return(none)
    m <- m[keep, , drop = FALSE]
    m <- m[order(-m$avg_log2FC), , drop = FALSE]
    i <- seq_len(min(15, nrow(m)))
    tibble(cluster = k, rank = i, gene = rownames(m)[i],
           avg_log2FC = round(m$avg_log2FC[i], 2),
           pct.1 = m$pct.1[i], pct.2 = m$pct.2[i])
  })
  save_tab(flagged_markers, "S1_01g_markers_of_flagged_clusters")

  for (k in needs_decision) {
    cat("\n", strrep("=", 78), "\n", sep = "")
    cat("CLUSTER ", k, "   n = ", qc_by_cluster$n_cells[qc_by_cluster$cluster == k],
        "\n", sep = "")
    cat(strrep("-", 78), "\n", sep = "")
    cat("reason: ", decision_table$reason[decision_table$cluster == k], "\n\n", sep = "")

    cat("panel scores (ranked)\n")
    v <- sort(panel_scores[k, ], decreasing = TRUE, na.last = NA)
    print(round(v[seq_len(min(5, length(v)))], 3))

    cat("\ndiagnostic genes (mean expression, non-zero only)\n")
    d <- diagnostics[, k]; names(d) <- rownames(diagnostics)
    print(sort(d[d > 0.05], decreasing = TRUE))

    cat("\nlineage sentinels above floor: ",
        multilineage$lineages[multilineage$cluster == k], "\n", sep = "")
    cat(multilineage$values[multilineage$cluster == k], "\n")

    cat("\nQC\n")
    print(as.data.frame(qc_by_cluster[qc_by_cluster$cluster == k, ]),
          row.names = FALSE)

    cat("\ntop markers vs all other cells\n")
    print(as.data.frame(flagged_markers[flagged_markers$cluster == k,
                                        c("gene", "avg_log2FC", "pct.1", "pct.2")]),
          row.names = FALSE)
  }
  cat("\n", strrep("=", 78), "\n\n", sep = "")
}

## --- fingerprint: detects that the manual block was written for another run ----
# Sizes of every cluster, in order. Any change upstream changes this string, and
# a CELLTYPE_MANUAL written against the old clustering is then rejected instead
# of being applied to differently numbered clusters.
cluster_fingerprint <- paste0(
  nrow(qc_by_cluster), "c:",
  paste(qc_by_cluster$n_cells[order(as.integer(qc_by_cluster$cluster))],
        collapse = ","))
message("Cluster fingerprint: ", cluster_fingerprint)

## >>> MANUAL CALLS ------------------------------------------------------------
## Fill these in from the read-out above, then re-run. Copy the fingerprint
## printed by the run you based the calls on.
##
## Guidance for the recurring cases in this dataset:
##   Low-quality  several mutually exclusive lineages elevated together (Alb with
##                Cdh5 with Rgs5 with Plp1). Must be excluded: these cells would
##                otherwise enter the fibroblast object.
##   Mural        Rgs5-high and Notch3-positive like a pericyte but Myh11-high
##                like smooth muscle, and Pdgfra-negative. The panel margin is
##                near zero because these markers genuinely do not resolve the
##                sub-identity; the descriptive label is the honest answer.
##   HSC          Lrat and Reln clearly above the fibroblast clusters. The margin
##                against Fibroblast is structurally small because both panels
##                share the collagen genes, so a small margin here is not doubt.
##
## Any label is allowed as long as it appears in CELLTYPE_ORDER below.

CELLTYPE_FINGERPRINT <-
  "16c:1065,716,463,458,402,363,305,278,276,219,184,149,92,64,42,29"
CELLTYPE_MANUAL <- c(
  # Lrat 4.04 and Reln 1.89, against Lrat 0.04-0.36 in the fibroblast clusters,
  # and Lrat is the sixth strongest marker of this cluster against all other
  # cells (log2FC 2.56, detected in 82% versus 15%). The margin against
  # Fibroblast is small because the two panels share the collagen genes, which is
  # a property of the panels rather than doubt about the call.
  "1" = "HSC",
  # Rgs5 21.77 and Notch3 4.22 like the pericytes, Myh11 43.16 like the vascular
  # smooth muscle, Cnn1 3.72 lower than either, Pdgfra 0.06. The panel margin of
  # 0.009 is the correct answer: these markers do not resolve the sub-identity,
  # so the label is descriptive. Its own markers are strong and specific
  # (Pgam2 10.19, Corin 7.43, Pln 5.06, all at pct.2 < 0.02), so this is a real
  # population, not a mixture.
  "5" = "Mural",
  # 598 detected genes and 966 UMIs against 2,400-4,200 and 8,100-17,500
  # elsewhere; Alb, Cdh5, Rgs5 and Plp1 all elevated together; and not one of its
  # fifteen strongest "markers" has pct.1 > pct.2, which is what a library too
  # small to normalise looks like rather than a cell type. Gm42418, a standard
  # low-quality indicator, is detected in 97.5% of these cells.
  "7" = "Low-quality")
## -----------------------------------------------------------------------------

if (length(needs_decision)) {
  if (is.na(CELLTYPE_FINGERPRINT) || !length(CELLTYPE_MANUAL))
    stop("\n", length(needs_decision), " cluster(s) need a manual call: ",
         paste(needs_decision, collapse = ", "),
         "\nRead the block printed above (also in S1_01f_assignment_decisions and ",
         "S1_01g_markers_of_flagged_clusters), then set CELLTYPE_MANUAL and ",
         "CELLTYPE_FINGERPRINT <- \"", cluster_fingerprint, "\" and re-run.")

  if (!identical(CELLTYPE_FINGERPRINT, cluster_fingerprint))
    stop("CELLTYPE_MANUAL was written for a different clustering.\n",
         "  recorded: ", CELLTYPE_FINGERPRINT, "\n",
         "  current:  ", cluster_fingerprint, "\n",
         "Cluster numbers are not comparable between the two. Re-read the ",
         "read-out above and rewrite the calls.")

  missing_calls <- setdiff(needs_decision, names(CELLTYPE_MANUAL))
  extra_calls   <- setdiff(names(CELLTYPE_MANUAL), needs_decision)
  if (length(missing_calls) || length(extra_calls))
    stop("CELLTYPE_MANUAL does not match the clusters that need a decision.\n",
         if (length(missing_calls))
           paste0("  no call for: ", paste(missing_calls, collapse = ", "), "\n"),
         if (length(extra_calls))
           paste0("  call given for cluster(s) that were assigned automatically: ",
                  paste(extra_calls, collapse = ", "),
                  " -- remove them, or lower AUTO_MARGIN if you mean to override ",
                  "an automatic call\n"))

  CELLTYPE[names(CELLTYPE_MANUAL)] <- unname(CELLTYPE_MANUAL)
  message("Manual calls applied: ",
          paste(sprintf("%s = %s", names(CELLTYPE_MANUAL), CELLTYPE_MANUAL),
                collapse = ", "))
}

if (anyNA(CELLTYPE))
  stop("unassigned cluster(s): ",
       paste(names(CELLTYPE)[is.na(CELLTYPE)], collapse = ", "))

# Several clusters usually map to the same lineage. Number them so that each
# cluster keeps its own label in the figures, while the lineage prefix stays
# greppable (for example "Fibroblast-1", "Fibroblast-2").
CELLTYPE <- ave(CELLTYPE, CELLTYPE, FUN = function(v)
  if (length(v) == 1) v else paste0(v, "-", seq_along(v)))
stopifnot(setequal(names(CELLTYPE), levels(Idents(sc))),
          !anyDuplicated(CELLTYPE))
print(CELLTYPE)
## -----------------------------------------------------------------------------

# Lineage order for every downstream figure. factor() would otherwise sort the
# labels alphabetically, which interleaves unrelated lineages (aEC next to the
# fibroblasts) and, worse, is locale-dependent: under a C locale the lower-case
# endothelial labels sort after the upper-case ones, so the same script would
# produce a different panel order on a different machine. Fixing the order here
# makes every dot plot, bar plot and legend reproducible.
CELLTYPE_ORDER <- c(
  "Fibroblast", "HSC",                       # collagen-producing mesenchyme
  "Pericyte", "Mural", "VSM",                # perivascular mesenchyme
  "sEC", "vEC", "aEC", "lyEC",               # endothelium
  "Mesothelium", "Hepatocyte",               # epithelium/mesothelium
  "Low-quality")

# Every label must come either from a panel name or from the manual overrides
# above. This catches a typo in CELLTYPE_ORDER at the point where it is written,
# rather than further down where the failure is harder to read.
stopifnot(all(CELLTYPE_ORDER %in%
              union(names(PANELS), c("Mural", "Low-quality"))))

# Expand the lineage order to the numbered labels actually present, keeping the
# -1/-2/-3 suffixes in cluster order within each lineage.
lineage_of  <- sub("-[0-9]+$", "", CELLTYPE)
stopifnot(all(lineage_of %in% CELLTYPE_ORDER))
celltype_levels <- unlist(lapply(CELLTYPE_ORDER,
                                 function(l) CELLTYPE[lineage_of == l]),
                          use.names = FALSE)

sc$celltype <- factor(unname(CELLTYPE[as.character(Idents(sc))]),
                      levels = celltype_levels)
stopifnot(!anyNA(sc$celltype))
Idents(sc)  <- sc$celltype
print(levels(sc$celltype))

save_fig(DimPlot(sc, label = TRUE, label.size = 3, repel = TRUE) +
           theme(legend.text = element_text(size = 8)),
         "explor_umap_allcells", 7, 4)

# Meflin across all populations, and the marker dot plot.
# Pdgfrb sits next to Pdgfra deliberately. It is used as a stromal marker in the
# histology and in situ hybridisation figures, and the Methods state that it was
# NOT used to define fibroblasts here because it is expressed more highly by
# pericytes. That justification has to be visible in a panel rather than asserted
# in the text, and the two receptors side by side make it self-evident: Pdgfra
# marks the fibroblast rows, Pdgfrb the mural ones. Omitting it would also leave
# a reader unable to connect this figure to the earlier ones.
s1b_genes <- c("Islr",
               "Col1a1", "Col1a2", "Pdgfra", "Pdgfrb", "Dcn", "Lum",
               "Thy1", "Pi16", "Mfap5", "Scara5",
               "Lrat", "Reln", "Lama1",
               "Rgs5", "Notch3", "Cspg4",
               "Acta2", "Tagln", "Cnn1", "Myh11", "Des", "Pln",
               "Pecam1", "Cdh5", "Stab2", "Clec4g", "Vwf", "Tek",
               "Sema3g", "Ly6c1", "Lyve1", "Prox1",
               "Msln", "Upk3b", "Lrrn4", "Krt19", "Krt7",
               "Alb", "Fabp1", "Apoa5", "Rgn",
               "Ptprc", "Cd3e")

# Same reading direction as the human panel in section 16: cell types follow the
# gene blocks, fibroblasts first, top to bottom. CELLTYPE_ORDER above already
# holds that order; the y axis is reversed because a discrete ggplot axis places
# the first level at the bottom.
save_fig(DotPlot(sc, features = present(s1b_genes, sc), dot.scale = 3,
                 cols = c("grey90", "red")) + RotatedAxis() +
           scale_y_discrete(limits = rev(levels(sc$celltype))) +
           theme(axis.title = element_blank(),
                 axis.text = element_text(size = 8),
                 legend.position = "top"),
         "explor_dotplot_allcells", 11, 3.5)

islr_by_celltype <- DotPlot(sc, features = "Islr")$data %>%
  as_tibble() %>%
  transmute(celltype = id, pct_expressing = round(pct.exp, 1),
            mean_expression = round(avg.exp, 3)) %>%
  arrange(desc(mean_expression))
save_tab(islr_by_celltype, "S1_02_Islr_by_celltype")
print(islr_by_celltype, n = Inf)

# Pdgfra and Pdgfrb adjacent, for the reason given above the dot plot.
save_fig(FeaturePlot(sc, features = present(c("Islr", "Thy1", "Col1a1",
                                              "Pdgfra", "Pdgfrb", "Acta2"), sc),
                     cols = c("grey90", "red"), order = TRUE, ncol = 3),
         "explor_featureplots_allcells", 9, 6)

saveRDS(sc, obj("02_allcells_annotated.rds"))


#@ STEP 7
# =============================================================================
# 7. REFERENCE EXPRESSION TABLES
# -----------------------------------------------------------------------------
# Two quantities the manuscript reports are ratios between fibroblasts and other
# stromal populations of the same tissue. They are computed here, on the full
# annotated object, because both require the excluded populations as reference.
# =============================================================================

reference_genes <- c("Col1a1", "Pdgfra", "Pdgfrb", "Islr", "Thy1",
                     "Rgs5", "Notch3", "Cspg4", "Kcnj8", "Abcc9",
                     "Acta2", "Tagln", "Cnn1", "Myh11", "Des",
                     "Lrat", "Reln", "Pecam1", "Cdh5", "Msln", "Krt19")

reference_table <- avg_expr(sc, reference_genes) %>%
  round(2) %>% as.data.frame()
save_tab(reference_table, "S1_03_reference_expression", row.names = TRUE)
print(reference_table)

## NOTE FOR INTERPRETATION
# Contractile genes in this table provide the evidence that the normal tissue
# contains no myofibroblasts: in the reported run, fibroblast clusters expressed
# Acta2, Cnn1 and Myh11 at 44- to 148-fold, 45- to 630-fold and 64- to 304-fold
# lower levels than the vascular smooth-muscle clusters. Tagln is the least
# discriminating of the four and should not be quoted alone.
#
# Pdgfrb, Kcnj8 and Abcc9 are included to document why they are not used to
# define fibroblasts: Pdgfrb is highest in pericytes, and Kcnj8 and Abcc9 are
# expressed at low absolute levels throughout.


## --- is a cholangiocyte population hidden inside the mesothelium? -------------
# The data published with GSE163777 annotate this cluster as biliary epithelium.
# Here it is annotated as mesothelium, on Msln, Upk3b and Lrrn4, and that is a
# disagreement with a published figure, so it has to be tested rather than
# asserted: a keratin-positive cluster could be mesothelium with a minority of
# cholangiocytes buried in it, and at whole-tissue resolution the two would
# merge. Mesothelium expresses Krt19 and Krt7, so the keratins cannot decide it;
# Epcam, Sox9, Cftr and Hnf1b can.
#
# RESULT OF THE REPORTED RUN (219 cells, resolution 0.5, three subclusters):
#
#              n   Msln    Epcam    Epcam+Msln-   median genes   median UMI
#   0         87   100%    1.1%          0             3634        11832
#   1         86   100%   14.0%          0             3430        11339
#   2         46    65%   15.2%          6 (13.0%)      1345         3566
#
# Three things follow.
#   - Subclusters 0 and 1 are unambiguous mesothelium: Msln in 100 % of cells
#     and NOT ONE Epcam+Msln- cell in 173 cells. The Epcam detected in
#     subcluster 1 is entirely in Msln+ cells and averages 0.13, against 2.58 in
#     subcluster 2 -- the same detection rate, a twentyfold lower level. That is
#     ambient, and detection rate alone would have missed the distinction.
#   - Subcluster 2 is low-complexity, not a different cell type: 1,345 genes
#     against ~3,500, every mesothelial marker down by about the same factor,
#     and percent.mt NOT raised (3.96 vs 3.61 and 3.20), so it is shallow
#     capture rather than dying cells. Dropout also explains its 21.7 %
#     Epcam-Msln- cells.
#   - Six cells (2.7 % of 219) are Epcam+Msln-. These are genuine epithelial
#     cells that escaped the EpCAM depletion used during isolation. Six cells
#     cannot be typed -- Cftr is LOWEST in subcluster 2 and Hnf1b is in three
#     cells -- so they are not called cholangiocytes here.
#
# The claim the manuscript can make is therefore "no cholangiocyte cluster was
# identified", not "no cholangiocytes". The annotation of the 219-cell cluster
# as mesothelium stands either way: 0 of the 173 high-quality cells are
# Epcam+Msln-. These cells are excluded with the mesothelium before the
# fibroblast analysis and touch no result.
meso_check <- function(x = sc, label = "Mesothelium", resolution = 0.5) {
  if (!label %in% levels(Idents(x))) {
    message("  no '", label, "' cluster; skipping the cholangiocyte check")
    return(invisible(NULL))
  }
  m <- subset(x, idents = label)
  m <- FindNeighbors(m, dims = 1:15, verbose = FALSE)
  m <- FindClusters(m, resolution = resolution, random.seed = CONFIG$seed,
                    verbose = FALSE)
  g <- present(c("Msln", "Upk3b", "Lrrn4", "Krt19", "Krt7",
                 "Epcam", "Sox9", "Cftr", "Hnf1b"), m)
  d <- FetchData(m, vars = g, layer = "data")
  ids <- Idents(m)
  out <- tibble(
    subcluster   = levels(ids),
    n            = as.integer(table(ids)),
    median_genes = round(tapply(m$nFeature_RNA, ids, median)),
    median_umi   = round(tapply(m$nCount_RNA,   ids, median)),
    pct_Msln     = round(tapply(d$Msln  > 0, ids, mean) * 100, 1),
    pct_Epcam    = round(tapply(d$Epcam > 0, ids, mean) * 100, 1),
    mean_Epcam   = round(tapply(expm1(d$Epcam), ids, mean), 2),
    # the number that decides it: epithelial cells that are not mesothelium
    n_Epcam_only = as.integer(tapply(d$Epcam > 0 & d$Msln == 0, ids, sum)))
  save_tab(out, "S1_03b_mesothelium_subclusters")
  print(as.data.frame(out))
  message("  Epcam+Msln- cells in total: ", sum(out$n_Epcam_only),
          " of ", ncol(m),
          ".  A cholangiocyte CLUSTER would show up as a subcluster that is ",
          "mostly Epcam+Msln-; scattered single cells are not one.")
  invisible(out)
}
meso_check()


#@ STEP 8
# =============================================================================
# 8. FIBROBLAST EXTRACTION
# =============================================================================

FIB_TYPES <- grep("^Fibroblast", levels(sc$celltype), value = TRUE)
stopifnot(length(FIB_TYPES) >= 1)

sc_fib <- subset(sc, idents = FIB_TYPES)
message("Fibroblasts: ", ncol(sc_fib), " cells")

# Purity, expressed as the ratio of mean expression in the retained fibroblasts
# to mean expression in the excluded source population. Detection rate is not
# used: ambient RNA inflates it and it is not interpretable here.
purity_genes <- c("Rgs5", "Notch3", "Cspg4", "Myh11", "Cnn1",
                  "Lrat", "Reln", "Pecam1", "Cdh5", "Ptprc", "Msln", "Krt19")

fib_mean <- rowMeans(avg_expr(sc_fib, purity_genes))
all_mean <- avg_expr(sc, purity_genes)
other    <- setdiff(colnames(all_mean), FIB_TYPES)

purity <- tibble(
  gene              = names(fib_mean),
  fibroblast_mean   = round(fib_mean, 3),
  max_other_mean    = round(apply(all_mean[names(fib_mean), other, drop = FALSE],
                                  1, max), 3),
  source_population = colnames(all_mean[, other, drop = FALSE])[
                        apply(all_mean[names(fib_mean), other, drop = FALSE],
                              1, which.max)]) %>%
  mutate(fold_lower = round(max_other_mean / pmax(fibroblast_mean, 1e-6), 1))
save_tab(purity, "S1_04_fibroblast_purity")
print(purity, n = Inf)

save_fig(DotPlot(sc_fib, features = present(c("Col1a1", "Pdgfra", "Islr",
                                              purity_genes), sc_fib),
                 cols = c("grey90", "red"), dot.scale = 4) + RotatedAxis(),
         "QC_04_fibroblast_purity", 9, 3)

# The extracted fibroblasts, before any re-clustering. This is the checkpoint
# section 9 starts from, and section 9 must not write to it: section 9 removes
# contaminated cells, so a section that both read and wrote this file would
# start each re-run from its own previous output and never give the same answer
# twice. Section 9 writes 03c_fibroblast_subsets.rds instead.
saveRDS(sc_fib, obj("03_fibroblasts.rds"))


#@ STEP 9
# =============================================================================
# 9. RE-CLUSTERING OF FIBROBLASTS
# =============================================================================

sc_fib <- FindVariableFeatures(sc_fib, selection.method = "vst",
                               nfeatures = CONFIG$n_hvg)
sc_fib <- ScaleData(sc_fib, features = rownames(sc_fib))
sc_fib <- RunPCA(sc_fib, features = VariableFeatures(sc_fib), npcs = 50)

save_fig(ElbowPlot(sc_fib, ndims = 50), "QC_05_elbow_fibroblast", 4, 3)

dims_fib <- seq_len(CONFIG$n_pcs_fib)
sc_fib   <- FindNeighbors(sc_fib, dims = dims_fib)

sc_fib <- FindClusters(sc_fib, resolution = CONFIG$resolution_fib_qc)
sc_fib$cluster_id <- factor(as.integer(as.character(Idents(sc_fib))) + 1L)
Idents(sc_fib)    <- sc_fib$cluster_id
message("Fibroblast clusters for the contamination check: ",
        nlevels(sc_fib$cluster_id))
print(table(sc_fib$cluster_id))

## --- second decision point: what is each fibroblast cluster? -------------------
# Two things can split a fibroblast compartment without being fibroblast biology,
# and both need to be seen before the subsets are named.
#
# 1. Ambient contamination that survived the first pass. A small cluster whose
#    markers belong to another lineage entirely (hepatocyte transcripts, say) is
#    carrying ambient RNA, not defining a fibroblast state. Merging it into a
#    subset silently contaminates that subset.
# 2. The dissociation response. Enzymatic digestion induces immediate-early and
#    chemokine genes in a subset of cells, and clustering will happily split a
#    homogeneous population along that axis (van den Brink et al., Nat Methods
#    2017). A split whose markers are Fos, Jun, Egr1, Cxcl1/2/3, Il6, Il11 is an
#    artefact of tissue handling and should be merged, not reported as a state.
#
# Both are reported here per cluster; neither is acted on automatically.

# The response to enzymatic digestion is broader than the AP-1 and heat-shock
# genes it is usually summarised by: it includes oxidative-stress, NF-kB and
# inflammatory-mediator genes, and a list restricted to the obvious members
# leaves most of it in place. A first version of this panel held 26 genes and
# removed only 19 of the 771 genes separating two clusters -- while the genes
# that remained at the top (Fosl1, Cyr61, Srxn1, Gch1, Txnip, Has1) were
# themselves stress-responsive. The set below is a superset informed by
# van den Brink et al., Nat Methods 2017.
#
# Genes that double as fibroblast identity markers are deliberately NOT included
# even though they appear in published dissociation lists: Dcn, Des, Tagln2 and
# Sox9 carry real information here, and removing them would hide biology rather
# than an artefact. This is the judgement the list encodes; change it knowingly.
DISSOCIATION_GENES <- c(
  # AP-1 and immediate early
  "Fos", "Fosb", "Fosl1", "Fosl2", "Jun", "Junb", "Jund", "Atf3", "Atf4",
  "Egr1", "Egr2", "Egr3", "Ier2", "Ier3", "Ier5", "Btg1", "Btg2", "Klf2",
  "Klf4", "Klf6", "Klf9", "Myc", "Nr4a1", "Nr4a2", "Nr4a3", "Per1", "Rhob",
  "Srf", "Trib1", "Csrnp1", "Gadd45b", "Gadd45g", "Ppp1r15a", "Phlda1",
  "Errfi1", "Skil", "Sertad1", "Midn", "Zfp36", "Zfp36l1", "Zfp36l2",
  # heat shock and chaperones
  "Hspa1a", "Hspa1b", "Hspa5", "Hspa8", "Hspb1", "Hsph1", "Hsp90aa1",
  "Hsp90ab1", "Dnaja1", "Dnajb1", "Dnajb4", "Bag3", "Hspe1",
  # NF-kB and inflammatory mediators
  "Nfkbia", "Nfkbiz", "Nfkbib", "Tnfaip3", "Tnfaip6", "Zc3h12a", "Sqstm1",
  "Socs3", "Stat3", "Irf1", "Litaf", "Sbno2", "Cebpb", "Cebpd",
  "Cxcl1", "Cxcl2", "Cxcl3", "Cxcl5", "Ccl2", "Ccl7", "Il6", "Il11", "Lif",
  "Ptgs2", "Timp1", "Serpine1", "Has1", "Sdc4", "Cd14",
  # oxidative stress and redox
  "Srxn1", "Txnip", "Gch1", "Hmox1", "Nqo1", "Gclm", "Mt1", "Mt2", "Sat1",
  # cytoskeletal and other reported members
  "Cyr61", "Ccn1", "Ankrd1", "Zyx", "Lmna", "Myadm", "Tiparp", "Ifrd1",
  "Slc38a2", "Slc2a3", "Pnrc1", "Odc1", "Mcl1", "Rrad", "Gem", "Pde4b",
  "Dusp1", "Dusp8", "Ubc", "Eif1", "Ccnl1")
diss <- present(DISSOCIATION_GENES, sc_fib)
sc_fib <- AddModuleScore(sc_fib, features = list(diss),
                         name = "dissociation_", seed = CONFIG$seed)

fib_cluster_qc <- tibble(
    cluster = as.character(sc_fib$cluster_id),
    nFeature = sc_fib$nFeature_RNA, nCount = sc_fib$nCount_RNA,
    percent_mt = sc_fib$percent.mt,
    dissociation = sc_fib$dissociation_1) %>%
  group_by(cluster) %>%
  summarise(n_cells = n(),
            nFeature_median = median(nFeature),
            nCount_median   = median(nCount),
            percent_mt_median = round(median(percent_mt), 2),
            dissociation_median = round(median(dissociation), 3),
            .groups = "drop")

# Lineage sentinels again, now within the fibroblast compartment: any cluster
# scoring on a non-fibroblast lineage is carrying ambient RNA from it.
fib_sentinel <- avg_expr(sc_fib, intersect(unique(unlist(LINEAGE_SENTINELS)),
                                           rownames(sc_fib)),
                         group.by = "cluster_id") %>% round(2)
fib_cluster_qc <- fib_cluster_qc %>%
  left_join(as.data.frame(t(fib_sentinel)) %>% rownames_to_column("cluster"),
            by = "cluster") %>%
  arrange(as.integer(cluster))
save_tab(fib_cluster_qc, "S1_05e_fibroblast_cluster_qc")
print(fib_cluster_qc, n = Inf, width = Inf)

# Markers of every fibroblast cluster against the rest, so that a cluster defined
# by another lineage or by the dissociation response is visible immediately.
fib_markers_first <- FindAllMarkers(sc_fib, only.pos = TRUE, min.pct = 0.25,
                                    logfc.threshold = 0.5, verbose = FALSE) %>%
  filter(!is.na(p_val_adj), p_val_adj < 0.05) %>%
  group_by(cluster) %>% slice_max(avg_log2FC, n = 15) %>% ungroup()
save_tab(fib_markers_first, "S1_05f_fibroblast_cluster_markers")

fib_summary <- fib_markers_first %>% group_by(cluster) %>%
  summarise(top15 = paste(gene, collapse = ", "), .groups = "drop") %>%
  left_join(fib_cluster_qc %>% select(cluster, n_cells, dissociation_median),
            by = "cluster")
print(fib_summary, n = Inf, width = Inf)

# For every pair of adjacent clusters, how much of the separation survives when
# the dissociation-response genes are set aside. A pair that separates on 800
# genes but on only a handful once those are removed is one population split by
# tissue handling; a pair that keeps a substantial, interpretable marker set is
# two states. This is the test that decides whether to merge.
pair_separability <- function(x, exclude, min_lfc = 0.5) {
  ids  <- Idents(x)
  cent <- t(sapply(levels(ids), function(k)
    colMeans(Embeddings(x, "pca")[ids == k, dims_fib, drop = FALSE])))
  dmat <- as.matrix(dist(cent)); diag(dmat) <- Inf
  pairs <- unique(t(apply(cbind(levels(ids),
                                colnames(dmat)[apply(dmat, 1, which.min)]),
                          1, sort)))
  map_dfr(seq_len(nrow(pairs)), function(i) {
    a <- pairs[i, 1]; b <- pairs[i, 2]
    m <- FindMarkers(x, ident.1 = a, ident.2 = b, logfc.threshold = min_lfc,
                     min.pct = 0.25, verbose = FALSE)
    m <- m[!is.na(m$p_val_adj) & m$p_val_adj < 0.05, , drop = FALSE]
    kept <- m[!rownames(m) %in% exclude, , drop = FALSE]
    kept <- kept[order(-abs(kept$avg_log2FC)), , drop = FALSE]
    top50 <- head(rownames(m)[order(-abs(m$avg_log2FC))], 50)
    tibble(cluster_a = a, cluster_b = b,
           n_genes = nrow(m),
           n_genes_excl_dissociation = nrow(kept),
           fraction_retained = round(nrow(kept) / max(nrow(m), 1), 3),
           # Counting genes is not enough: two clusters split by handling and two
           # split by identity can differ by a similar NUMBER of genes. What
           # differs is which genes. This is the fraction of the 50 strongest
           # discriminating genes that are stress-responsive.
           stress_fraction_top50 = round(mean(top50 %in% exclude), 3),
           top_remaining = paste(head(rownames(kept), 10), collapse = ", "))
  })
}

pair_sep <- pair_separability(sc_fib, exclude = diss)
save_tab(pair_sep, "S1_05g_pair_separability_without_dissociation")
print(pair_sep, n = Inf, width = Inf)

diss_hi <- fib_cluster_qc$cluster[
  fib_cluster_qc$dissociation_median >
    min(fib_cluster_qc$dissociation_median) + 0.15]
if (length(diss_hi))
  message("Cluster(s) with an elevated dissociation score: ",
          paste(diss_hi, collapse = ", "),
          " -- check S1_05f: if the markers are immediate-early and chemokine ",
          "genes, this split is a handling artefact and the clusters should be ",
          "merged rather than reported as separate states")

fib_fingerprint <- paste0(nlevels(sc_fib$cluster_id), "c:",
                          paste(as.integer(table(sc_fib$cluster_id)),
                                collapse = ","))
message("Fibroblast cluster fingerprint: ", fib_fingerprint)

## >>> FIBROBLAST CLUSTERS TO DROP ----------------------------------------------
## Clusters removed before the subsets are defined, as character cluster ids.
## character(0) keeps everything. Read S1_05e, S1_05f and S1_05g first.
##
## Cluster 5 is hepatocyte-contaminated: Alb 11.92 against 0.87-1.14 in the other
## clusters, all fifteen of its top markers are liver transcripts (Pigr, Ces3a,
## Kng1, Mug1, Apoh, Hamp, Fga, Ttr, Apoa1), and it is separated from its nearest
## neighbour by only 66 genes where every other pair differs by 771-895. Its
## complexity is normal (nFeature 3273), so these are fibroblasts carrying a
## heavy ambient load rather than damaged cells -- but their transcriptome is
## dominated by that load and cannot be used for differential expression.
FIBROBLAST_DROP     <- c("5")
FIBROBLAST_DROP_FOR <- "5c:717,449,317,301,49"
## -----------------------------------------------------------------------------

if (length(FIBROBLAST_DROP)) {
  # Same discipline as the cell-type calls: cluster numbers are not stable, so
  # the run these were chosen against is recorded and checked. Without this, a
  # change anywhere upstream would silently drop a different population.
  if (!identical(FIBROBLAST_DROP_FOR, fib_fingerprint))
    stop("FIBROBLAST_DROP was chosen against a different clustering.\n",
         "  recorded: ", FIBROBLAST_DROP_FOR, "\n",
         "  current:  ", fib_fingerprint, "\n",
         "Re-read S1_05e / S1_05f / S1_05g and choose again.")
  stopifnot(all(FIBROBLAST_DROP %in% levels(sc_fib$cluster_id)))
  message("Dropping fibroblast cluster(s) ", paste(FIBROBLAST_DROP, collapse = ", "),
          " (", sum(sc_fib$cluster_id %in% FIBROBLAST_DROP), " cells)")
  sc_fib <- subset(sc_fib, cells = colnames(sc_fib)[
    !sc_fib$cluster_id %in% FIBROBLAST_DROP])
  sc_fib$cluster_id <- droplevels(sc_fib$cluster_id)
  Idents(sc_fib) <- sc_fib$cluster_id
  message("Fibroblasts retained: ", ncol(sc_fib))
}


## --- second pass: embedding recomputed on the retained cells -------------------
# Dropping a contaminated cluster changes which genes are variable and therefore
# where the principal components point. Re-running the embedding on the retained
# cells means the subsets are defined on the cleaned data rather than on an
# embedding partly shaped by the cells that were removed.
sc_fib <- FindVariableFeatures(sc_fib, selection.method = "vst",
                               nfeatures = CONFIG$n_hvg)
sc_fib <- ScaleData(sc_fib, features = rownames(sc_fib))
sc_fib <- RunPCA(sc_fib, features = VariableFeatures(sc_fib), npcs = 50)
sc_fib <- FindNeighbors(sc_fib, dims = dims_fib)

res_scan <- map_dfr(seq(0.1, 1.2, by = 0.1), function(r)
  tibble(resolution = r,
         n_clusters = nlevels(FindClusters(sc_fib, resolution = r,
                                           verbose = FALSE)$seurat_clusters)))
save_tab(res_scan, "S1_05_resolution_scan")
print(res_scan, n = Inf)

## --- justifying the resolution -------------------------------------------------
# The number of clusters must not be chosen from how many blobs the UMAP appears
# to show. UMAP optimises local neighbourhoods; the size, separation and even the
# number of visible blobs depend on n_neighbors and min_dist, and a single
# continuous population is routinely drawn as two or three. Nor should the
# resolution be lowered until the expected number appears, which selects the
# answer rather than testing it.
#
# The test used here is separability: at a given resolution, can each cluster be
# distinguished from its most similar neighbour by genes? A partition in which
# some cluster has almost no genes separating it from its neighbour is splitting
# a continuum, and the lower resolution is correct. A partition in which every
# cluster has a solid marker set is real, and collapsing it hides a population.
separability <- function(x, resolution, min_lfc = 0.5, min_genes_call = 10) {
  y <- FindClusters(x, resolution = resolution, random.seed = CONFIG$seed,
                    verbose = FALSE)
  ids <- Idents(y)
  if (nlevels(ids) < 2) return(tibble())
  # nearest neighbour of each cluster in PCA space, by centroid distance
  cent <- t(sapply(levels(ids), function(k)
    colMeans(Embeddings(y, "pca")[ids == k, dims_fib, drop = FALSE])))
  dmat <- as.matrix(dist(cent)); diag(dmat) <- Inf
  map_dfr(levels(ids), function(k) {
    nb <- colnames(dmat)[which.min(dmat[k, ])]
    m  <- FindMarkers(y, ident.1 = k, ident.2 = nb, logfc.threshold = min_lfc,
                      min.pct = 0.25, verbose = FALSE)
    m  <- m[m$p_val_adj < 0.05, , drop = FALSE]
    tibble(resolution = resolution, cluster = k, n_cells = sum(ids == k),
           nearest_neighbour = nb, n_genes_vs_neighbour = nrow(m),
           top_genes = paste(head(rownames(m[order(-abs(m$avg_log2FC)), ]), 5),
                             collapse = ", "))
  })
}

# Compare the chosen resolution with the next step up. Edit this vector to the
# two or three values that bracket the decision in your own run.
sep <- map_dfr(CONFIG$resolution_candidates, separability, x = sc_fib)
save_tab(sep, "S1_05d_cluster_separability")
print(sep, n = Inf)

weak <- sep %>% filter(n_genes_vs_neighbour < 10)
if (nrow(weak))
  message("Clusters not separable from their neighbour (<10 genes): ",
          paste(sprintf("res %.1f cluster %s (%d genes)", weak$resolution,
                        weak$cluster, weak$n_genes_vs_neighbour),
                collapse = "; "),
          " -- these splits are not supported by expression")
## -----------------------------------------------------------------------------

sc_fib <- FindClusters(sc_fib, resolution = CONFIG$resolution_fib)
sc_fib <- RunUMAP(sc_fib, dims = dims_fib)

sc_fib$cluster_id <- factor(as.integer(as.character(Idents(sc_fib))) + 1L)
Idents(sc_fib)    <- sc_fib$cluster_id
message("Fibroblast subsets: ", nlevels(sc_fib$cluster_id),
        " at resolution ", CONFIG$resolution_fib)
print(table(sc_fib$cluster_id))

## --- is the subset structure compositional? -----------------------------------
# This is the check that matters. Ribosomal transcripts are the largest single
# block of genes separating the peribiliary from the periportal subset, and the
# separation survives UMI downsampling, so it is a genuine difference in library
# composition. The question is whether that difference DEFINED the clusters or
# merely accompanies them.
#
# Three outcomes:
#   ARI close to 1  the partition does not depend on these genes; report the
#                   comparison and keep the original clustering
#   ARI 0.7-0.9     the partition is mostly stable; check which cells moved
#   ARI below 0.7   the clustering was partly compositional; the version without
#                   ribosomal and mitochondrial variable features becomes the
#                   primary analysis and the rest of the script must be re-run
#                   from here on those labels
message(sprintf("Ribosomal/mitochondrial genes among the %d variable features: %d",
                length(VariableFeatures(sc_fib)),
                sum(grepl(CONFIG$exclude_pattern, VariableFeatures(sc_fib)))))
print(grep(CONFIG$exclude_pattern, VariableFeatures(sc_fib), value = TRUE))

pc_cor_fib <- pc_covariate_cor(sc_fib, dims_fib, "fibroblast")

alt_fib <- cluster_without_housekeeping(sc_fib, CONFIG$n_hvg, dims_fib,
                                        CONFIG$resolution_fib, CONFIG$seed)
ari_fib <- adj_rand(Idents(sc_fib), alt_fib)
message(sprintf("Fibroblast clustering without ribosomal/mitochondrial HVGs: %d clusters, ARI = %.3f",
                nlevels(alt_fib), ari_fib))

sc_fib$cluster_no_housekeeping <- alt_fib
crosstab_fib <- as.data.frame.matrix(
  table(original = sc_fib$cluster_id, no_housekeeping = alt_fib))
save_tab(crosstab_fib, "S1_05b_clustering_without_housekeeping", row.names = TRUE)
print(crosstab_fib)

# Ribosomal fraction per cluster. A cluster defined by composition shows a large
# separation here; clusters defined by identity do not have to.
ribo_by_cluster <- tibble(cluster = sc_fib$cluster_id,
                          percent_ribo = sc_fib$percent.ribo,
                          percent_mt   = sc_fib$percent.mt) %>%
  group_by(cluster) %>%
  summarise(n = n(),
            percent_ribo_median = round(median(percent_ribo), 2),
            percent_mt_median   = round(median(percent_mt), 2), .groups = "drop")
save_tab(ribo_by_cluster, "S1_05c_ribosomal_fraction_by_cluster")
print(ribo_by_cluster, n = Inf)

if (ari_fib < 0.7)
  warning("fibroblast clustering changes substantially when ribosomal and ",
          "mitochondrial genes are removed from the variable features (ARI = ",
          round(ari_fib, 3), "); the original partition is partly compositional")
## -----------------------------------------------------------------------------

# --- identification of the subsets -------------------------------------------
# PBF  peribiliary   Thy1-high, basement-membrane collagens, morphogen antagonists
# PPF  periportal    Thy1-low, Pi16/Mfap5/Scara5 adventitial signature
# CxF  Cxcl12-high   Cxcl12-high with low-level vascular-associated transcripts
#
# CxF is named descriptively. Its position in the tissue has not been determined
# in situ, so it is not called perivascular here.

SUBSET_PANELS <- list(
  PBF = c("Thy1", "Ptn", "Col4a5", "Col4a6", "Col7a1", "Notum", "Wif1",
          "Dkk2", "Hhip", "Fbn2", "Ntrk2"),
  PPF = c("Pi16", "Mfap5", "Scara5", "Clec3b", "Pcolce2", "Fbln2", "Efemp1",
          "Wt1", "Lrrn4cl"),
  CxF = c("Cxcl12", "Ednrb", "Gja4", "Nrxn1", "Notch3", "Fbln5", "Adamtsl2",
          "Colec11", "Tinagl1")
)

subset_scores <- score_panels(sc_fib, SUBSET_PANELS)
save_tab(as.data.frame(round(subset_scores, 3)),
         "Fig5_01_subset_panel_scores", row.names = TRUE)
print(round(subset_scores, 2))

subset_call <- setNames(
  colnames(subset_scores)[apply(subset_scores, 1, which.max)],
  rownames(subset_scores))
print(subset_call)

## >>> REVIEW ------------------------------------------------------------------
# Confirm the assignment against the defining markers below before proceeding.
# If two clusters are assigned the same label, either lower the resolution so
# that they merge, or give the additional cluster its own descriptive name.
SUBSET <- subset_call
stopifnot(all(c("PBF", "PPF") %in% SUBSET))
## -----------------------------------------------------------------------------

# CHANGED 2026-10-01: the level order was PBF, PPF, CxF and is now CxF, PPF,
# PBF. This single line fixes the reading order of every panel that groups by
# subset -- the numbered key on the Fig. 3A UMAP, the GO columns in Fig. 3D,
# the facets in Fig. 3E -- and of every table keyed on the subset. The main
# text lists the subsets in the same order.
#
# The flipped dot plot (Fig. 3B) is the one exception: dot_plot(flip = TRUE)
# used to reverse the identities itself, which is why that panel already read
# CxF, PPF, PBF while the UMAP key read the other way. That rev() has been
# removed in figures_final.R, so the panel now takes its order from here like
# everything else. If Fig. 3B comes out backwards, that rev() is what to look
# at -- do not fix it by reordering this line.
sc_fib$subset <- factor(unname(SUBSET[as.character(sc_fib$cluster_id)]),
                        levels = intersect(c("CxF", "PPF", "PBF"), SUBSET))
Idents(sc_fib) <- sc_fib$subset

# PBF, PPF and CxF stay as the internal codes: they key the output file names and
# every table written so far, and renaming them would break that trail. What the
# figures show is a separate question. The manuscript no longer asserts the
# anatomical correspondence in the cluster names -- the scRNA-seq describes
# phenotypes and the correspondence with the in situ compartments is made in the
# text -- so the figures carry phenotype labels.
#
# PPF was labelled "Pi16-high" until the numbers were checked, and that label was
# not supportable. Detection rate and mean normalised counts, by subset:
#
#              Islr   Thy1   Pi16  Cxcl12  Mfap5  Scara5      Islr Thy1 Pi16 Cxcl12 Mfap5 Scara5
#   PBF       94.2%  90.0%  49.8%  85.5%   8.7%   1.3%        2.95 4.51 2.87   6.27  0.18   0.01
#   PPF       91.2%  19.5%  36.0%  72.1%  55.2%  44.8%        7.05 0.36 4.28   4.56  7.14   1.44
#   CxF       95.1%  34.7%   1.2%  99.8%   1.6%   1.0%        3.56 0.46 0.03  26.54  0.02   0.01
#
# Pi16 is detected in HALF of PBF and only a third of PPF. It separates PBF and
# PPF together from CxF (about 40 % against 1.2 %) and does not separate PBF from
# PPF at all -- detection ranks PBF higher, mean counts rank PPF higher, by 1.4-
# and 1.5-fold respectively. Naming one of the two subsets after a gene that is
# no more expressed in it than in its neighbour misdescribes the data.
#
# Mfap5 does separate them: 55.2 % against 8.7 % and 1.6 %, and 7.14 against 0.18
# and 0.02 counts. It is expressed in the majority of the subset's cells, which
# Scara5 (44.8 %) is not, though Scara5 is the more specific of the two.
#
# Meflin is detected in over 91 % of every subset, so "Meflin-high" ALONE would
# read as "the Meflin-positive subset" and imply the other two are negative --
# the opposite of what this paper argues, and directly at odds with the Islr in
# situ hybridisation in Figure 1. Its level is nonetheless genuinely highest
# here: 7.05 against 2.95 and 3.56 counts, 2.0- to 2.4-fold.
#
# The compound label states both facts and cannot be misread as positivity:
# Mfap5 is the discriminator (55.2 % against 8.7 % and 1.6 %), Meflin carries
# "-high" because that is a claim about level within cells that are all
# positive.
#
# CHANGED 2026-10-01: Mfap5 carried a "+" until now, on the argument that a
# detection difference deserves a different mark from a level difference. The
# reader does not make that distinction -- what they see is three labels of
# which one is punctuated differently for no visible reason -- so all three
# discriminators now carry "-high" and only Meflin carries "+". The figure legend must say so:
#   "Meflin (Islr) was detected in over 91 % of cells in all three subsets;
#    -high denotes relative transcript abundance, not positivity."
#
# Note the asymmetry with the other two labels is deliberate. Thy1 and Cxcl12
# each separate their subset on their own; this one needs two genes because no
# single gene does both jobs.
# Written in a plain form that carries the superscripts as markup rules rather
# than as characters: "-high" and a trailing "+" are turned into superscripts by
# figures.R when it draws them, and stay readable here and in every CSV. Do not
# put Unicode superscripts in these strings -- they would end up in the tables,
# and the fonts used for the figures do not carry a superscript "high".
SUBSET_DISPLAY <- c(PBF = "Meflin+/Thy1-high",
                    PPF = "Meflin-high/Mfap5-high",
                    CxF = "Meflin+/Cxcl12-high")
stopifnot(all(levels(sc_fib$subset) %in% names(SUBSET_DISPLAY)))
sc_fib$subset_display <- factor(
  unname(SUBSET_DISPLAY[as.character(sc_fib$subset)]),
  levels = unname(SUBSET_DISPLAY[levels(sc_fib$subset)]))
message("Figure labels: ",
        paste(sprintf("%s -> %s", levels(sc_fib$subset),
                      levels(sc_fib$subset_display)), collapse = ", "))

# Every table written by this script keys its rows on the internal code, and
# every figure is labelled with the display name. That is deliberate -- renaming
# the codes would break the file names and the trail through the output -- but it
# leaves a reader of the supplementary tables with "PBF" in a column and
# "Meflin+/Thy1-high" in the figure beside it and no way to connect the two.
#
# The key is therefore written as its own table, generated from the same factor
# levels as the labels, so it cannot fall out of step with either. Cite it in the
# supplementary materials; it is two lines of legend and it removes the whole
# question.
subset_key <- tibble(code = levels(sc_fib$subset),
                     label = levels(sc_fib$subset_display),
                     n_cells = as.integer(table(sc_fib$subset)[
                       levels(sc_fib$subset)]))
save_tab(subset_key, "Fig5_00_subset_label_key")
print(subset_key)

defining <- avg_expr(sc_fib, c("Islr", "Thy1", "Pi16", "Mfap5", "Scara5",
                               "Cxcl12", "Ednrb", "Gja4", "Col1a1", "Pdgfra",
                               "Acta2", "Col4a5", "Notum"),
                     group.by = "subset") %>% round(2) %>% as.data.frame()
save_tab(defining, "Fig5_02_defining_markers", row.names = TRUE)
print(defining)

# The defining contrast must hold, or the assignment is wrong.
stopifnot(defining["Thy1", "PBF"] > defining["Thy1", "PPF"],
          defining["Pi16", "PPF"] > defining["Pi16", "PBF"])

save_fig(DimPlot(sc_fib, label = TRUE, label.size = 3) + NoLegend(),
         "explor_umap_fibroblasts", 3.2, 3)
save_fig(FeaturePlot(sc_fib, features = present(c("Islr", "Thy1", "Pi16",
                                                  "Cxcl12", "Col1a1", "Acta2"),
                                                sc_fib),
                     cols = c("grey90", "red"), order = TRUE, ncol = 3),
         "explor_featureplots_fibroblast", 9, 6)
save_fig(VlnPlot(sc_fib, features = present(c("Islr", "Thy1", "Pi16", "Cxcl12"),
                                            sc_fib),
                 pt.size = 0, ncol = 4) & theme(axis.title.x = element_blank()),
         "explor_violin_defining_genes", 10, 2.6)

# Sequencing depth by subset, reported because it bears on section 10.
depth <- sc_fib@meta.data %>% group_by(subset) %>%
  summarise(n = n(),
            nFeature_median = median(nFeature_RNA),
            nCount_median   = median(nCount_RNA),
            percent_mt_median = round(median(percent.mt), 2), .groups = "drop")
save_tab(depth, "Fig5_03_depth_by_subset")
print(depth)

# Fibroblasts after contamination removal, re-embedding and subset assignment.
# Deliberately a different file from the one this section reads (see section 8).
saveRDS(sc_fib, obj("03c_fibroblast_subsets.rds"))


#@ STEP 10
# =============================================================================
# 10. DEPTH-CONTROLLED MARKER GENES
# -----------------------------------------------------------------------------
# Marker genes are computed twice: at native depth, and after downsampling every
# cell to a common number of UMIs. Only genes significant in both are used for
# interpretation. Without this control, genes detected at trace levels appear
# enriched in whichever cluster was sequenced more deeply.
# =============================================================================

markers_native <- FindAllMarkers(sc_fib, only.pos = TRUE,
                                 min.pct = 0.25, logfc.threshold = 0.25)
save_tab(markers_native, "Fig5_04_markers_native")

group_median_umi <- tapply(sc_fib$nCount_RNA, sc_fib$subset, median)
target_umi <- if (is.na(CONFIG$downsample_umi))
  floor(min(group_median_umi)) else CONFIG$downsample_umi
message("Library size by subset (median): ",
        paste(names(group_median_umi), round(group_median_umi),
              sep = " = ", collapse = ", "))
message("Downsampling all cells to ", target_umi, " UMIs")
if (target_umi < 0.25 * median(sc_fib$nCount_RNA))
  warning("downsampling target is less than a quarter of the median library ",
          "size; detection will be lost and real markers may be discarded")

# ADDED 2026-10-01, two changes, both about reproducibility of the downsampling.
#
# 1. set.seed() immediately before SampleUMI. SampleUMI draws from the GLOBAL
#    RNG and ignores every Seurat seed argument, so without it the downsampled
#    matrix depends on how much randomness ran earlier in the session.
#    set.seed() at the top of the script is not enough -- any intervening
#    random call moves the stream. Section 4 already follows this rule for
#    scDblFinder; it was missing here.
#
# 2. The downsampled object is CACHED, exactly as the doublet calls are. Three
#    numbers in the Methods come through this matrix -- the depth-robust marker
#    counts (779 / 252 / 527), the post-downsampling median UMI (9,785-9,853)
#    and the whole GSEA, whose ranked lists are built from sc_ds -- and those
#    numbers are already written into the manuscript. The cache keeps them
#    attached to the object that produced them instead of regenerating a
#    different draw.
#    DO NOT DELETE 03b_fibroblasts_downsampled.rds. If it is deleted the
#    downsampling is redone; with set.seed() in place that is at least
#    reproducible from then on, but the three numbers above must be re-read
#    from the new run and the Methods updated.
#    Set CONFIG$redo_downsample <- TRUE to force a fresh draw on purpose.
ds_file <- obj("03b_fibroblasts_downsampled.rds")
if (file.exists(ds_file) && !isTRUE(CONFIG$redo_downsample)) {
  message("Reusing the cached downsampled object from ", ds_file)
  sc_ds <- readRDS(ds_file)
  stopifnot(identical(colnames(sc_ds), colnames(sc_fib)))
  Idents(sc_ds) <- sc_fib$subset
} else {
  sc_ds <- sc_fib
  set.seed(CONFIG$seed)
  ds_counts <- SampleUMI(get_counts(sc_fib), max.umi = target_umi, upsample = FALSE)
  dimnames(ds_counts) <- dimnames(get_counts(sc_fib))
  sc_ds <- tryCatch(SetAssayData(sc_ds, layer = "counts", new.data = ds_counts),
                    error = function(e)
                      SetAssayData(sc_ds, slot = "counts", new.data = ds_counts))
  sc_ds <- NormalizeData(sc_ds, verbose = FALSE)
  Idents(sc_ds) <- sc_fib$subset
}

markers_ds <- FindAllMarkers(sc_ds, only.pos = TRUE,
                             min.pct = 0.1, logfc.threshold = 0.25)
save_tab(markers_ds, "Fig5_05_markers_downsampled")

# sc_ds is used again in section 12 for the depth-controlled differential
# expression, so it is written to disk: re-running section 12 alone after
# restarting R then costs a readRDS rather than a second downsampling pass.
# Written only when it was just generated -- writing back an object that was
# read from this same file would be a no-op at best and, if anything upstream
# had modified it, would silently replace the cache the reported numbers came
# from.
if (!exists("ds_file") || !file.exists(ds_file) || isTRUE(CONFIG$redo_downsample))
  saveRDS(sc_ds, obj("03b_fibroblasts_downsampled.rds"))

markers_robust <- inner_join(
  markers_native %>% select(gene, cluster, avg_log2FC, pct.1, pct.2, p_val_adj),
  markers_ds     %>% select(gene, cluster, avg_log2FC, pct.1, pct.2, p_val_adj),
  by = c("gene", "cluster"), suffix = c("_native", "_downsampled")) %>%
  filter(p_val_adj_downsampled < 0.05, avg_log2FC_downsampled > 0.25) %>%
  arrange(cluster, desc(avg_log2FC_downsampled))
save_tab(markers_robust, "Fig5_06_markers_depth_robust")

message(sprintf("Markers: %d native, %d depth-robust (%.0f%% retained)",
                nrow(markers_native), nrow(markers_robust),
                100 * nrow(markers_robust) / nrow(markers_native)))

# How much detection the downsampling cost. If the median ratio is far below 1
# the control is destroying signal rather than removing a confound, and
# CONFIG$downsample_umi should be raised. A ratio near 1 means the two libraries
# were already comparable and the control was cheap.
detection_retained <- median(markers_robust$pct.1_downsampled /
                             markers_robust$pct.1_native)
message(sprintf("Detection retained after downsampling: %.2f (median pct.1 ratio)",
                detection_retained))
if (detection_retained < 0.5)
  warning("downsampling removed more than half of all detection events; ",
          "raise CONFIG$downsample_umi before interpreting the depth-robust lists")

top_robust <- markers_robust %>% group_by(cluster) %>%
  slice_max(avg_log2FC_downsampled, n = 10)
save_fig(DoHeatmap(sc_fib, features = top_robust$gene, size = 2.5) +
           theme(axis.text.y = element_text(size = 5)),
         "explor_marker_heatmap_top10", 7, 6)


#@ STEP 11
# =============================================================================
# 11. TRANSCRIPTIONAL PROGRAMME SCORES
# -----------------------------------------------------------------------------
# Programmes marked data_derived = TRUE were built from the marker genes of this
# dataset. They describe the clusters but cannot serve as independent validation
# of them, and are reported as descriptive only.
# =============================================================================

# WHY universal_Pi16 AND universal_Col15a1 ARE TWO SCORES AND NOT ONE
#
# Buechler et al., Nature 2021 did not describe one universal fibroblast. They
# described two, both Dpt+, present in essentially every tissue they surveyed
# and in both steady state and perturbation:
#
#   Pi16+     Dpt+ Pi16+ Ly6c1+ Sema3c+, adventitial in location
#   Col15a1+  Dpt+ Col15a1+ Hspg2+ Vwa1+ Cd34+, parenchymal, basement-membrane
#             associated
#
# The two are ordered, not interchangeable: their Dpt lineage tracing and
# trajectory analysis place Pi16+ upstream of Col15a1+ and both upstream of the
# tissue-specialised states. Collapsing them into a single "universal score"
# would throw away the only part of the framework that speaks to direction,
# which is precisely the part this manuscript wants to invoke. It would also
# turn two observations into one: the two panels share no gene, so a subset
# that scores highest on both has been shown to be universal-like twice
# independently, and AddModuleScore bins each panel against its own controls.
#
# Keep them separate. If universal_Col15a1 turns out to have a small
# median_range in Fig5_08_programme_scores, the honest report is that the
# Col15a1+ type is not resolved in this dataset -- not a merged score that
# hides which of the two axes carried the result.
PROGRAMMES <- list(
  universal_Pi16    = list(genes = c("Pi16", "Dpt", "Ly6c1", "Sema3c", "Ackr3"),
                           data_derived = FALSE),
  universal_Col15a1 = list(genes = c("Col15a1", "Hspg2", "Vwa1", "Cd34"),
                           data_derived = FALSE),
  # CHANGED 2026-09-29: Scara5 and Mfap5 removed. They are the two genes the
  # manuscript uses to NAME the periportal subset, so a panel containing them
  # cannot say anything the name has not already said.
  #
  # CHECKED 2026-09-29, and the result matters for how this panel is described.
  # All four remaining genes are depth-robust periportal markers:
  #   Fbln2    log2FC 4.24  p.adj 3.8e-132
  #   Efemp1   log2FC 2.77  p.adj 2.2e-46   (also a weak peribiliary marker, 0.58)
  #   Pcolce2  log2FC 2.29  p.adj 4.5e-34
  #   Clec3b   log2FC 2.06  p.adj 1.1e-91
  # So removing Scara5 and Mfap5 removed the naming genes but did NOT make the
  # score independent of the marker analysis, and an earlier comment here
  # claiming it did was wrong.
  #
  # The panel is kept anyway, because of what it is being asked. It is not
  # asked whether the periportal subset has a distinct programme -- the marker
  # analysis answers that, and this score would only restate it. It is asked
  # whether the cluster the manuscript CALLS periportal is anatomically
  # adventitial, which is an inference from a dissociated library and needs
  # independent support. Four adventitial genes from prior knowledge, none of
  # them used in the naming, all enriched in that cluster, is evidence for the
  # assignment. The figure legend must say that these genes are also among the
  # subset's markers, and must not present the P value as a test of a
  # programme.
  #
  # Redundancy to watch: in Buechler et al. the Pi16+ universal fibroblast IS
  # the adventitial population, so this panel and universal_Pi16 bear on the
  # same claim through non-overlapping gene sets. That is corroboration, but
  # only if the text says so once rather than counting it twice.
  #
  # CHANGED TO TRUE 2026-09-29. All four genes are depth-robust markers of the
  # Mfap5+ subset, so the score restates the marker analysis rather than testing
  # it, and the Results paragraph never used it: the claim the paragraph makes
  # is that the Mfap5+ subset is closest to the Pi16+ universal type, which
  # universal_Pi16 carries on its own. A panel in the figure that no sentence in
  # the text uses is a panel a reviewer will ask about, and the honest answer
  # would be that its genes are the subset's own markers. With the flag TRUE the
  # existing machinery moves it to explor_programme_scores_data_derived, where
  # it stays available if it is ever needed, and Fig. 3E becomes the two
  # externally defined universal programmes only.
  adventitial       = list(genes = c("Fbln2", "Pcolce2", "Clec3b", "Efemp1"),
                           data_derived = TRUE),
  # REMOVED 2026-09-29: myofibroblast, apCAF and iCAF. This dataset is not a
  # tumour and contains no CAFs, so a CAF classification cannot be tested in it;
  # scoring these panels here invited the reader to map three normal-tissue
  # fibroblast subsets onto myCAF/apCAF/iCAF, which is a claim the experiment
  # cannot make. The CAF question belongs to the tumour figures, not to Fig. 5.
  # Their definitions and the reasoning behind them are kept here so that the
  # decision is visible and reversible:
  #   myofibroblast  Acta2, Cthrc1, Postn, Tnc
  #                  Tagln and Thbs2 had already been removed from this panel:
  #                  both are depth-robust markers of the Cxcl12-high subset and
  #                  both are expressed by mural cells, so with them included the
  #                  score ranked the Cxcl12-high subset highest even though
  #                  Acta2 is almost threefold higher in the peribiliary subset
  #                  and Postn is a peribiliary marker -- the panel was reading
  #                  perivascular identity, not myofibroblastic differentiation.
  #   apCAF          H2-Ab1, H2-Aa, Cd74, Saa3
  #   iCAF           Il6, Cxcl1, Cxcl2, Ccl2, Lif, Has1, Has2  (added 2026-09-24
  #                  from Ohlund 2017 and Biffi 2019; note that Cxcl1 is also in
  #                  the chemotaxis leading edge of the Cxcl12-high subset, so a
  #                  high score there would have been partly shared genes rather
  #                  than independent agreement)
  # Four programmes remain, which is also what lets Fig. 5C be drawn 2 x 2.
  #
  # REMOVED 2026-09-24: mesothelial. It was there to rule out mesothelium for
  # the Mfap5+/Scara5+ subset, but the text no longer discusses mesothelium in
  # the fibroblast subsets, and a panel the text never mentions only prompts the
  # question it was meant to answer. The whole-tissue mesothelium call is a
  # separate matter and is evidenced in section 7 and Fig. S1.
  # PARTLY derived: Col4a5, Col4a6 and Col7a1 were noticed in the Thy1-high
  # marker list, the other five are canonical basement-membrane genes taken from
  # prior knowledge. This was left at FALSE on the argument that the majority of
  # the panel is independent of these data.
  #
  # CHANGED TO TRUE 2026-09-29. The argument does not survive the numbers.
  # Against Fig5_06_markers_depth_robust, the eight genes split like this:
  #   Col7a1  PBF 3.91   Col4a6  PBF 3.89   Col4a5  PBF 3.16   <- the derived three
  #   Nid1    PBF 0.41, PPF 0.80
  #   Col4a1  CxF 0.51   Col4a2  CxF 0.77
  #   Lama1, Lama5  not depth-robust markers of any subset
  # The three genes that came from the Thy1-high marker list are the only ones
  # that point at the Thy1-high subset, and they do so an order of magnitude
  # more strongly than anything else in the panel. The five independent genes
  # point at the Cxcl12-high subset or nowhere. So the score is not a panel that
  # happens to favour the peribiliary subset; it is three of its own markers
  # plus noise, and the majority-independent argument was counting genes rather
  # than weighting them.
  #
  # With the flag TRUE the existing machinery moves this out of Fig. 5C and into
  # the descriptive panel, which is what that machinery is for. The claim it was
  # supporting -- a position adjacent to epithelium -- is better made directly:
  # Col4a5, Col4a6 and Col7a1 are in the subset marker dot plot, where they are
  # presented as markers, which is what they are.
  basement_membrane = list(genes = c("Col4a1", "Col4a2", "Col4a5", "Col4a6",
                                     "Col7a1", "Lama1", "Lama5", "Nid1"),
                           data_derived = TRUE)
  # REMOVED 2026-09-23: morphogen_antagonist. The panel was written after Notum,
  # Wif1, Dkk2, Nkd2 and Hhip were seen among the top markers of the Thy1-high
  # subset, so scoring it and reporting that the Thy1-high subset ranks highest
  # (P = 7.8e-148) tested nothing. Figure S1F was dropped for this reason on
  # 2026-09-16; the programme is removed here so no score, table or panel is
  # regenerated. The genes remain in the subset marker dot plot, where they are
  # presented as markers rather than as a predicted programme.
  #
  # REMOVED 2026-09-24: vascular_associated. Four of its eight genes are
  # depth-robust markers of the Cxcl12-high subset, so that subset scored
  # highest by construction. It was kept as a descriptive supplementary panel;
  # with the panel dropped there is nothing for it to feed, and a circular score
  # in the record invites someone to quote it later. An externally published
  # perivascular signature would make the point testable; none is used here.
  #
  # As of 2026-09-29 one programme is data_derived: basement_membrane. Three
  # remain in Fig. 5C -- universal_Pi16, universal_Col15a1 and adventitial --
  # and the flag does exactly the job it was kept for, which is to move a panel
  # out of the main figure the moment the question is asked of it rather than
  # after a reviewer asks it.
)

prog_genes <- lapply(PROGRAMMES, function(p) present(p$genes, sc_fib))
prog_genes <- prog_genes[lengths(prog_genes) >= 3]

prog_derived <- map_lgl(PROGRAMMES[names(prog_genes)], "data_derived")

prog_info <- tibble(programme    = names(prog_genes),
                    n_genes      = lengths(prog_genes),
                    data_derived = prog_derived,
                    genes = map_chr(prog_genes, paste, collapse = ", "))
save_tab(prog_info, "Fig5_07_programme_definitions")
print(prog_info %>% select(-genes), n = Inf)

# A module score is a single number over several genes and can be driven by one
# or two of them. Write out the per-gene means so that every score in the figure
# can be checked against its constituents. This is the table that showed the
# myofibroblast panel was being carried by Tagln and Thbs2 rather than by Acta2.
prog_gene_expr <- avg_expr(sc_fib, unique(unlist(prog_genes)),
                           group.by = "subset") %>%
  round(3) %>% as.data.frame() %>% rownames_to_column("gene") %>%
  left_join(tibble(gene      = unlist(prog_genes),
                   programme = rep(names(prog_genes), lengths(prog_genes))),
            by = "gene") %>%
  relocate(programme) %>% arrange(programme, gene)
save_tab(prog_gene_expr, "Fig5_07b_programme_gene_expression")

# Drop any programme_* columns left by an earlier run. PROGRAMMES has changed
# length more than once (mesothelial and vascular_associated were removed,
# iCAF added), and AddModuleScore only overwrites the columns it writes, so a
# shortened list would otherwise leave a stale programme_N behind that
# figures.R would mislabel against the current definitions table.
stale <- grep("^programme_[0-9]+$", colnames(sc_fib@meta.data), value = TRUE)
if (length(stale)) sc_fib@meta.data[stale] <- NULL

sc_fib <- AddModuleScore(sc_fib, features = prog_genes, name = "programme_",
                         seed = CONFIG$seed)
score_cols <- paste0("programme_", seq_along(prog_genes))

scores_long <- sc_fib@meta.data[, score_cols, drop = FALSE] %>%
  setNames(names(prog_genes)) %>%
  mutate(subset = sc_fib$subset) %>%
  pivot_longer(-subset, names_to = "programme", values_to = "score")

score_summary <- scores_long %>% group_by(programme, subset) %>%
  summarise(mean = mean(score), median = median(score), .groups = "drop") %>%
  pivot_wider(names_from = subset, values_from = c(mean, median))

# Effect size, reported alongside the means. With this many cells the P values
# below are astronomically small for differences of any size, so the spread
# between the highest and lowest subset is the quantity that should be read.
score_summary <- scores_long %>% group_by(programme) %>%
  summarise(median_range = max(tapply(score, subset, median)) -
                           min(tapply(score, subset, median)),
            top_subset   = names(which.max(tapply(score, subset, median))),
            .groups = "drop") %>%
  left_join(score_summary, by = "programme") %>%
  arrange(desc(median_range))
save_tab(score_summary, "Fig5_08_programme_scores")
print(score_summary, n = Inf)

kruskal <- scores_long %>% group_by(programme) %>%
  summarise(p = kruskal.test(score ~ subset)$p.value, .groups = "drop") %>%
  mutate(p_adj = p.adjust(p, "BH"))
dunn <- scores_long %>% group_by(programme) %>%
  dunn_test(score ~ subset, p.adjust.method = "BH")
save_tab(kruskal, "Fig5_09_programme_kruskal")
save_tab(dunn,    "Fig5_10_programme_dunn")
print(kruskal, n = Inf)

## NOTE ON THE STATISTICS -------------------------------------------------------
# The unit of observation in these tests is the cell, and GSE163777 is a single
# 10x library (see section 3). There is therefore no biological replication in
# this dataset, and the P values quantify how many cells were sequenced rather
# than how reproducible the difference is between animals. Per-animal
# aggregation is not possible here. State this in the figure legend, report
# median_range as the primary quantity, and treat the tests as descriptive.
## -----------------------------------------------------------------------------

# The main panel shows only the programmes defined from external literature.
# Data-derived programmes are plotted separately so that a reader cannot mistake
# a circular score for independent evidence.
plot_scores <- function(df, name, ncol_wrap = 4, w = 9, h = 5)
  save_fig(
    ggplot(df, aes(subset, score, fill = subset)) +
      geom_boxplot(outlier.size = 0.2, linewidth = 0.3) +
      facet_wrap(~ programme, scales = "free_y", ncol = ncol_wrap) +
      theme_classic(base_size = 8) + NoLegend() +
      theme(axis.title.x = element_blank()),
    name, w, h)

derived_names <- names(prog_genes)[prog_derived]
plot_scores(scores_long %>% filter(!programme %in% derived_names),
            "explor_programme_scores", 4, 9, 5)
if (length(derived_names))
  plot_scores(scores_long %>% filter(programme %in% derived_names),
              "explor_programme_scores_data_derived",
              length(derived_names), 3 * length(derived_names), 3)

# Persist the scores (added 2026-09-29). Until now this section computed
# programme_1 ... programme_n into the in-session object and never wrote it
# back, so the checkpoint that figures.R loads -- 03c_fibroblast_subsets.rds,
# written at the end of section 9, before this section runs -- carried no
# programme_* columns and programme_long() returned NULL. Writing back to the
# same file is safe: this section only adds metadata columns, and the stale-
# column purge above keeps the column set in step with PROGRAMMES. Sections
# 12-15 read this file and are unaffected by the extra columns.
saveRDS(sc_fib, obj("03c_fibroblast_subsets.rds"))


#@ STEP 12
# =============================================================================
# 12. DIFFERENTIAL EXPRESSION, PBF vs PPF
# -----------------------------------------------------------------------------
# The peribiliary subset is sequenced about 30 % more deeply than the periportal
# subset (section 9). At native depth this shows up directly in the result: every
# ribosomal and mitochondrial transcript that reaches significance does so in the
# peribiliary direction and none in the periportal direction, which is not a
# biological pattern. Excluding those genes at the enrichment step hides the
# symptom without removing the cause, because any depth-sensitive gene is still
# in the list. The contrast is therefore repeated on the UMI-downsampled matrix
# built in section 10, and only genes significant in both are carried forward.
#
# Three tables are written:
#   _01  significant genes at native depth (the supplementary DEG table)
#   _01b every gene tested, unfiltered (provenance; not a DEG list)
#   _01c the depth-robust intersection, which is what feeds sections 13 and 14
# =============================================================================

Idents(sc_fib) <- sc_fib$subset

deg_raw <- FindMarkers(sc_fib, ident.1 = "PBF", ident.2 = "PPF",
                       min.pct = 0.25, logfc.threshold = 0.25) %>%
  rownames_to_column("gene") %>% as_tibble() %>% arrange(desc(avg_log2FC))
save_tab(deg_raw, "S1D_01b_DEG_PBF_vs_PPF_all_tested")

# The supplementary table holds significant genes only. Saving the full
# FindMarkers output under a name containing "DEG" invites the gene count to be
# quoted from the wrong table: here that would be 2502 rather than 765.
deg <- deg_raw %>% filter(p_val_adj < 0.05)
save_tab(deg, "S1D_01_DEG_PBF_vs_PPF")

deg_all <- FindMarkers(sc_fib, ident.1 = "PBF", ident.2 = "PPF",
                       min.pct = 0.1, logfc.threshold = 0) %>%
  rownames_to_column("gene") %>% as_tibble()
save_tab(deg_all, "S1D_02_DEG_PBF_vs_PPF_unfiltered")

# Depth control. sc_ds carries the same cells with every library downsampled to
# a common UMI count; Idents were set to subset in section 10.
deg_ds <- FindMarkers(sc_ds, ident.1 = "PBF", ident.2 = "PPF",
                      min.pct = 0.1, logfc.threshold = 0.25) %>%
  rownames_to_column("gene") %>% as_tibble()
save_tab(deg_ds, "S1D_01d_DEG_PBF_vs_PPF_downsampled")

deg_robust <- inner_join(
  deg     %>% select(gene, avg_log2FC, pct.1, pct.2, p_val_adj),
  deg_ds  %>% select(gene, avg_log2FC, pct.1, pct.2, p_val_adj),
  by = "gene", suffix = c("_native", "_downsampled")) %>%
  filter(p_val_adj_downsampled < 0.05,
         sign(avg_log2FC_native) == sign(avg_log2FC_downsampled)) %>%
  arrange(desc(avg_log2FC_downsampled))
save_tab(deg_robust, "S1D_01c_DEG_PBF_vs_PPF_depth_robust")

message(sprintf(
  "DEG PBF vs PPF: %d up in PBF, %d up in PPF at native depth; %d and %d depth-robust",
  sum(deg$avg_log2FC > 0), sum(deg$avg_log2FC < 0),
  sum(deg_robust$avg_log2FC_native > 0), sum(deg_robust$avg_log2FC_native < 0)))

# Direct measurement of the confound, for the Methods. A ratio far from 1:1 is
# evidence that depth, not biology, is driving part of the native-depth result.
housekeeping <- deg %>% filter(grepl(CONFIG$exclude_pattern, gene))
message(sprintf(
  "Ribosomal/mitochondrial transcripts among native-depth DEGs: %d up in PBF, %d up in PPF (%d depth-robust)",
  sum(housekeeping$avg_log2FC > 0), sum(housekeeping$avg_log2FC < 0),
  sum(grepl(CONFIG$exclude_pattern, deg_robust$gene))))
save_tab(housekeeping, "S1D_01e_housekeeping_genes_in_DEG")

save_fig(
  deg_all %>%
    mutate(lp = -log10(pmax(p_val_adj, 1e-300)),
           direction = case_when(p_val_adj < 0.05 & avg_log2FC >  0.25 ~ "PBF",
                                 p_val_adj < 0.05 & avg_log2FC < -0.25 ~ "PPF",
                                 TRUE ~ "ns")) %>%
    ggplot(aes(avg_log2FC, lp, colour = direction)) +
    geom_point(size = 0.4, alpha = 0.6) +
    scale_colour_manual(values = c(PBF = "#C1443C", PPF = "#3C6EC1",
                                   ns = "grey80")) +
    labs(x = "log2 fold change (PBF / PPF)", y = "-log10 adjusted P") +
    theme_classic(base_size = 8),
  "S1D_volcano", 4, 3.4)


#@ STEP 13
# =============================================================================
# 13. GENE ONTOLOGY OVER-REPRESENTATION ANALYSIS
# -----------------------------------------------------------------------------
# Run on the depth-robust gene lists from section 12, not on the native-depth
# lists. Ribosomal and mitochondrial transcripts are excluded from both the
# query and the universe; without that exclusion the peribiliary result is
# dominated by cytoplasmic-translation terms carried by the same 29 ribosomal
# genes, which is a consequence of sequencing depth rather than of phenotype.
# The raw enrichGO output returns only terms that already pass the cutoff, so
# its size is not a measure of anything; read the simplified table, in which
# redundant parent and child terms have been collapsed.
# =============================================================================

## --- inputs that are data frames, not checkpoints ------------------------------
# run.R restores the Seurat objects a step needs, because those are saved as
# .rds checkpoints. deg_robust (section 12) and markers_robust (section 10) are
# plain data frames that only ever existed in the session, so run_steps(13) in a
# fresh session failed on them with "object not found".
#
# They are already written as tables, so the fix is to read them back rather than
# to add two more checkpoint files. Re-run sections 10 and 12 instead if the
# clustering has changed: these tables are only as current as the run that wrote
# them, and nothing here can tell.
restore_tab <- function(name, file) {
  if (exists(name, envir = globalenv())) return(invisible(FALSE))
  f <- tab(paste0(file, ".csv"))
  if (!file.exists(f))
    stop(name, " is missing and so is its table:\n  ", f,
         "\nRun the section that produces it first.")
  message("  restoring ", name, " <- ", basename(f),
          " (written ", format(file.mtime(f), "%Y-%m-%d %H:%M"), ")")
  assign(name, tibble::as_tibble(utils::read.csv(f, stringsAsFactors = FALSE)),
         envir = globalenv())
  invisible(TRUE)
}
restore_tab("deg_robust",     "S1D_01c_DEG_PBF_vs_PPF_depth_robust")
restore_tab("markers_robust", "Fig5_06_markers_depth_robust")
restore_tab("deg_raw",        "S1D_01b_DEG_PBF_vs_PPF_all_tested")

drop_genes <- grepl(CONFIG$exclude_pattern, rownames(sc_fib))
universe <- rownames(sc_fib)[
  rowSums(get_counts(sc_fib) > 0) >= CONFIG$min_cells_gene & !drop_genes]

up_pbf <- deg_robust %>% filter(avg_log2FC_native >  0.25) %>% pull(gene) %>%
  setdiff(rownames(sc_fib)[drop_genes])
up_ppf <- deg_robust %>% filter(avg_log2FC_native < -0.25) %>% pull(gene) %>%
  setdiff(rownames(sc_fib)[drop_genes])
message(sprintf("GO input: %d genes up in PBF, %d up in PPF (depth-robust)",
                length(up_pbf), length(up_ppf)))

run_go <- function(genes, label) {
  if (length(genes) < 10) {
    warning("too few genes for ", label, " (", length(genes), "); skipped")
    return(NULL)
  }
  res <- enrichGO(genes, universe = universe, OrgDb = org.Mm.eg.db,
                  keyType = "SYMBOL", ont = "BP", pAdjustMethod = "BH",
                  pvalueCutoff = 0.05, qvalueCutoff = 0.05)
  if (is.null(res) || nrow(as.data.frame(res)) == 0) {
    warning("no enriched terms for ", label)
    return(NULL)
  }
  save_tab(as.data.frame(res), paste0("S1D_03_GO_", label))
  simp <- clusterProfiler::simplify(res, cutoff = 0.7, by = "p.adjust",
                                    select_fun = min)
  save_tab(as.data.frame(simp), paste0("S1D_04_GO_", label, "_simplified"))
  message(sprintf("GO %s: %d terms, %d after simplify()",
                  label, nrow(as.data.frame(res)), nrow(as.data.frame(simp))))
  list(full = res, simplified = simp)
}

go_pbf_res <- run_go(up_pbf, "up_in_PBF")
go_ppf_res <- run_go(up_ppf, "up_in_PPF")
go_pbf <- go_pbf_res$simplified
go_ppf <- go_ppf_res$simplified

## --- every subset, not only the two being contrasted ---------------------------
# The PBF-versus-PPF contrast above answers one question and leaves the largest
# population of the three uncharacterised: the Cxcl12-high subset is 1,165 of
# 1,784 cells and never appears in it. That is a hole in the Results as written,
# because the manuscript groups the two Thy1-low subsets together as candidates
# for a periportal identity, and a reader cannot judge that grouping without
# seeing what the second of them expresses.
#
# One-versus-rest rather than three pairwise contrasts: it gives each subset a
# description on the same footing, and it is one testing family per subset
# instead of three overlapping ones.
#
# The input is the DEPTH-ROBUST marker list from section 10, which already
# exists, so this costs one enrichGO call per subset and no new differential
# expression. Using markers_native here instead would reintroduce the depth
# asymmetry that section 10 was written to control.
# markers_robust is already filtered on the downsampled statistics; requiring
# native significance as well is what makes it depth-robust in both directions.
# levels() would return NULL if FindAllMarkers had given a character column, so
# the subsets are taken from the object's own factor and kept in its order.
go_subsets <- list()
for (cl in levels(sc_fib$subset)) {
  g <- markers_robust %>%
    filter(as.character(cluster) == cl, p_val_adj_native < 0.05) %>%
    pull(gene) %>% setdiff(rownames(sc_fib)[drop_genes])
  message(sprintf("GO input for %s (one vs rest): %d depth-robust genes",
                  cl, length(g)))
  go_subsets[[cl]] <- run_go(g, paste0("up_in_", cl, "_vs_rest"))
}

# Saved so that the figures can be redrawn without repeating the enrichment,
# which is the slowest step in this section. figures.R reads this file.
saveRDS(list(pbf = go_pbf_res, ppf = go_ppf_res, subsets = go_subsets),
        obj("05_go_results.rds"))

# REMOVED 2026-09-23: the Wnt/Hedgehog antagonist analysis.
#
# This block searched the peribiliary enrichment result for Wnt, hedgehog, BMP
# and smoothened terms, and then reported Notum, Wif1, Dkk2, Nkd1/2, Hhip and
# Kremen1 at the gene level as a fallback. It was written as a test of a
# "pre-specified" hypothesis, but the hypothesis was not pre-specified: the gene
# panel was assembled on 2026-09-14 from the top markers of the Thy1-high subset
# itself, so showing that those genes are restricted to that subset tests
# nothing. Figure S1F was dropped on 2026-09-16 for the same reason and the code
# is removed here so that a later run does not regenerate the figure and the two
# CSVs (S1D_04b, S1D_04c) that no longer belong to any panel.
#
# The genes themselves remain legitimate as markers and are shown that way in
# the subset marker dot plot. What is gone is the claim that they constitute an
# independently predicted programme. Testing that claim would need an externally
# published Wnt-antagonist signature, or a measurement in tissue; neither is
# done here.

go_panel <- function(x, title)
  if (is.null(x)) NULL else
    dotplot(x, showCategory = 12, x = "GeneRatio") +
      ggtitle(title) + theme(text = element_text(size = 8))

go_panels <- Filter(Negate(is.null),
                    list(go_panel(go_pbf, "Upregulated in PBFs"),
                         go_panel(go_ppf, "Upregulated in PPFs")))
if (length(go_panels))
  # Exploratory only. The pairwise peribiliary-versus-periportal contrast is no
  # longer a manuscript panel -- Figure S1D now carries the one-versus-rest
  # enrichment, which describes all three subsets on the same footing. Renamed
  # from FigS1D_GO on 2026-09-23 so that the exploratory folder cannot be
  # mistaken for the figure set.
  save_fig(Reduce(`|`, go_panels), "explor_GO_pairwise_PBF_vs_PPF",
           4.5 * length(go_panels), 4)


#@ STEP 14
# =============================================================================
# 14. GENE SET ENRICHMENT ANALYSIS
# -----------------------------------------------------------------------------
# Threshold-free, and therefore complementary to section 13 rather than a
# repetition of it. Ribosomal and mitochondrial transcripts are removed from the
# ranked list for the reason given in the header.
#
# REWRITTEN 2026-09-29. Until now this section ran ONE contrast, PBF versus PPF,
# on the native-depth ranked list. Two things were wrong with that once the
# figure changed:
#
#   1. It is not the contrast in the figure. Section 13 became a one-versus-rest
#      enrichment over the three subsets on 2026-09-23 and that is what the
#      manuscript shows. The Cxcl12-high subset -- 1,165 of 1,784 cells, 65 % of
#      the fibroblasts -- never entered the GSEA at all. A stated negative about
#      a contrast that is no longer shown is not the negative a reviewer asking
#      "did you do a threshold-free analysis of your three subsets?" wants.
#   2. The one significant set in the old run, HALLMARK_MYC_TARGETS_V1, had a
#      leading edge of initiation and elongation factors, chaperonins,
#      spliceosome and proteasome components. CONFIG$exclude_pattern removes the
#      ribosomal PROTEINS and nothing else, and the subsets differ in ribosomal
#      content (17.0 / 13.5 / 12.8 %), so that hit was the depth confounder
#      under a biological heading. Excluding gene families after the fact only
#      renames it.
#
# Both are fixed the same way: the ranked lists are built ONE VERSUS REST, one
# per subset, FROM THE DOWNSAMPLED OBJECT (sc_ds, section 10), where every
# library carries the same number of UMIs. That removes the depth asymmetry at
# source instead of by exclusion, and it puts the GSEA on the same footing as
# the enrichment in the figure.
#
# PRE-SPECIFIED, BEFORE LOOKING. Three subsets x four collections (GO-BP,
# hallmark, Reactome, WikiPathways) is TWELVE Benjamini-Hochberg families. Every
# one of them is reported, significant or not, and the count of families is
# reported with them. Picking the collection that gave something is what this
# design exists to prevent.
#
# WHAT IS NOT DONE HERE, deliberately: the native-depth ranked list is not run
# alongside as a sensitivity check. The downsampled matrix IS the control;
# running both and reporting the intersection would double the families to
# twenty-four for no gain, since a set significant only at native depth is
# exactly the thing the downsampling was built to discard.
# =============================================================================

# Printed first so that "is this the current file?" is answered by the run
# itself. Bump it with every change to this section.
SEC14_REVISION <- "2026-09-29b (one-vs-rest, effect-size ranking)"
message("Section 14 revision: ", SEC14_REVISION)

# THE RANKING STATISTIC, and why it is avg_log2FC and not a signed P value.
#
# The first version of this section ranked genes by
# sign(avg_log2FC) * -log10(p_val), matching the pairwise block below. That was
# a mistake, and it is worth recording exactly how it failed, because the result
# looked like a finding.
#
# It returned 193 significant sets across the twelve families, and the NET
# DIRECTION of each subset's result was ordered precisely by that subset's
# number of detected genes:
#
#   subset   nFeature (downsampled)   sets up   sets down
#   CxF      2929                     43        7
#   PBF      2833                     4         6
#   PPF      2716                     7         126
#
# The mechanism: at equal UMI depth the three subsets still differ in what their
# transcriptomes are made of -- ribosomal transcripts are 17.0 / 13.5 / 12.8 %
# of counts (peribiliary / periportal / Cxcl12-high, unchanged by downsampling,
# which equalises depth and not composition) and the periportal subset detects
# the fewest genes (2,716 against 2,833 and 2,929). Whatever is displaced is
# displaced by a few percent, which is an avg_log2FC of a few hundredths and is
# nothing. But a shift of that size applied CONSISTENTLY to thousands of
# housekeeping genes, tested on 308 cells against 1,476, produces P values
# around 1e-20. The signed -log10 P statistic then places every one of those
# genes at the bottom of the ranking. Eighty Reactome sets -- proteasome,
# ubiquitin, cell cycle, NF-kB, stress response -- came out depleted in the
# periportal subset for that reason; re-ranking by effect size left one.
#
# Do not compute these percentages with PercentageFeatureSet() on sc_ds. It
# divides by the stored nCount_RNA column, which SetAssayData left at its native
# value, so the numerator is downsampled and the denominator is not. That gave
# 12.8 / 12.1 / 10.0 % on 2026-09-29 and the error was briefly used to argue
# that HALLMARK_MYC_TARGETS could not be compositional. With the correct
# figures it plainly can: weighting each subset's complement by cell number, the
# peribiliary subset sits at 17.0 % against 13.0 % for the rest, the periportal
# subset at 13.5 % against 13.7 %, and the Cxcl12-high subset at 12.8 % against
# 15.3 % -- which predicts the enrichment observed in the first and its absence
# in the other two. manuscript_numbers.R computes them from the matrix.
#
# Ranking by avg_log2FC removes it. The same data, re-ranked, drops the
# periportal depletion from 39 GO-BP terms to 21, destroys the ordering by
# nFeature, and eliminates the proteasome and cell-cycle block entirely. What
# survives is two axes: a type I interferon / antiviral response high in the
# Cxcl12-high subset (defense response to virus NES +2.25, response to
# interferon-beta +2.02), and extracellular matrix organization high in the
# Mfap5+ subset (NES +2.09, adjusted P 3.9e-5), which is also what the
# over-representation analysis in section 13 finds.
#
# This is the single-library problem in a new place. The script already says
# elsewhere that with cells as the unit of observation a P value measures how
# many cells were sequenced; here that P value had been built into the ranking
# itself, which is the one place it does the most damage.
#
# The pct filter is the second half of the fix: avg_log2FC is unstable for genes
# detected in a handful of cells, and ranking BY it promotes exactly those genes
# unless they are excluded. 0.1 in either group is the usual floor.
rank_one_vs_rest <- function(obj, cl) {
  Idents(obj) <- obj$subset
  FindMarkers(obj, ident.1 = cl, ident.2 = NULL,
              logfc.threshold = 0, min.pct = 0, verbose = FALSE) %>%
    rownames_to_column("gene") %>% as_tibble() %>%
    filter(!grepl(CONFIG$exclude_pattern, gene)) %>%
    filter(pmax(pct.1, pct.2) >= 0.1) %>%
    filter(is.finite(avg_log2FC)) %>%
    arrange(desc(avg_log2FC))
}

gsea_file <- obj("06b_gsea_one_vs_rest_all.rds")
go_objects <- list()
if (file.exists(gsea_file) && !isTRUE(CONFIG$redo_gsea)) {
  gsea_all <- readRDS(gsea_file)
  message(sprintf("GSEA: reusing the cached full result (%d rows). ",
                  nrow(gsea_all)),
          "Set CONFIG$redo_gsea <- TRUE to recompute.")
} else {
gsea_all <- map_dfr(levels(sc_fib$subset), function(cl) {
  message("GSEA, ", cl, " versus rest (downsampled)")
  r <- rank_one_vs_rest(sc_ds, cl)
  save_tab(r, paste0("Fig5_11_GSEA_ranked_list_", cl))
  gl <- setNames(r$avg_log2FC, r$gene)
  message("  ", length(gl), " genes ranked by avg_log2FC")

  # eps = 0 because fgsea otherwise floors the P value at 1e-10 and warns that
  # the true values are smaller. It costs time and buys an honest number.
  #
  # pvalueCutoff = 1 keeps every tested set, significant or not. clusterProfiler
  # applies Benjamini-Hochberg across all tested sets BEFORE this cutoff is
  # used, so no adjusted P value changes; the cutoff only decides how many rows
  # come back. seed = TRUE is what makes fgsea's permutations reproducible, and
  # it was already set -- without it this change would have re-rolled every P
  # value in the manuscript.
  go <- gseGO(gl, OrgDb = org.Mm.eg.db, keyType = "SYMBOL", ont = "BP",
              minGSSize = 15, maxGSSize = 500, pAdjustMethod = "BH",
              pvalueCutoff = 1, seed = TRUE, verbose = FALSE, eps = 0)
  go_df <- as.data.frame(go)
  message(sprintf("  GO-BP: %d sets tested, %d significant",
                  nrow(go_df), sum(go_df$p.adjust < 0.05)))
  if (nrow(go_df)) {
    go_df$collection <- "GO_BP"; go_df$subset <- cl
    go_df$n_sets_in_collection <- nrow(go_df)
  }

  msig <- if (is.null(CONFIG$msigdb_dir)) tibble() else
    map_dfr(names(CONFIG$msigdb_collections), function(nm) {
      gmt <- clusterProfiler::read.gmt(
        file.path(CONFIG$msigdb_dir, CONFIG$msigdb_collections[[nm]]))
      res <- GSEA(gl, TERM2GENE = gmt, minGSSize = 15, maxGSSize = 500,
                  pAdjustMethod = "BH", pvalueCutoff = 1, seed = TRUE,
                  verbose = FALSE, eps = 0)
      df <- as.data.frame(res)
      message(sprintf("  %-13s %5d sets in the GMT, %d tested, %d significant",
                      nm, length(unique(gmt$term)), nrow(df),
                      sum(df$p.adjust < 0.05)))
      if (nrow(df)) {
        df$collection <- nm
        df$subset <- cl
        # The denominator that matters is the number of sets actually TESTED
        # (those passing minGSSize/maxGSSize against this ranked list), not the
        # number in the GMT file. The previous version recorded the latter.
        df$n_sets_in_collection <- nrow(df)
      }
      df
    })
  # The gseaResult OBJECT, not just its data frame. A running-score enrichment
  # curve (enrichplot::gseaplot2) needs the object, and a table cannot be
  # turned back into one. Kept per subset for GO BP only, which is where the
  # sets under discussion live; the MSigDB objects would double the file for
  # no current use.
  go_objects[[cl]] <<- go
  bind_rows(if (nrow(go_df)) go_df else tibble(), msig)
})
  if (nrow(gsea_all)) saveRDS(gsea_all, gsea_file)
  if (length(go_objects)) saveRDS(go_objects, obj("06c_gsea_go_objects.rds"))
}

if (nrow(gsea_all)) {
  gsea_all <- gsea_all %>% as_tibble() %>%
    arrange(subset, collection, desc(abs(NES)))
  save_tab(gsea_all, "Fig5_12c_GSEA_one_vs_rest_all")
}

# The significant table is a FILTER of the full result, not a second run.
gsea_subsets <- if (nrow(gsea_all))
  dplyr::filter(gsea_all, p.adjust < 0.05) else gsea_all

# Guard against a silent change to a number the manuscript quotes. If the
# previously saved table is still on disk, compare the significant-set count.
.prev <- tab("Fig5_12_GSEA_one_vs_rest.csv")
if (file.exists(.prev)) {
  .n <- nrow(utils::read.csv(.prev, stringsAsFactors = FALSE))
  if (.n != nrow(gsea_subsets))
    warning(sprintf(paste("GSEA significant sets changed from %d to %d.",
                          "Every GSEA number in the manuscript must be",
                          "re-checked before this result is used."),
                    .n, nrow(gsea_subsets)))
}
rm(.prev)

n_families <- length(levels(sc_fib$subset)) *
  (1L + length(CONFIG$msigdb_collections))
message(sprintf("GSEA one-versus-rest: %d testing families, %d significant sets in total",
                n_families, nrow(gsea_subsets)))
if (nrow(gsea_subsets)) {
  # as_tibble() is not cosmetic. map_dfr over as.data.frame(gseaResult) returns
  # a base data.frame, and print.data.frame passes n = Inf into print.default,
  # where "n" partially matches na.print and the call fails with
  # "invalid 'na.print' specification". That error aborted the run on
  # 2026-09-29 AFTER the table had been written but BEFORE the plot.
  gsea_subsets <- gsea_subsets %>% as_tibble() %>%
    arrange(subset, collection, desc(abs(NES)))
  save_tab(gsea_subsets, "Fig5_12_GSEA_one_vs_rest")

  # The plot is written before anything is printed, so a formatting problem in
  # the console output can never again cost the figure.
  top_sets <- gsea_subsets %>%
    mutate(label = paste0(Description, "  [", collection, "]")) %>%
    group_by(subset) %>% slice_max(abs(NES), n = 10) %>% ungroup()
  save_fig(
    ggplot(top_sets, aes(NES, reorder(label, NES), colour = p.adjust,
                         size = setSize)) +
      geom_point() + facet_wrap(~ subset, scales = "free_y") +
      geom_vline(xintercept = 0, linewidth = 0.3) +
      scale_colour_gradient(low = "#C1443C", high = "#3C6EC1") +
      labs(x = "Normalised enrichment score", y = NULL) +
      theme_classic(base_size = 8),
    "explor_GSEA_one_vs_rest", 12, 6)

  # BASE R ON PURPOSE. This section runs after library(clusterProfiler) and
  # library(org.Mm.eg.db), which mask dplyr verbs; count() in particular
  # resolves to something that fails here with "Argument 'x' is not a vector:
  # list". Two runs were lost to that and to print(data.frame, n = Inf) before
  # it. Console summaries use table() and base printing, which nothing masks.
  # The saved tables still use dplyr, and those have always worked.
  message("Significant sets by family:")
  print(as.data.frame(table(subset     = gsea_subsets$subset,
                            collection = gsea_subsets$collection)))

  # THE GUARD THAT CAUGHT THE RANKING ARTEFACT, kept so it cannot come back
  # unnoticed. If the net direction (up minus down) of each subset's result is
  # ordered by that subset's number of detected genes, the analysis is
  # measuring detection, not biology. Print both and let the reader compare.
  cds_  <- get_counts(sc_ds)
  nfeat <- tapply(Matrix::colSums(cds_ > 0), sc_ds$subset, median)
  dtab  <- table(gsea_subsets$subset, ifelse(gsea_subsets$NES > 0, "up", "down"))
  up    <- if ("up"   %in% colnames(dtab)) as.integer(dtab[, "up"])   else 0L
  down  <- if ("down" %in% colnames(dtab)) as.integer(dtab[, "down"]) else 0L
  guard <- data.frame(subset = rownames(dtab), up = up, down = down,
                      net = up - down,
                      nFeature = as.integer(nfeat[rownames(dtab)]))
  message("Net direction against detected-gene count -- if these rank together, ",
          "suspect the ranking, not the biology:")
  print(guard[order(-guard$nFeature), ], row.names = FALSE)

  # THE FIRST THING TO LOOK AT, not the top hits. Two different confounders
  # have produced apparently biological results in this analysis already:
  # HALLMARK_MYC_TARGETS_V1 from ribosomal composition, and an eighty-set
  # proteasome and cell-cycle block from the ranking statistic. Both were in
  # the translation/proteostasis family. If anything from that family is still
  # here, read it before reading any hit.
  depth_like <- gsea_subsets %>%
    filter(grepl(paste0("ribosom|translation|MYC|spliceosom|proteasom|rRNA|",
                        "peptide chain|nonsense-mediated|oxidative phosphor"),
                 Description, ignore.case = TRUE))
  if (nrow(depth_like)) {
    message("Depth-confounder candidates still present -- read these first:")
    print(as.data.frame(depth_like)[, c("subset", "collection", "Description",
                                        "NES", "p.adjust")], row.names = FALSE)
  } else {
    message("No translation/ribosome/MYC set is significant; the downsampling ",
            "removed the confounder that the pairwise run reported.")
  }

  message("Top 5 by |NES| in each family:")
  d_ <- as.data.frame(gsea_subsets)
  top5 <- do.call(rbind, lapply(
    split(d_, list(d_$subset, d_$collection), drop = TRUE),
    function(x) head(x[order(-abs(x$NES)), ], 5)))
  print(top5[, c("subset", "collection", "Description", "NES", "p.adjust",
                 "setSize")], row.names = FALSE)
} else {
  # This is a result, not a failure. Write the empty table so that the run
  # leaves a record of the analysis having been done.
  save_tab(tibble(subset = character(), collection = character(),
                  Description = character(), NES = numeric(),
                  p.adjust = numeric()),
           "Fig5_12_GSEA_one_vs_rest")
  message("No gene set reached significance in any of the ",
          n_families, " families. Report this in the Methods as stated.")
}


## --- THE ORIGINAL PAIRWISE GSEA, kept and switched off -------------------------
# PBF versus PPF, native depth. Superseded by the block above on 2026-09-29 and
# left runnable because the pairwise DEG table is still a supplementary table,
# so someone may want its threshold-free counterpart.
#
# It is OFF by default for one reason only: it adds four more testing families
# to the twelve above, for a contrast the figure no longer shows. If it is
# switched on, the Methods must report sixteen families, not twelve.
#
# For the record, what it gave on 2026-09-23: GO-BP 20 terms, adjusted P 0.013
# to 0.044, three of them photoreceptor terms carried by the same five genes
# (Thy1, Ntrk2, Ptn, Cfh, Vegfa) and four of them 234-482-gene terms; hallmark
# 1 of 50 (HALLMARK_MYC_TARGETS_V1, the depth confounder); Reactome 0;
# WikiPathways 0.
RUN_PAIRWISE_GSEA <- FALSE

if (RUN_PAIRWISE_GSEA) {

ranks <- deg_all %>%
  filter(!grepl(CONFIG$exclude_pattern, gene)) %>%
  mutate(stat = sign(avg_log2FC) * -log10(pmax(p_val, 1e-300))) %>%
  filter(is.finite(stat)) %>%
  arrange(desc(stat))
gene_list <- setNames(ranks$stat, ranks$gene)
save_tab(ranks, "S1D_05_GSEA_ranked_list")

gsea_go <- gseGO(gene_list, OrgDb = org.Mm.eg.db, keyType = "SYMBOL",
                 ont = "BP", minGSSize = 15, maxGSSize = 500,
                 pAdjustMethod = "BH", pvalueCutoff = 0.05, seed = TRUE)
save_tab(as.data.frame(gsea_go), "S1D_06_GSEA_GO_BP")
saveRDS(gsea_go, obj("06_gsea_go.rds"))

gsea_df <- as.data.frame(gsea_go)
if (nrow(gsea_df) > 0) {
  top <- gsea_df %>% mutate(direction = ifelse(NES > 0, "PBF", "PPF")) %>%
    group_by(direction) %>% slice_max(abs(NES), n = 12) %>% ungroup()
  save_fig(
    ggplot(top, aes(NES, reorder(Description, NES), colour = p.adjust,
                    size = setSize)) +
      geom_point() +
      scale_colour_gradient(low = "#C1443C", high = "#3C6EC1") +
      geom_vline(xintercept = 0, linewidth = 0.3) +
      labs(x = "Normalised enrichment score  (positive, PBF)", y = NULL) +
      theme_classic(base_size = 8),
    # Exploratory only. Renamed from FigS1E_GSEA_GO on 2026-09-23: the GSEA is
    # reported as a negative result in the Methods and has no manuscript panel,
    # so it must not hold a figure letter. The analysis itself is kept -- a
    # stated negative forestalls a reviewer asking for exactly this.
    "explor_GSEA_GO_BP", 7, 5)
}

# --- MSigDB collections, from local GMT files ---------------------------------
# One analysis per pre-specified collection. Results are written per collection
# and also combined, with the collection recorded, so that a hit can never be
# quoted without saying which family it came from and how large that family was.

if (is.null(CONFIG$msigdb_dir)) {
  message("CONFIG$msigdb_dir not set; MSigDB enrichment skipped.")
  gsea_msig_all <- NULL
} else {
  stopifnot(dir.exists(CONFIG$msigdb_dir))
  gmt_files <- file.path(CONFIG$msigdb_dir, CONFIG$msigdb_collections)
  names(gmt_files) <- names(CONFIG$msigdb_collections)
  absent <- gmt_files[!file.exists(gmt_files)]
  if (length(absent))
    stop("GMT file(s) not found:\n  ", paste(absent, collapse = "\n  "),
         "\nCheck CONFIG$msigdb_dir and CONFIG$msigdb_collections.")
  if (any(grepl("entrez", gmt_files)))
    stop("entrez GMT selected; the ranked list is keyed on symbols and would ",
         "match nothing. Use the *.symbols.gmt files.")

  gsea_msig_all <- map_dfr(names(gmt_files), function(nm) {
    gmt <- clusterProfiler::read.gmt(gmt_files[[nm]])
    res <- GSEA(gene_list, TERM2GENE = gmt, minGSSize = 15, maxGSSize = 500,
                pAdjustMethod = "BH", pvalueCutoff = 0.05, seed = TRUE,
                verbose = FALSE)
    df  <- as.data.frame(res)
    message(sprintf("MSigDB %s: %d sets tested, %d significant",
                    nm, length(unique(gmt$term)), nrow(df)))
    save_tab(df, paste0("S1D_07_GSEA_MSigDB_", nm))
    if (nrow(df)) {
      save_fig(dotplot(res, showCategory = 20, split = ".sign") +
                 facet_grid(. ~ .sign) + theme(text = element_text(size = 8)),
               # Exploratory only; renamed from FigS1F_ on 2026-09-23 (see above).
               paste0("explor_GSEA_MSigDB_", nm), 8, 5)
      df$collection  <- nm
      df$n_sets_in_collection <- length(unique(gmt$term))
      df
    } else tibble()
  })
  if (nrow(gsea_msig_all))
    save_tab(gsea_msig_all %>% arrange(collection, desc(abs(NES))),
             "S1D_07_GSEA_MSigDB_all_collections")
}

}  # end if (RUN_PAIRWISE_GSEA)


# NOTE. Sections 15 (cross-species signature scoring) and 16 (human annotation
# audit) of the working script analysed a human dCCA single-cell dataset that is
# not part of the published study; they are omitted here. Section numbering is
# left unchanged so that it matches the working copy.

#@ STEP 17
# =============================================================================
# 17. PROVENANCE
# =============================================================================

saveRDS(sc_fib, obj("04_fibroblasts_final.rds"))

parameters <- tibble(
  parameter = c("seed", "min genes per cell", "max percent mitochondrial",
                "max genes per cell", "highly variable genes",
                "principal components, all cells", "resolution, all cells",
                "cells raw", "cells after QC", "cells after doublet removal",
                "principal components, fibroblasts", "resolution, fibroblasts",
                "fibroblasts analysed", "fibroblast clusters",
                "downsampling target (UMI)", "genes excluded from enrichment",
                "pseudotime"),
  value = c(CONFIG$seed, CONFIG$min_features, CONFIG$max_percent_mt,
            ifelse(is.na(CONFIG$max_features), "none", CONFIG$max_features),
            CONFIG$n_hvg, CONFIG$n_pcs_all, CONFIG$resolution_all,
            n_raw, n_qc, n_singlet,
            CONFIG$n_pcs_fib, CONFIG$resolution_fib,
            ncol(sc_fib), nlevels(sc_fib$cluster_id),
            target_umi, CONFIG$exclude_pattern, "not performed"))
save_tab(parameters, "00_analysis_parameters")
print(parameters, n = Inf)

writeLines(capture.output(sessionInfo()), out("sessionInfo.txt"))
message("\nFiles written under ", normalizePath(OUT_DIR))
print(list.files(OUT_DIR, recursive = TRUE))


# =============================================================================
# METHODS TEXT
# =============================================================================
# Mouse bilio-vascular tree single-cell RNA sequencing data (GSE163777) were
# analysed using Seurat (v4.3.0) in R (v4.2.0). Cells with fewer than 200
# detected genes or more than 10% mitochondrial transcripts were excluded, and
# doublets were removed using scDblFinder (v1.12.0). Data were log-normalised
# (LogNormalize; scale factor, 10,000) and scaled, 2,000 highly variable genes
# were identified using variance-stabilising transformation, and principal
# component analysis was performed on these genes; the first 30 principal
# components were used in downstream analyses. Cells were clustered using a
# shared-nearest-neighbour graph (resolution = 0.5) and visualised using uniform
# manifold approximation and projection. Populations were annotated by scoring
# each cluster against panels of canonical markers, using mean expression per
# cluster rather than detection rate, because ambient RNA inflates detection
# rate for abundant transcripts. Fibroblasts were defined by high expression of
# Col1a1 and Pdgfra together with absence of the mural markers Rgs5 and Notch3,
# which distinguished them from the pericyte and vascular smooth-muscle clusters
# of the same tissue; Pdgfrb was not used for this purpose because it was
# expressed more highly by pericytes than by fibroblasts. Lrat-positive hepatic
# stellate cells were excluded. The remaining fibroblasts were extracted in
# silico and re-clustered, and subsets were identified using Meflin (Islr), Thy1,
# Pi16 and Cxcl12. Because transcriptome complexity differed between clusters,
# marker genes were re-tested after downsampling all cells to a common
# sequencing depth, and only genes significant in both analyses were used for
# interpretation. Published fibroblast and cancer-associated fibroblast
# programmes were scored per cell using AddModuleScore and compared between
# subsets using the Kruskal-Wallis test followed by Dunn's post hoc test with
# Benjamini-Hochberg adjustment. Because this dataset comprises a single library
# and the unit of observation in these tests is therefore the cell rather than
# the animal, the tests are reported as descriptive and the difference in median
# score between subsets is given as the effect size. Peribiliary and periportal
# fibroblasts were compared using FindMarkers at native sequencing depth and
# again after downsampling to a common depth, and only genes significant in both
# comparisons and concordant in direction were taken forward to enrichment
# analysis. Ribosomal protein and mitochondrially encoded transcripts were
# excluded before enrichment analysis, because their contribution reflects
# differences in transcriptome complexity between clusters rather than a
# biological programme. Gene Ontology over-representation analysis of the
# depth-robust differentially expressed genes and gene set enrichment analysis of
# the full ranked gene list were performed using clusterProfiler (v4.6.0) and
# org.Mm.eg.db (v3.16.0), with redundant Gene Ontology terms collapsed using the
# simplify function. Pseudotime analysis was not performed: the normal
# extrahepatic biliary stroma contained no myofibroblast population, so no
# terminal state was available in this dataset, and myofibroblastic conversion
# was instead assessed experimentally by lineage tracing after bile duct
# ligation and in orthotopic tumours.
