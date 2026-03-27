# ============================================================================
# Figure 5B Supplementary: IntNMF Euclidean-Based Silhouette + DOCX Report
# ============================================================================
# Distance metric : Euclidean distance on the HVG VST matrix (dat_rna)
# Clustering      : IntNMF fit$clusters (K = 2)
# Output          : DOCX summary saved to the Rscripts directory
#
# ⚠ PREREQUISITE: Run Figure5.patients_bulkRNA_analysis.R first so that the
#   following objects exist in your R session:
#     fit, meta, cpi_100, cpi_300, dat_rna, K_final
# ============================================================================

suppressPackageStartupMessages({
  library(cluster)   # silhouette()
  library(officer)   # DOCX generation
  library(dplyr)
  library(tidyr)
})

# Output directory (save DOCX alongside the R scripts)
SCRIPT_DIR <- "/Users/daehwankim/Desktop/KIST_folder/2025/Progress/Cachexia model/Serpina3_Paper/Nature communications/Rscripts_natcom"

# ----------------------------------------------------------------------------
# 1. Euclidean Distance  +  Silhouette
# ----------------------------------------------------------------------------
cat("Computing Euclidean distance matrix...\n")
dist_euc <- dist(dat_rna, method = "euclidean")

cat("Computing silhouette widths (Euclidean)...\n")
sil_euc     <- silhouette(fit$clusters, dist_euc)
sil_sum     <- summary(sil_euc)

cluster_sizes <- sil_sum$clus.sizes      # named: 1 -> n_C1, 2 -> n_C2
cluster_avg   <- sil_sum$clus.avg.widths # per-cluster average sil width
mean_sil      <- sil_sum$avg.width
indiv_stats   <- sil_sum$si.summary     # Min / 1st Qu. / Median / Mean / 3rd Qu. / Max.

# Map IntNMF numeric clusters to C1 / C2
meta$cluster_hc <- paste0("C", fit$clusters)

# Cancer type × cluster crosstab
breakdown      <- as.data.frame(table(Cluster = meta$cluster_hc, CancerType = meta$group))
breakdown_wide <- pivot_wider(breakdown, names_from = CancerType, values_from = Freq)
c1_pan  <- breakdown_wide$Pancreas[breakdown_wide$Cluster    == "C1"]
c1_col  <- breakdown_wide$Colorectal[breakdown_wide$Cluster  == "C1"]
c2_pan  <- breakdown_wide$Pancreas[breakdown_wide$Cluster    == "C2"]
c2_col  <- breakdown_wide$Colorectal[breakdown_wide$Cluster  == "C2"]

# Print to console
cat("\n--- Silhouette Results (Euclidean) ---\n")
cat("print(mean_cpi_100):\n"); print(rowMeans(cpi_100))
cat("print(mean_cpi_300):\n"); print(rowMeans(cpi_300))
cat("\n")
print(sil_euc)
cat("\nIndividual silhouette widths:\n")
print(indiv_stats)
cat(sprintf("------------------------------------------------\nMean silhouette width =  %.7f\n", mean_sil))
cat("\nCluster × Cancer Type:\n"); print(breakdown_wide)

# ----------------------------------------------------------------------------
# 2. Build DOCX
# ----------------------------------------------------------------------------
cat("\nGenerating DOCX report...\n")

add_heading <- function(doc, text, level = 1) {
  style <- if (level == 1) "heading 1" else "heading 2"
  body_add_par(doc, text, style = style)
}
add_para <- function(doc, text) body_add_par(doc, text, style = "Normal")
add_blank <- function(doc)      body_add_par(doc, "",   style = "Normal")

doc <- read_docx()

# ── Title
doc <- doc %>%
  add_heading("IntNMF Euclidean-Based Clustering Summary", level = 1) %>%
  add_para(paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))) %>%
  add_blank()

# ── Section 1: Silhouette Output
doc <- doc %>%
  add_heading("1. Silhouette Results", level = 2) %>%
  add_para("CPI (100 iterations):") %>%
  add_para(paste(sprintf("k%d", 2:5), collapse = "        ")) %>%
  add_para(paste(round(rowMeans(cpi_100), 7), collapse = " ")) %>%
  add_para("CPI (300 iterations):") %>%
  add_para(paste(sprintf("k%d", 2:5), collapse = "        ")) %>%
  add_para(paste(round(rowMeans(cpi_300), 7), collapse = " ")) %>%
  add_blank() %>%
  add_para(sprintf(
    "Silhouette of %d units in %d clusters from silhouette.default(x = fit$clusters, dmatrix = distance) :",
    nrow(meta), K_final
  )) %>%
  add_para("Cluster sizes and average silhouette widths:") %>%
  add_para(paste(cluster_sizes, collapse = "        ")) %>%
  add_para(paste(round(cluster_avg, 7), collapse = " ")) %>%
  add_blank() %>%
  add_para("Individual silhouette widths:") %>%
  add_para(sprintf(
    "Min. 1st Qu.  Median    Mean 3rd Qu.    Max.\n%.4f   %.4f   %.4f   %.4f   %.4f   %.4f",
    indiv_stats["Min."],   indiv_stats["1st Qu."],
    indiv_stats["Median"], indiv_stats["Mean"],
    indiv_stats["3rd Qu."], indiv_stats["Max."]
  )) %>%
  add_para(strrep("-", 48)) %>%
  add_para(sprintf("Mean silhouette width =  %.7f", mean_sil)) %>%
  add_blank()

# ── Section 2: Results
c1_n   <- cluster_sizes["1"]
c2_n   <- cluster_sizes["2"]
c1_sil <- round(cluster_avg[1], 4)
c2_sil <- round(cluster_avg[2], 4)

doc <- doc %>%
  add_heading("2. Results (English)", level = 2) %>%
  add_para(sprintf(paste0(
    "Unsupervised clustering using an integrative non-negative matrix factorization ",
    "(IntNMF) framework identified two sample groups (K = 2) across the %d skeletal ",
    "muscle RNA-seq samples. Cluster coherence was assessed with silhouette width on the ",
    "Euclidean distance matrix of the HVG-by-sample VST expression matrix. ",
    "The mean silhouette width was %.4f ",
    "(cluster C1 [n = %d]: %.4f; cluster C2 [n = %d]: %.4f; ",
    "minimum individual silhouette = %.4f), indicating moderate-to-good separation between ",
    "the two groups. Both cancer types were represented in each cluster ",
    "(C1: Pancreas %d, Colorectal %d; C2: Pancreas %d, Colorectal %d), ",
    "suggesting that the clustering structure was not driven solely by cancer type."
  ),
    nrow(meta), mean_sil,
    c1_n, c1_sil, c2_n, c2_sil,
    indiv_stats["Min."],
    c1_pan, c1_col, c2_pan, c2_col
  )) %>%
  add_blank()

# ── Section 3: Methods
doc <- doc %>%
  add_heading("3. Methods (English)", level = 2) %>%
  add_para(paste0(
    "Gene-level counts were imported from kallisto quantifications using tximport and ",
    "analyzed in DESeq2. Genes were pre-filtered by retaining those with at least 10 counts ",
    "in at least 10 samples, and variance-stabilizing transformation (VST; blind = TRUE) was ",
    "applied to obtain a normalized expression matrix. Highly variable genes (HVGs) were ",
    "selected as the top 2,000 genes by variance across samples in VST space."
  )) %>%
  add_blank() %>%
  add_para(paste0(
    "For IntNMF, the HVG matrix was arranged as a samples-by-genes matrix and transformed ",
    "to satisfy the non-negativity requirement by shifting values to be \u2265 0 and rescaling ",
    "by the global maximum. IntNMF clustering was performed with K = 2 using the multiplicative ",
    "update/NNLS-based solver (nmf.mnnals; n.ini = 30 initializations, maxiter = 300, ",
    "st.count = 20, seed = TRUE). Cluster assignments were taken from the fitted model output ",
    "and labeled as C1 and C2."
  )) %>%
  add_blank() %>%
  add_para(paste0(
    "Cluster coherence was evaluated using silhouette width computed with Euclidean distances ",
    "on the HVG VST expression matrix, as implemented in the cluster R package. ",
    "The optimal number of clusters (K = 2) was selected based on the Cophenetic Prediction ",
    "Index (CPI) computed across K = 2\u20135 with 5-fold cross-validation and 30 runs per ",
    "fold at both 100 and 300 iterations."
  ))

# ── Save
docx_path <- file.path(
  SCRIPT_DIR,
  "Figure5B_IntNMF_euclidean_clustering_summary.docx"
)
print(doc, target = docx_path)
cat("\n-> DOCX saved to:", docx_path, "\n")
cat("============================================\n")
cat("Done!\n")
cat("============================================\n")
