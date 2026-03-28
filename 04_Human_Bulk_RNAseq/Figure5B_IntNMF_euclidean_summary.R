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

