# ============================================================================
# Figure 3: Proteomics Meta-Analysis for Cancer Cachexia
# ============================================================================
# This script performs:
#   - Figure 3A  : Meta-analysis volcano plot (Weighted Stouffer)
#   - Figure 3A' : Y-axis capped volcano plot (prettier version)
#   - Figure 3B  : Sensitivity analysis — rank correlation heatmap
#                  (Weighted Stouffer vs. Fisher vs. Random Effects)
#   - Figure 3C  : Leave-One-Out (LOO) meta-analysis stability scatter plots
#
# Input file:
#   - Weighted_Stouffer_Meta_Analysis.xlsx
#     Columns required: Gene, meta_log2FC, meta_padj, meta_pval,
#                       Log2FC_C26, Log2FC_KIC, Log2FC_LLC,
#                       pval_C26, pval_KIC, pval_LLC
#
# Author : Dae-Hwan Kim, Jebeom Ko
# Paper  : "Diagnostic Platform for Cancer Cachexia: Multi-Omics Based Novel Biomarker Discovery"
#
# Public datasets used (LC-MS/MS proteomics, Mus musculus, skeletal muscle):
#   - PXD027490 : C26 model   | CD2F1, male, ~6 weeks
#   - PXD036752 : KIC model   | Mixed background (C57BL/6 + FVB), male, 6-7 weeks
#   - PXD016474 : LLC model   | C57BL/6, male, 6 weeks
# ============================================================================

# ----------------------------------------------------------------------------
# 0. Load Required Packages
# ----------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(readxl)
  library(ggplot2)
  library(ggrepel)
  library(dplyr)
  library(tidyr)
  library(pheatmap)
  library(patchwork)
})

set.seed(1234)

# ----------------------------------------------------------------------------
# 1. Data Paths Configuration
# ----------------------------------------------------------------------------
# NOTE: Modify these paths according to your local environment.
#
# Input:  Weighted Stouffer meta-analysis result table (xlsx)
# Output: PDF and PNG figures saved to OUTPUT_DIR

INPUT_FILE <- "/Users/daehwankim/Desktop/KIST_folder/Progress/Cachexia model/Serpina3_Paper/260119/Weighted_Stouffer_Meta_Analysis.xlsx"
OUTPUT_DIR <- "/Users/daehwankim/Desktop/KIST_folder/Progress/Cachexia model/Serpina3_Paper/260119"

# Create output directory if needed
if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)

# Helper: build output file path
out <- function(filename) file.path(OUTPUT_DIR, filename)

# ----------------------------------------------------------------------------
# 2. Load Data
# ----------------------------------------------------------------------------
cat("Loading dataset...\n")
df <- read_excel(INPUT_FILE)

# Candidate secretome genes to highlight across all plots
candidate_genes <- c(
  "Serpina3m", "Col19a1", "Serpina3n", "C4b", "Lox",
  "Gpx3", "Fibin", "Igfbp3", "Tnfrsf12a", "Ctsl", "Vegfd"
)

# Model labels and sample-size-based weights (for Stouffer LOO)
# Weights are proportional to sample sizes: C26 (n=10, PXD027490), KIC (n=6, PXD036752), LLC (n=7, PXD016474)
models <- c("C26", "KIC", "LLC")
weights_map <- c(C26 = 10, KIC = 6, LLC = 7) # adjust to actual n per model

# ----------------------------------------------------------------------------
# 3. Figure 3A: Meta-Analysis Volcano Plot (full y-axis)
# ----------------------------------------------------------------------------
cat("Generating Figure 3A — Meta-Analysis Volcano Plot...\n")

df_vol <- df %>%
  select(Gene, meta_log2FC, meta_padj) %>%
  drop_na() %>%
  mutate(
    Significance = ifelse(
      meta_log2FC >= 1 & meta_padj <= 0.05,
      "Significant", "Not Significant"
    )
  )

# Candidate genes that pass significance threshold (to label)
label_genes <- df_vol %>%
  filter(
    tolower(Gene) %in% tolower(candidate_genes),
    Significance == "Significant"
  )

cat("  Labeled candidate genes:", nrow(label_genes), "\n")
if (nrow(label_genes) > 0) {
  cat("  ", paste(label_genes$Gene, collapse = ", "), "\n")
}

fig3A <- ggplot(
  df_vol,
  aes(
    x = meta_log2FC, y = -log10(meta_padj),
    color = Significance, size = Significance
  )
) +
  geom_point(alpha = 0.5) +
  scale_color_manual(
    values = c("Not Significant" = "grey80", "Significant" = "red")
  ) +
  scale_size_manual(
    values = c("Not Significant" = 1.5, "Significant" = 2.5),
    guide  = "none"
  ) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey50") +
  geom_text_repel(
    data = label_genes,
    aes(label = Gene),
    size = 4,
    color = "black",
    fontface = "italic",
    box.padding = 0.7,
    max.overlaps = Inf,
    min.segment.length = 0
  ) +
  labs(
    title = "Meta-Analysis Volcano Plot",
    x     = bquote(~Meta ~ Log[2] ~ "Fold Change"),
    y     = bquote(~ -Log[10] ~ "(Meta Adjusted p-value)")
  ) +
  theme_minimal() +
  theme(
    legend.position = "none",
    plot.title      = element_text(hjust = 0.5, face = "bold", size = 14),
    axis.title      = element_text(size = 12)
  )

print(fig3A)
ggsave(out("Figure3A_MetaVolcano.pdf"), fig3A, width = 6, height = 5)
ggsave(out("Figure3A_MetaVolcano.png"), fig3A, width = 6, height = 5, dpi = 300)
cat("  -> Saved Figure3A_MetaVolcano.pdf / .png\n")

# ----------------------------------------------------------------------------
# 4. Figure 3A': Y-axis Capped Volcano Plot
# ----------------------------------------------------------------------------
cat("Generating Figure 3A' — Y-axis capped Volcano Plot...\n")

Y_CAP <- 8 # cap -log10(padj) at this value (triangles mark capped points)

df_vol_cap <- df_vol %>%
  mutate(
    neg_log10_padj = -log10(meta_padj),
    is_capped      = neg_log10_padj > Y_CAP,
    neg_log10_padj = ifelse(is_capped, Y_CAP, neg_log10_padj),
    point_shape    = ifelse(is_capped, "capped", "normal")
  )

label_genes_cap <- df_vol_cap %>%
  filter(
    tolower(Gene) %in% tolower(candidate_genes),
    Significance == "Significant"
  )

fig3A_cap <- ggplot(
  df_vol_cap,
  aes(
    x = meta_log2FC, y = neg_log10_padj,
    color = Significance, shape = point_shape, size = Significance
  )
) +
  geom_point(alpha = 0.5) +
  scale_color_manual(
    values = c("Not Significant" = "grey80", "Significant" = "red")
  ) +
  scale_shape_manual(
    values = c("normal" = 16, "capped" = 17),
    guide  = "none"
  ) +
  scale_size_manual(
    values = c("Not Significant" = 1.5, "Significant" = 2.5),
    guide  = "none"
  ) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey50") +
  geom_text_repel(
    data = label_genes_cap,
    aes(label = Gene),
    size = 4,
    color = "black",
    fontface = "italic",
    box.padding = 0.7,
    max.overlaps = Inf,
    min.segment.length = 0
  ) +
  scale_y_continuous(
    limits = c(0, Y_CAP + 0.3),
    breaks = seq(0, Y_CAP, by = 2)
  ) +
  labs(
    title = "Meta-Analysis Volcano Plot",
    x     = bquote(~Meta ~ Log[2] ~ "Fold Change"),
    y     = bquote(~ -Log[10] ~ "(Meta Adjusted p-value)")
  ) +
  theme_classic() +
  theme(
    legend.position = "none",
    plot.title      = element_text(hjust = 0.5, face = "bold", size = 14),
    axis.title      = element_text(size = 12)
  )

print(fig3A_cap)
ggsave(out("Figure3A_MetaVolcano_capped.pdf"), fig3A_cap, width = 6, height = 5)
ggsave(out("Figure3A_MetaVolcano_capped.png"), fig3A_cap, width = 6, height = 5, dpi = 300)
cat("  -> Saved Figure3A_MetaVolcano_capped.pdf / .png\n")
cat("     (▲ = points with -log10(padj) >", Y_CAP, "capped)\n")

# ----------------------------------------------------------------------------
# 5. Figure 3B: Sensitivity Analysis — Rank Correlation Heatmap
# ----------------------------------------------------------------------------
cat("Generating Figure 3B — Sensitivity Analysis...\n")

# -- 5-1. Helper: Fisher's combined probability method
fisher_meta <- function(pvals, log2fcs) {
  valid <- !is.na(pvals) & !is.na(log2fcs)
  pvals <- pvals[valid]
  log2fcs <- log2fcs[valid]
  k <- length(pvals)
  if (k == 0) {
    return(list(pval = NA, log2fc = NA))
  }

  # Fisher: -2 * sum(log(p)) ~ chi-sq(2k)
  pvals[pvals == 0] <- .Machine$double.xmin
  chi_sq <- -2 * sum(log(pvals))
  fisher_p <- pchisq(chi_sq, df = 2 * k, lower.tail = FALSE)

  list(pval = fisher_p, log2fc = mean(log2fcs))
}

# -- 5-2. Helper: DerSimonian-Laird random-effects meta-analysis
random_effects_meta <- function(log2fcs, pvals) {
  valid <- !is.na(log2fcs) & !is.na(pvals)
  log2fcs <- log2fcs[valid]
  pvals <- pvals[valid]
  k <- length(log2fcs)
  if (k < 2) {
    return(list(
      log2fc = ifelse(k == 1, log2fcs, NA),
      pval   = ifelse(k == 1, pvals, NA)
    ))
  }

  # Estimate SE from p-value
  z_obs <- abs(qnorm(pvals / 2))
  z_obs[is.infinite(z_obs)] <- 10
  z_obs[z_obs == 0] <- 0.001
  se <- abs(log2fcs) / z_obs
  se[se == 0] <- 0.001

  # Fixed-effect weights & estimate
  w <- 1 / se^2
  theta_fe <- sum(w * log2fcs) / sum(w)

  # Cochran's Q → DerSimonian-Laird tau^2
  Q <- sum(w * (log2fcs - theta_fe)^2)
  df_val <- k - 1
  c_val <- sum(w) - sum(w^2) / sum(w)
  tau2 <- max(0, (Q - df_val) / c_val)

  # Random-effects estimate and p-value
  w_re <- 1 / (se^2 + tau2)
  theta_re <- sum(w_re * log2fcs) / sum(w_re)
  se_re <- sqrt(1 / sum(w_re))
  z_re <- theta_re / se_re
  p_re <- 2 * pnorm(-abs(z_re))

  list(log2fc = theta_re, pval = p_re)
}

# -- 5-3. Compute all three methods for every gene
cat("  Computing Fisher and Random-effects results for each gene...\n")

method_results <- lapply(seq_len(nrow(df)), function(i) {
  row <- df[i, ]
  log2fcs <- c(row$Log2FC_C26, row$Log2FC_KIC, row$Log2FC_LLC)
  pvals <- c(row$pval_C26, row$pval_KIC, row$pval_LLC)

  fish <- fisher_meta(pvals, log2fcs)
  re <- random_effects_meta(log2fcs, pvals)

  data.frame(
    Gene = row$Gene,
    Stouffer_log2FC = row$meta_log2FC,
    Stouffer_pval = row$meta_pval,
    Stouffer_padj = row$meta_padj,
    Fisher_log2FC = fish$log2fc,
    Fisher_pval = fish$pval,
    RE_log2FC = re$log2fc,
    RE_pval = re$pval,
    stringsAsFactors = FALSE
  )
})

method_df <- bind_rows(method_results)
method_df$Fisher_padj <- p.adjust(method_df$Fisher_pval, method = "BH")
method_df$RE_padj <- p.adjust(method_df$RE_pval, method = "BH")

# -- 5-4. Rank correlation heatmap
cat("  Computing Spearman rank correlation matrix...\n")

method_df_complete <- method_df %>%
  filter(!is.na(Stouffer_pval) & !is.na(Fisher_pval) & !is.na(RE_pval))

rank_stouffer <- rank(
  -method_df_complete$Stouffer_log2FC / method_df_complete$Stouffer_pval,
  na.last = "keep"
)
rank_fisher <- rank(
  -method_df_complete$Fisher_log2FC / method_df_complete$Fisher_pval,
  na.last = "keep"
)
rank_re <- rank(
  -method_df_complete$RE_log2FC / method_df_complete$RE_pval,
  na.last = "keep"
)

cor_mat <- cor(
  cbind(
    `Weighted\nStouffer` = rank_stouffer,
    `Fisher`             = rank_fisher,
    `Random\nEffects`    = rank_re
  ),
  use = "complete.obs",
  method = "spearman"
)

heatmap_args <- list(
  cor_mat,
  color           = colorRampPalette(c("white", "steelblue", "darkblue"))(50),
  breaks          = seq(0.5, 1, length.out = 51),
  display_numbers = TRUE,
  number_format   = "%.3f",
  cluster_rows    = FALSE,
  cluster_cols    = FALSE,
  main            = "Spearman Rank Correlation between Methods",
  fontsize_number = 14,
  angle_col       = 0
)

pdf(out("Figure3B_Sensitivity_RankCorrelation.pdf"), width = 5, height = 4.5)
do.call(pheatmap, heatmap_args)
dev.off()

png(out("Figure3B_Sensitivity_RankCorrelation.png"),
  width = 5 * 300, height = 4.5 * 300, res = 300
)
do.call(pheatmap, heatmap_args)
dev.off()

cat("  -> Saved Figure3B_Sensitivity_RankCorrelation.pdf / .png\n")

# ----------------------------------------------------------------------------
# 6. Figure 3C: Leave-One-Out (LOO) Meta-Analysis Stability
# ----------------------------------------------------------------------------
cat("Generating Figure 3C — Leave-One-Out Meta-Analysis...\n")

# -- 6-1. Helper: Weighted Stouffer for a subset of models
stouffer_loo <- function(log2fcs, pvals, ws) {
  valid <- !is.na(log2fcs) & !is.na(pvals)
  log2fcs <- log2fcs[valid]
  pvals <- pvals[valid]
  ws <- ws[valid]
  if (length(log2fcs) == 0) {
    return(list(log2fc = NA, z = NA, pval = NA))
  }

  meta_log2fc <- sum(log2fcs * ws) / sum(ws)

  z_scores <- qnorm(1 - pvals / 2) * sign(log2fcs)
  z_scores[is.infinite(z_scores)] <- sign(z_scores[is.infinite(z_scores)]) * 10
  z_scores[is.nan(z_scores)] <- 0

  meta_z <- sum(ws * z_scores) / sqrt(sum(ws^2))
  meta_p <- 2 * pnorm(-abs(meta_z))

  list(log2fc = meta_log2fc, z = meta_z, pval = meta_p)
}

# -- 6-2. Compute LOO for each gene
loo_results <- lapply(seq_len(nrow(df)), function(i) {
  row <- df[i, ]
  log2fcs <- c(row$Log2FC_C26, row$Log2FC_KIC, row$Log2FC_LLC)
  pvals <- c(row$pval_C26, row$pval_KIC, row$pval_LLC)
  ws <- c(weights_map["C26"], weights_map["KIC"], weights_map["LLC"])

  results <- list()
  for (j in seq_along(models)) {
    res <- stouffer_loo(log2fcs[-j], pvals[-j], ws[-j])
    results[[paste0("LOO_", models[j])]] <- res
  }

  data.frame(
    Gene             = row$Gene,
    meta_log2FC      = row$meta_log2FC,
    meta_padj        = row$meta_padj,
    LOO_noC26_log2FC = results$LOO_C26$log2fc,
    LOO_noC26_pval   = results$LOO_C26$pval,
    LOO_noKIC_log2FC = results$LOO_KIC$log2fc,
    LOO_noKIC_pval   = results$LOO_KIC$pval,
    LOO_noLLC_log2FC = results$LOO_LLC$log2fc,
    LOO_noLLC_pval   = results$LOO_LLC$pval,
    stringsAsFactors = FALSE
  )
})

loo_df <- bind_rows(loo_results)
loo_df$LOO_noC26_padj <- p.adjust(loo_df$LOO_noC26_pval, method = "BH")
loo_df$LOO_noKIC_padj <- p.adjust(loo_df$LOO_noKIC_pval, method = "BH")
loo_df$LOO_noLLC_padj <- p.adjust(loo_df$LOO_noLLC_pval, method = "BH")

cat("  LOO analysis complete for", nrow(loo_df), "genes.\n")

# -- 6-3. Helper: build one LOO scatter panel
make_loo_panel <- function(loo_col, exclude_model) {
  plot_df <- loo_df %>%
    filter(!is.na(meta_log2FC), !is.na(.data[[loo_col]])) %>%
    mutate(
      Significance = ifelse(
        meta_log2FC >= 1 & meta_padj <= 0.05,
        "Significant", "Not Significant"
      )
    )

  r_val <- cor(plot_df$meta_log2FC, plot_df[[loo_col]], use = "complete.obs")

  ggplot(plot_df, aes(x = meta_log2FC, y = .data[[loo_col]], color = Significance)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
    geom_point(alpha = 0.4, size = 1.2) +
    scale_color_manual(
      values = c("Not Significant" = "grey80", "Significant" = "red")
    ) +
    annotate(
      "text",
      x = min(plot_df$meta_log2FC, na.rm = TRUE) + 0.3,
      y = max(plot_df[[loo_col]], na.rm = TRUE) - 0.3,
      label = paste0("r = ", round(r_val, 3)),
      size = 4, fontface = "bold", color = "steelblue"
    ) +
    labs(
      title = paste0("Without ", exclude_model),
      x     = bquote("All Models " ~ Log[2] ~ "FC"),
      y     = bquote("LOO " ~ Log[2] ~ "FC")
    ) +
    theme_classic() +
    theme(
      legend.position = "none",
      plot.title      = element_text(hjust = 0.5, face = "bold", size = 12),
      axis.title      = element_text(size = 10)
    )
}

# -- 6-4. Combine three panels with patchwork
p_loo1 <- make_loo_panel("LOO_noC26_log2FC", "C26")
p_loo2 <- make_loo_panel("LOO_noKIC_log2FC", "KIC")
p_loo3 <- make_loo_panel("LOO_noLLC_log2FC", "LLC")

fig3C <- p_loo1 + p_loo2 + p_loo3 +
  plot_annotation(
    title = "Leave-One-Out Meta-Analysis Stability",
    subtitle = "Each panel excludes one model. Points on diagonal = no change.",
    theme = theme(
      plot.title    = element_text(hjust = 0.5, face = "bold", size = 14),
      plot.subtitle = element_text(hjust = 0.5, size = 10, color = "grey40")
    )
  )

print(fig3C)
ggsave(out("Figure3C_LOO_Scatter.pdf"), fig3C, width = 14, height = 5)
ggsave(out("Figure3C_LOO_Scatter.png"), fig3C, width = 14, height = 5, dpi = 300)
cat("  -> Saved Figure3C_LOO_Scatter.pdf / .png\n")

# ----------------------------------------------------------------------------
# 7. Save Session Info
# ----------------------------------------------------------------------------
writeLines(
  capture.output(sessionInfo()),
  out("session_info_Figure3.txt")
)

cat("\n============================================\n")
cat("Figure 3 analysis completed!\n")
cat("Output files saved to:", OUTPUT_DIR, "\n")
cat("============================================\n")
