# ============================================================================
# Figure 5: Human Patient Bulk RNA-seq Analysis (GSE254877)
# ============================================================================
# This script performs:
#   - Figure 5B : IntNMF clustering (K=2) + Silhouette plot
#   - Figure 5C : DESeq2 (C2 vs C1, adjusted for cancer type) +
#                 GO BP ORA bubble plot (inflammation/immune terms)
#   - Figure 5D : Correlation plots — CTSL & SERPINA3 vs. atrophy markers
#                 (EDA2R, FBXO32 only; 4 scatter plots total)
#
# Dataset : GSE254877
#           Pancreas (n=38) + Colorectal (n=46) cancer skeletal muscle
# Normalization : DESeq2 VST (blind = TRUE, design ~ 1)
#
# Author : Dae-Hwan Kim, Jebeom Ko
# Paper  : "Diagnostic Platform for Cancer Cachexia: Multi-Omics Based Novel Biomarker Discovery"
#
# Public dataset used (Bulk RNA-seq, Homo sapiens, Rectus abdominis):
#   - GSE254877 : Colorectal (n=46) + Pancreatic (n=38) cancer patients, male and female
# ============================================================================

# ----------------------------------------------------------------------------
# 0. Load Required Packages
# ----------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(tximport)
  library(biomaRt)
  library(DESeq2)
  library(IntNMF)
  library(cluster)
  library(readxl)
  library(pheatmap)
  library(dplyr)
  library(tibble)
  library(tidyr)
  library(matrixStats)
  library(ggplot2)
  library(ggalluvial)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(fgsea)
  library(msigdbr)
  library(sva)
  library(apeglm)
  library(rhdf5)
})

set.seed(1234)

# ----------------------------------------------------------------------------
# 1. Data Paths Configuration
# ----------------------------------------------------------------------------
# NOTE: Modify these paths according to your local environment.
#
# Required:
#   - meta_data.xlsx  : Sample metadata (Sheet2: srr_id, sample_name, group)
#   - Kallisto results: <results_dir>/<srr_id>_abundance.tsv per sample

BASE_DIR    <- ""
RESULTS_DIR <- file.path(BASE_DIR, "Kallisto_Results")
META_FILE   <- file.path(BASE_DIR, "GSE254877", "meta_data.xlsx")
OUTPUT_DIR  <- file.path(BASE_DIR, "output", "kmeans_K2")

if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)

# Target genes for correlation analysis (Fig 5D)
# NOTE: TRIM63 excluded per analysis design (EDA2R + FBXO32 only)
TARGET_GENES <- c("CTSL", "FBXO32", "EDA2R", "SERPINA3")

# ----------------------------------------------------------------------------
# 2. Load Metadata
# ----------------------------------------------------------------------------
meta <- read_excel(META_FILE, sheet = "Sheet2")
colnames(meta) <- c("srr_id", "sample_name", "group")
meta$group <- factor(meta$group, levels = c("Pancreas", "Colorectal"))


# ----------------------------------------------------------------------------
# 3. Build Transcript-to-Gene Mapping via biomaRt
# ----------------------------------------------------------------------------

tx2gene <- tryCatch(
  {
    ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
    t2g <- getBM(
      attributes = c("ensembl_transcript_id_version", "external_gene_name"),
      mart = ensembl
    )
    colnames(t2g) <- c("TXNAME", "GENEID")
    t2g <- t2g[t2g$GENEID != "", ]
    t2g$TXNAME <- gsub("\\.\\d+$", "", t2g$TXNAME)  # strip version suffix
    t2g <- t2g[!duplicated(t2g), ]
    cat("  biomaRt: retrieved", nrow(t2g), "transcript-gene mappings\n")
    t2g
  },
  error = function(e) {
    stop("biomaRt failed. Check internet connection: ", conditionMessage(e))
  }
)

# ----------------------------------------------------------------------------
# 4. Import Kallisto Results with tximport
# ----------------------------------------------------------------------------

files <- file.path(RESULTS_DIR, paste0(meta$srr_id, "_abundance.tsv"))
names(files) <- meta$srr_id

missing <- !file.exists(files)
if (any(missing)) {
  cat("  WARNING: Missing files for:", paste(names(files)[missing], collapse = ", "), "\n")
  meta  <- meta[!missing, ]
  files <- files[!missing]
}

txi <- tximport(
  files,
  type                = "kallisto",
  tx2gene             = tx2gene,
  txOut               = FALSE,
  ignoreTxVersion     = TRUE,
  countsFromAbundance = "no"
)

# ----------------------------------------------------------------------------
# 5. DESeq2 VST Normalization (design ~ 1, for unsupervised clustering)
# ----------------------------------------------------------------------------

dds <- DESeqDataSetFromTximport(txi, colData = meta, design = ~1)
dds <- dds[rowSums(counts(dds) >= 10) >= 10, ]  # pre-filter low-count genes
vsd <- vst(dds, blind = TRUE)
vst_mat <- assay(vsd)
cat("  VST complete:", nrow(vst_mat), "genes x", ncol(vst_mat), "samples\n")

# ----------------------------------------------------------------------------
# 6. Figure 5B: IntNMF Clustering (Optimal K Selection + K=2)
# ----------------------------------------------------------------------------

# -- 6-1. Select highly variable genes (top 2000)
rv     <- rowVars(vst_mat)
top_ix <- order(rv, decreasing = TRUE)[seq_len(min(2000, length(rv)))]

# -- 6-2. Prepare data for IntNMF (samples x genes, non-negative)
dat_rna <- t(vst_mat[top_ix, ])
if (!all(dat_rna >= 0)) {
  dat_rna <- pmax(dat_rna + abs(min(dat_rna)), .Machine$double.eps)
}
dat_rna   <- dat_rna / max(dat_rna)
dat_list  <- list(rna = dat_rna)

# -- 6-3. Find optimal K (CPI metric, 2:5)
set.seed(1234)
cpi_100 <- nmf.opt.k(
  dat = dat_list, k.range = 2:5, n.fold = 5, n.runs = 30,
  maxiter = 100, st.count = 10, make.plot = TRUE, progress = TRUE
)
set.seed(1234)
cpi_300 <- nmf.opt.k(
  dat = dat_list, k.range = 2:5, n.fold = 5, n.runs = 30,
  maxiter = 300, st.count = 10, make.plot = TRUE, progress = TRUE
)


# -- 6-4. Final clustering (K = 2)
K_final <- 2
fit <- nmf.mnnals(
  dat = dat_list, k = K_final, maxiter = 300, st.count = 20,
  n.ini = 30, ini.nndsvd = TRUE, seed = TRUE
)
meta$cluster_intNMF <- paste0("C", fit$clusters)
cat("  Cluster distribution:\n")
print(table(meta$cluster_intNMF, meta$group))

# -- 6-5. Silhouette plot (ggplot2)
dist_mat <- dist(dat_rna, method = "euclidean")
sil      <- silhouette(fit$clusters, dist_mat)
sil_df   <- as.data.frame(sil[, 1:3])
colnames(sil_df) <- c("cluster", "neighbor", "sil_width")

sil_df$cluster <- factor(paste0("C", sil_df$cluster))
sil_df         <- sil_df[order(sil_df$cluster, -sil_df$sil_width), ]
sil_df$index   <- seq_len(nrow(sil_df))
avg_sil        <- mean(sil_df$sil_width)

fig5B <- ggplot(sil_df, aes(x = index, y = sil_width, fill = cluster)) +
  geom_bar(stat = "identity", width = 1, color = NA) +
  geom_hline(
    yintercept = avg_sil, linetype = "dashed",
    color = "grey30", linewidth = 0.6
  ) +
  annotate(
    "text",
    x     = nrow(sil_df) / 2,
    y     = avg_sil + 0.02,
    label = sprintf("Mean = %.3f", avg_sil),
    color = "grey30", fontface = "bold", vjust = 0
  ) +
  scale_fill_manual(values = c(C1 = "#56B4E9", C2 = "#E69F00"), name = "Cluster") +
  scale_y_continuous(limits = c(-0.1, 1.1), expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  coord_cartesian(ylim = c(0, 1.1)) +
  labs(
    title    = paste0("Silhouette Plot (K = ", K_final, ")"),
    subtitle = "Euclidean Distance",
    x        = "Samples",
    y        = "Silhouette Width"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title        = element_text(face = "bold", hjust = 0.5),
    plot.subtitle     = element_text(hjust = 0.5, color = "grey40"),
    axis.text.x       = element_blank(),
    axis.ticks.x      = element_blank(),
    axis.ticks.y      = element_line(color = "grey30", linewidth = 0.6),
    axis.line.y       = element_line(color = "grey30", linewidth = 0.6),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    legend.position   = "bottom"
  )


# ----------------------------------------------------------------------------
# 7. Figure 5C: DESeq2 (C2 vs C1) + GO BP ORA Bubble Plot
# ----------------------------------------------------------------------------

# -- 7-1. DESeq2 (adjusted for cancer type group)
meta_de          <- as.data.frame(meta)
rownames(meta_de) <- meta_de$srr_id
stopifnot(all(colnames(txi$counts) %in% rownames(meta_de)))
meta_de          <- meta_de[colnames(txi$counts), , drop = FALSE]
meta_de$cluster_raw <- factor(meta_de$cluster_intNMF, levels = c("C1", "C2"))

dds_de <- DESeqDataSetFromTximport(txi, colData = meta_de, design = ~ group + cluster_raw)
dds_de <- dds_de[rowSums(counts(dds_de) >= 10) >= 10, ]
dds_de <- DESeq(dds_de, quiet = TRUE)

# Result: cluster_raw_C2_vs_C1 (positive = higher in C2)
res_shr <- lfcShrink(dds_de, coef = "cluster_raw_C2_vs_C1", type = "apeglm")

res_df <- as.data.frame(res_shr) %>%
  rownames_to_column("gene") %>%
  arrange(padj)


# -- 7-2. GO BP ORA (C1_up vs C2_up)
padj_cut <- 0.05
lfc_cut  <- 1

sig_df <- res_df %>%
  filter(!is.na(padj), padj < padj_cut,
         !is.na(log2FoldChange), abs(log2FoldChange) >= lfc_cut)

# C2_vs_C1: positive = C2_up, negative = C1_up
genes_c2_up <- sig_df %>% filter(log2FoldChange >= lfc_cut)  %>% pull(gene) %>% unique()
genes_c1_up <- sig_df %>% filter(log2FoldChange <= -lfc_cut) %>% pull(gene) %>% unique()
cat("  C2_up genes:", length(genes_c2_up), " | C1_up genes:", length(genes_c1_up), "\n")

ego_c2 <- enrichGO(
  gene = genes_c2_up, OrgDb = org.Hs.eg.db, keyType = "SYMBOL",
  ont = "BP", pAdjustMethod = "BH", readable = TRUE
)
ego_c1 <- enrichGO(
  gene = genes_c1_up, OrgDb = org.Hs.eg.db, keyType = "SYMBOL",
  ont = "BP", pAdjustMethod = "BH", readable = TRUE
)

df_c2 <- as.data.frame(ego_c2)
df_c1 <- as.data.frame(ego_c1)

# -- 7-3. Select representative terms for bubble plot
custom_terms <- c(
  "regulation of autophagy",
  "regulation of protein catabolic process",
  "response to oxidative stress",
  "regulation of innate immune response",
  "canonical NF-kappaB signal transduction"
)

pick_terms <- function(df, direction) {
  df %>%
    dplyr::filter(Description %in% custom_terms) %>%
    dplyr::transmute(
      direction  = direction,
      term       = Description,
      padj       = p.adjust,
      GeneRatio  = GeneRatio
    )
}

plot_df <- dplyr::bind_rows(
  pick_terms(df_c1, "C1_up"),
  pick_terms(df_c2, "C2_up")
) %>%
  tidyr::complete(
    direction = c("C1_up", "C2_up"),
    term      = custom_terms
  ) %>%
  dplyr::mutate(
    term      = factor(term, levels = rev(custom_terms)),
    direction = factor(direction, levels = c("C1_up", "C2_up")),
    GeneRatio_num = dplyr::if_else(
      is.na(GeneRatio), NA_real_,
      as.numeric(sub("/.*", "", GeneRatio)) / as.numeric(sub(".*/", "", GeneRatio))
    )
  )

fig5C <- ggplot(plot_df, aes(x = direction, y = term)) +
  geom_point(
    data  = dplyr::filter(plot_df, !is.na(padj)),
    aes(size = GeneRatio_num, color = padj),
    alpha = 0.9
  ) +
  scale_x_discrete(drop = FALSE, labels = c(C1_up = "C1", C2_up = "C2")) +
  scale_color_gradient(low = "#B2182B", high = "#2166AC") +
  scale_size_continuous(range = c(2, 10)) +
  labs(title = "GO BP ORA", color = "adj p", size = "GeneRatio") +
  theme_classic(base_size = 13) +
  theme(
    axis.title   = element_blank(),
    axis.text.x  = element_text(face = "bold"),
    axis.text.y  = element_text(size = 11),
    legend.position = "right"
  )


# ----------------------------------------------------------------------------
# 8. Figure 5D: Correlation Plots (CTSL & SERPINA3 vs EDA2R & FBXO32)
# ----------------------------------------------------------------------------
# 4 scatter plots:
#   CTSL    vs EDA2R    |  CTSL    vs FBXO32
#   SERPINA3 vs EDA2R   |  SERPINA3 vs FBXO32
# NOTE: TRIM63 excluded from this analysis.
# ----------------------------------------------------------------------------
cat("\nGenerating Figure 5D — Correlation Plots...\n")

# -- 8-1. Re-run VST with design ~1 (blind) for correlation analysis
sampleTable <- data.frame(row.names = meta$srr_id, group = meta$group)
dds_corr    <- DESeqDataSetFromTximport(txi, colData = sampleTable, design = ~1)
keep        <- rowSums(counts(dds_corr) >= 10) >= 10
dds_corr    <- dds_corr[keep, ]
vsd_corr    <- vst(dds_corr, blind = TRUE)
vst_corr    <- assay(vsd_corr)


# -- 8-3. Build expression data frame
expr_df <- data.frame(
  srr_id   = colnames(vst_corr),
  CTSL     = as.numeric(vst_corr["CTSL",     ]),
  FBXO32   = as.numeric(vst_corr["FBXO32",   ]),
  EDA2R    = as.numeric(vst_corr["EDA2R",     ]),
  SERPINA3 = as.numeric(vst_corr["SERPINA3",  ]),
  stringsAsFactors = FALSE
)
expr_df <- merge(expr_df, meta, by = "srr_id")

# -- 8-4. Group-adjusted residuals (robustness check)
cat("\nComputing group-adjusted residuals...\n")
expr_df$CTSL_resid     <- residuals(lm(CTSL     ~ group, data = expr_df))
expr_df$FBXO32_resid   <- residuals(lm(FBXO32   ~ group, data = expr_df))
expr_df$EDA2R_resid    <- residuals(lm(EDA2R    ~ group, data = expr_df))
expr_df$SERPINA3_resid <- residuals(lm(SERPINA3 ~ group, data = expr_df))

# -- 8-5. Correlation scatter plot helper
make_corr_plot <- function(df, x_col, y_col, x_label, y_label,
                           title, color_col = "group") {
  cor_p <- cor.test(df[[x_col]], df[[y_col]], method = "pearson")
  cor_s <- cor.test(df[[x_col]], df[[y_col]], method = "spearman")

  subtitle <- sprintf(
    "Pearson r = %.3f (p = %.2e)  |  Spearman rho = %.3f (p = %.2e)",
    cor_p$estimate, cor_p$p.value,
    cor_s$estimate, cor_s$p.value
  )

  p <- ggplot(df, aes(x = .data[[x_col]], y = .data[[y_col]])) +
    geom_smooth(
      method   = "lm", se = TRUE, color = "grey40",
      linetype = "dashed", linewidth = 0.8
    ) +
    labs(title = title, subtitle = subtitle, x = x_label, y = y_label) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title    = element_text(face = "bold", size = 16),
      plot.subtitle = element_text(size = 10, color = "grey30"),
      legend.position = "bottom"
    )

  if (!is.null(color_col) && color_col %in% colnames(df)) {
    p <- p +
      geom_jitter(
        aes(color = .data[[color_col]]),
        size = 2, alpha = 0.7, shape = 16, width = 0.15, height = 0.15
      ) +
      scale_color_manual(
        values = c("Pancreas" = "#E74C3C", "Colorectal" = "#3498DB"),
        name   = "Cancer Type"
      )
  } else {
    p <- p +
      geom_jitter(
        size = 2, alpha = 0.7, color = "#2C3E50",
        shape = 16, width = 0.15, height = 0.15
      )
  }
  p
}

# -- 8-6. Generate 4 scatter plots (2 reference genes × 2 atrophy markers)
p_ctsl_eda2r    <- make_corr_plot(expr_df, "CTSL",     "EDA2R",
                                  "CTSL (VST)",     "EDA2R (VST)",
                                  "CTSL vs EDA2R")
p_ctsl_fbxo32   <- make_corr_plot(expr_df, "CTSL",     "FBXO32",
                                  "CTSL (VST)",     "FBXO32 (VST)",
                                  "CTSL vs FBXO32")
p_ser_eda2r     <- make_corr_plot(expr_df, "SERPINA3", "EDA2R",
                                  "SERPINA3 (VST)", "EDA2R (VST)",
                                  "SERPINA3 vs EDA2R")
p_ser_fbxo32    <- make_corr_plot(expr_df, "SERPINA3", "FBXO32",
                                  "SERPINA3 (VST)", "FBXO32 (VST)",
                                  "SERPINA3 vs FBXO32")




