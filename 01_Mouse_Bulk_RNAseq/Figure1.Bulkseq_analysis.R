# ============================================================================
# Figure 1: Integrated Bulk RNA-seq Analysis for Cancer Cachexia
# ============================================================================
# This script performs:
#   - Figure 1B: PCA plot
#   - Figure 1C: Volcano plot  
#   - Figure 1D: K-means clustering heatmap (GC1-GC5)
#   - Figure 1E: GO Over-Representation Analysis (ORA) dot plot
#
# Author: Dae-Hwan Kim, Jibeom Ko
# Paper: "Diagnostic Platform for Cancer Cachexia: Multi-Omics Based Novel Biomarker Discovery"
# ============================================================================

# ----------------------------------------------------------------------------
# 0. Load Required Packages
# ----------------------------------------------------------------------------
library(dplyr)
library(ggrepel)
library(tximport)
library(DESeq2)
library(pheatmap)
library(openxlsx)
library(sva)
library(enrichplot)
library(ggplot2)
library(tibble)
library(cluster)
library(clusterProfiler)
library(org.Mm.eg.db)
library(biomaRt)

set.seed(1234)

# ----------------------------------------------------------------------------
# 1. Data Paths Configuration
# ----------------------------------------------------------------------------
# NOTE: Modify these paths according to your local environment
#
# Required data files:
#   - TXNAME: Transcript ID file (tab-delimited, no header)
#   - SYMBOL: Gene symbol file (tab-delimited, no header)
#   - Sample .tsv files: Kallisto quantification outputs
#
# Public datasets used (Bulk RNA-seq, Mus musculus, skeletal muscle):
#   - GSE65936  : CD2F1, male, ~6 weeks,   C26 tumor (Gastrocnemius)
#   - GSE123310 : C57BL/6J, male, 8 weeks, LLC tumor (Tibialis anterior)
#   - GSE138464 : C57BL/6, male, 24 weeks, KIC pancreatic cancer (Quadriceps)
#   - GSE142455 : Mixed background (C57BL/6 + FVB), male, 6-7 weeks, KIC model

DATA_DIR <- "path/to/your/data"  # Change this to your data directory
TXNAME_PATH <- file.path("/Users/daehwankim/Desktop/Desktop_MacBook Air/sequencing data/Bulk_seq_practice/TXNAME")
SYMBOL_PATH <- file.path("/Users/daehwankim/Desktop/Desktop_MacBook Air/sequencing data/Bulk_seq_practice/SYMBOL")
KALLISTO_DIR <- file.path("/Users/daehwankim/Desktop/Seqencing_Practicing/Bulk_seq_analysis_practice/kallisto")
OUTPUT_DIR <- file.path(DATA_DIR, "results")

# Create output directory if not exists
if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)

# ----------------------------------------------------------------------------
# 2. Load Transcript-to-Gene Mapping
# ----------------------------------------------------------------------------
TXNAME <- read.delim(TXNAME_PATH, header = FALSE)
SYMBOL <- read.delim(SYMBOL_PATH, header = FALSE)
tx2gene <- data.frame(TXNAME, SYMBOL)
colnames(tx2gene) <- c("TXNAME", "SYMBOL")

# ----------------------------------------------------------------------------
# 3. Filter for Protein-Coding Genes Only
# ----------------------------------------------------------------------------
# Connect to Ensembl database (requires internet connection)
mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")

# Retrieve gene biotype information
gene_info <- getBM(
  attributes = c("external_gene_name", "gene_biotype"),
  filters = "external_gene_name",
  values = unique(tx2gene$SYMBOL),
  mart = mart
)

# Filter for protein-coding genes
coding_genes <- gene_info %>%
  filter(gene_biotype == "protein_coding") %>%
  pull(external_gene_name)

tx2gene_coding <- tx2gene %>%
  filter(SYMBOL %in% coding_genes)

cat("Transcripts before filtering:", nrow(tx2gene), "\n")
cat("Transcripts after filtering (protein-coding only):", nrow(tx2gene_coding), "\n")

# ----------------------------------------------------------------------------
# 4. Define Sample Information
# ----------------------------------------------------------------------------
# Sample files from 4 independent cachexia studies
sample_names <- c(
  # Control samples (n=20)
  'JCI(Joshua)_Cont1.tsv','JCI(Joshua)_Cont2.tsv','JCI(Joshua)_Cont3.tsv','JCI(Joshua)_Cont4.tsv',
  "JCI(Joshua)_Sham1.tsv","JCI(Joshua)_Sham2.tsv","JCI(Joshua)_Sham3.tsv","JCI(Joshua)_Sham4.tsv",
  "JNCI(Tseng)_Cont1.tsv","JNCI(Tseng)_Cont2.tsv","JNCI(Tseng)_Cont3.tsv",
  "JEM(Rupert)_Sham1.tsv", "JEM(Rupert)_Sham2.tsv", "JEM(Rupert)_Sham3.tsv",
  "Embo(sophia)_Cont1.tsv","Embo(sophia)_Cont2.tsv","Embo(sophia)_Cont3.tsv",
  "Embo(sophia)_Cont4.tsv","Embo(sophia)_Cont5.tsv","Embo(sophia)_Cont6.tsv",
  # Cachexia samples (n=19)
  "JCI(Joshua)_SC1.tsv","JCI(Joshua)_SC2.tsv","JCI(Joshua)_SC3.tsv","JCI(Joshua)_SC4.tsv",
  "JCI(Joshua)_SPC1.tsv","JCI(Joshua)_SPC2.tsv","JCI(Joshua)_SPC3.tsv","JCI(Joshua)_SPC4.tsv",
  "JNCI(Tseng)_C261.tsv","JNCI(Tseng)_C262.tsv","JNCI(Tseng)_C263.tsv",
  "JEM(Rupert)_KPC1.tsv","JEM(Rupert)_KPC2.tsv","JEM(Rupert)_KPC3.tsv","JEM(Rupert)_KPC4.tsv",
  "Embo(sophia)_CX1.tsv","Embo(sophia)_CX2.tsv","Embo(sophia)_CX3.tsv", "Embo(sophia)_CX4.tsv"
)

# Define file paths
files <- file.path("/Users/daehwankim/Desktop/Desktop_MacBook Air/sequencing data/Integrated data", sample_names)
names(files) <- sample_names

# Define batch information for batch effect correction
batch <- c(
  # Control batches
  rep("JCI_Joshua", 4), rep("JCI_Joshua2", 4), rep("JNCI_Tseng", 3),
  rep("JEM_Rupert", 3), rep("Embo_sophia", 6),
  # Cachexia batches
  rep("JCI_Joshua", 4), rep("JCI_Joshua2", 4), rep("JNCI_Tseng", 3),
  rep("JEM_Rupert", 4), rep("Embo_sophia", 4)
)

# Create sample table
sampleTable <- data.frame(
  condition = factor(rep(c("Control", "Cachexia"), times = c(20, 19))),
  batch = factor(batch)
)

# ----------------------------------------------------------------------------
# 5. Import Kallisto Quantification Data
# ----------------------------------------------------------------------------
txi.kallisto <- tximport(
  files, 
  type = 'kallisto', 
  tx2gene = tx2gene_coding, 
  ignoreAfterBar = TRUE, 
  ignoreTxVersion = TRUE
)

rownames(sampleTable) <- colnames(txi.kallisto$counts)

# ----------------------------------------------------------------------------
# 6. DESeq2 Analysis with Batch Correction
# ----------------------------------------------------------------------------
# Create DESeq2 dataset with batch + condition design
dds_batch <- DESeqDataSetFromTximport(
  txi.kallisto, 
  sampleTable,
  design = ~ batch + condition
)

# Filter low-count genes (keep genes with >= 10 counts in >= 19 samples)
smallestGroupSize <- 19
keep <- rowSums(counts(dds_batch) >= 10) >= smallestGroupSize
dds_batch <- dds_batch[keep, ]

cat("Genes before filtering:", nrow(txi.kallisto$counts), "\n")
cat("Genes after filtering:", nrow(dds_batch), "\n")

# Set Control as reference level
dds_batch$condition <- relevel(dds_batch$condition, ref = "Control")

# Run DESeq2
deseq2.dds <- DESeq(dds_batch)
deseq2.res <- results(deseq2.dds)
deseq2.res <- deseq2.res[order(rownames(deseq2.res)), ]

# Convert to data frame
res_df <- data.frame(deseq2.res)

# Save DESeq2 results
write.xlsx(res_df, file.path("/Users/daehwankim/Desktop/KIST_folder/Progress/Cachexia model/Serpina3_Paper/260119/DESeq2_results.xlsx"), rowNames = TRUE)


# ----------------------------------------------------------------------------
# 7. Figure 1B: PCA Plot (with ComBat batch correction for visualization)
# ----------------------------------------------------------------------------
# VST normalization
vsd <- vst(dds_batch, blind = FALSE)
vst_mat <- assay(vsd)

# ComBat batch correction for visualization
mod_combat <- model.matrix(~ condition, data = sampleTable)
batch_var <- sampleTable$batch

combat_mat <- ComBat(
  dat = vst_mat, 
  batch = batch_var, 
  mod = mod_combat, 
  par.prior = TRUE, 
  prior.plots = FALSE
)

# PCA calculation
pca_res <- prcomp(t(combat_mat))
percentVar <- round(100 * (pca_res$sdev^2 / sum(pca_res$sdev^2)), 1)

# Create PCA data frame
pcaData <- data.frame(
  PC1 = pca_res$x[, 1],
  PC2 = pca_res$x[, 2],
  condition = sampleTable$condition,
  batch = sampleTable$batch,
  name = colnames(dds_batch)
)

# Color scheme
pca_colors <- c("Control" = "#373838", "Cachexia" = "#20416C") 
pca_fill_colors <- c("Control" = "lightgray", "Cachexia" = "skyblue")

# Plot PCA (Figure 1B)
pca_plot <- ggplot(pcaData, aes(x = PC1, y = PC2, color = condition, fill = condition)) +
  geom_point(size = 2, alpha = 1, stroke = 0) +
  stat_ellipse(
    geom = "polygon", level = 0.95, alpha = 0.1,  
    linetype = "dashed", linewidth = 0.5, show.legend = FALSE
  ) +
  scale_color_manual(values = pca_colors) +
  scale_fill_manual(values = pca_fill_colors) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  labs(title = "PCA Plot") +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5)
  ) +
  coord_cartesian(xlim = c(-60, 100), ylim = c(-35, 35))

pca_plot <- pca_plot +
  theme(legend.position = "none")

print(pca_plot)
ggsave(file.path(OUTPUT_DIR, "Figure1B_PCA.pdf"), pca_plot, width = 8, height = 6)

# ----------------------------------------------------------------------------
# 8. Figure 1C: Volcano Plot
# ----------------------------------------------------------------------------
# Prepare data for volcano plot
plot_data <- res_df %>%
  as.data.frame() %>%
  rownames_to_column("gene_symbol") %>%
  filter(!is.na(padj)) %>%
  mutate(gene = case_when(
    log2FoldChange >= 1 & padj < 0.05 ~ "Upregulated",
    log2FoldChange <= -1 & padj < 0.05 ~ "Downregulated",
    TRUE ~ "Non significant"
  ))



table(plot_data$gene)

# Count DEGs
cat("\nDEG counts:\n")
print(table(plot_data$gene))

# Color scheme
volcano_colors <- c(
  "Downregulated" = '#1B3361',
  "Upregulated" = '#9F1D1F',
  "Non significant" = 'grey90'
)

# Highlight atrophy markers (validation of cachexia signature)
atrophy_genes_list <- c("Fbxo32", "Trim63","Eda2r")
atrophy_subset <- plot_data %>% filter(gene_symbol %in% atrophy_genes_list)

# Plot Volcano (Figure 1C)
volcano_plot <- ggplot(plot_data, aes(x = log2FoldChange, y = -log10(padj), col = gene)) +
  geom_point(
    shape = 16, stroke = 0, alpha = 0.8, size = 1, 
    position = position_jitter(width = 0.05, height = 0.05)
  ) +
  scale_color_manual(values = volcano_colors) +
  labs(x = "Log2 fold change", y = "-Log10P") +
  theme_minimal() +
  theme(
    axis.title.x = element_text(face = "bold", size = 10),
    axis.title.y = element_text(face = "bold", size = 10),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5), 
    panel.grid.minor = element_blank(), 
    panel.grid.major = element_blank(), 
    legend.position = "right"
  ) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", linewidth = 0.4, col = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", linewidth = 0.4, col = "black") +
  # Highlight atrophy markers
  geom_point(
    data = atrophy_subset, 
    aes(x = log2FoldChange, y = -log10(padj)),
    color = "black", fill = "#9F1D1F", size = 1.5, shape = 21, stroke = 0
  ) +
  geom_text_repel(
    data = atrophy_subset, 
    aes(label = gene_symbol),
    box.padding = 0.5, segment.color = 'black', segment.size = 0.4,
    size = 3, fontface = "bold.italic", color = "black",
    max.overlaps = Inf, force = 5, nudge_y = 5
  )

volcano_plot <- volcano_plot +
  theme(legend.position = "none")


print(volcano_plot)

combined_plot <- pca_plot + volcano_plot 
combined_plot <- pca_plot / volcano_plot

print(combined_plot)

ggsave(file.path(OUTPUT_DIR, "Figure1C_Volcano.pdf"), volcano_plot, width = 8, height = 6)

# ----------------------------------------------------------------------------
# 9. Figure 1D: K-means Clustering Heatmap
# ----------------------------------------------------------------------------
# Filter significant DEGs
sig_genes <- deseq2.res %>%
  as.data.frame() %>%
  rownames_to_column("gene") %>%
  filter(!is.na(padj), padj < 0.05, abs(log2FoldChange) >= 1) %>%
  pull(gene)

cat("\nSignificant DEGs:", length(sig_genes), "\n")

# Subset and scale expression matrix
mat_deg <- combat_mat[intersect(rownames(combat_mat), sig_genes), , drop = FALSE]
keep_sd <- apply(mat_deg, 1, sd) > 0
mat_deg <- mat_deg[keep_sd, , drop = FALSE]

# Z-score normalization
z <- t(scale(t(mat_deg)))

# Determine optimal k using elbow and silhouette methods
wss <- sapply(2:10, function(k) {
  kmeans(z, centers = k, nstart = 100, iter.max = 100)$tot.withinss
})

# K-means clustering (k = 5)
set.seed(1234)
k <- 5
km <- kmeans(z, centers = k, nstart = 200, iter.max = 100)

# Create cluster mapping table
cluster_tbl <- tibble(
  gene = names(km$cluster),
  cluster_num = km$cluster,
  cluster_label = paste0("GC", km$cluster)
) %>%
  arrange(cluster_num) %>%
  mutate(cluster_label = factor(cluster_label, levels = paste0("GC", 1:k)))

# Order expression matrix by cluster
z_ord <- z[cluster_tbl$gene, , drop = FALSE]

# Prepare column annotation
ann_col <- sampleTable %>% 
  dplyr::select(condition, batch) %>%
  mutate(
    condition = factor(condition, levels = c("Control", "Cachexia")),
    batch = factor(batch)
  ) %>%
  arrange(condition, batch)

z_ord <- z_ord[, rownames(ann_col), drop = FALSE]

# Prepare row annotation
ann_row <- data.frame(Cluster = cluster_tbl$cluster_label)
rownames(ann_row) <- cluster_tbl$gene

# Calculate gaps between clusters
row_gaps <- cumsum(table(cluster_tbl$cluster_label))
row_gaps <- row_gaps[-length(row_gaps)]

# Color scheme for heatmap
bk <- seq(-1, 1, by = 0.005)
cols <- colorRampPalette(c("#1a3664", "white", "red4"))(length(bk) - 1)

ann_colors <- list(
  condition = c(Control = "gray", Cachexia = "#1a3664"),
  batch = c(
    Embo_sophia = "lightgray", JCI_Joshua = "gray40",
    JCI_Joshua2 = "skyblue2", JEM_Rupert = "lightpink2", JNCI_Tseng = "red4"
  ),
  Cluster = c(GC1 = "gray", GC2 = "lightpink", GC3 = "red4", GC4 = "skyblue2", GC5 = "#1a3664")
)

# Plot heatmap (Figure 1D)
pdf(file.path(OUTPUT_DIR, "Figure1D_Heatmap.pdf"), width = 10, height = 12)
pheatmap(
  z_ord,           
  scale = "none",
  cluster_rows = FALSE,    
  cluster_cols = FALSE,
  color = cols, 
  breaks = bk,
  annotation_col = ann_col,
  annotation_row = ann_row,
  annotation_colors = ann_colors,
  show_colnames = FALSE,
  show_rownames = FALSE,
  gaps_row = row_gaps,
  main = "K-means Clustering of DEGs (GC1-GC5)",
  border_color = NA
)
dev.off()

# Save cluster assignments
write.xlsx(cluster_tbl, file.path(OUTPUT_DIR, "kmeans_clusters.xlsx"), rowNames = FALSE)

# ----------------------------------------------------------------------------
# 10. Figure 1E: GO Over-Representation Analysis (ORA)
# ----------------------------------------------------------------------------
# Run compareCluster for all gene clusters
ora_all <- compareCluster(
  gene ~ cluster_label,
  data = cluster_tbl,
  fun = "enrichGO",
  OrgDb = org.Mm.eg.db,
  keyType = "SYMBOL",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2
)

# Save full ORA results
df_ora <- as.data.frame(ora_all)
write.xlsx(df_ora, file.path(OUTPUT_DIR, "ORA_full_results.xlsx"))

# Select representative GO terms for visualization
picked <- tribble(
  ~Cluster, ~Description,
  "GC1", "axon guidance",
  "GC1", "neuron projection guidance",
  "GC1", "amide metabolic process",
  "GC2", "acute inflammatory response",
  "GC2", "leukocyte migration",
  "GC2", "regulation of inflammatory response",
  "GC3", "regulation of cellular catabolic process",
  "GC3", "regulation of protein ubiquitination",
  "GC3", "regulation of autophagy",
  "GC4", "muscle system process",
  "GC4", "muscle contraction",
  "GC4", "muscle cell differentiation",
  "GC5", "extracellular matrix organization",
  "GC5", "external encapsulating structure organization",
  "GC5", "extracellular structure organization"
)

# Filter ORA results for selected terms
cmpdf_sub <- df_ora %>%
  inner_join(picked, by = c("Cluster", "Description"))

target_cluster_order <- c("GC1", "GC2", "GC3", "GC4", "GC5")
cmpdf_sub$Cluster <- factor(cmpdf_sub$Cluster, levels = target_cluster_order)

real_desc_order <- picked$Description[picked$Description %in% cmpdf_sub$Description]
cmpdf_sub$Description <- factor(cmpdf_sub$Description, levels = rev(real_desc_order))

# Create subset object for plotting
cmp_sub <- ora_all
cmp_sub@compareClusterResult <- cmpdf_sub

# Plot ORA dot plot (Figure 1E)
ora_plot <- dotplot(
  cmp_sub, 
  x = "Cluster", 
  color = "p.adjust", 
  showCategory = nrow(cmpdf_sub),
  label_format = 60
) + 
  labs(title = "GO Enrichment Analysis by Gene Cluster") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 10, face = "bold", color = "black"),
    axis.text.y = element_text(size = 10, color = "black"),
    panel.grid.major = element_line(color = "grey90", linetype = "dashed")
  ) +
  scale_color_gradient(low = "red", high = "blue", trans = "log10", guide = guide_colorbar(reverse = TRUE))

print(ora_plot)
ggsave(file.path(OUTPUT_DIR, "Figure1E_ORA_dotplot.pdf"), ora_plot, width = 10, height = 8)

# ----------------------------------------------------------------------------
# 11. Save Session Info and Key Objects
# ----------------------------------------------------------------------------
# Save key objects for downstream analysis
saveRDS(res_df, file.path(OUTPUT_DIR, "DESeq2_results.rds"))
saveRDS(cluster_tbl, file.path(OUTPUT_DIR, "kmeans_cluster_table.rds"))
saveRDS(combat_mat, file.path(OUTPUT_DIR, "combat_corrected_matrix.rds"))

# Session info
writeLines(capture.output(sessionInfo()), file.path(OUTPUT_DIR, "session_info.txt"))

cat("\n============================================\n")
cat("Figure 1 analysis completed!\n")
cat("Output files saved to:", OUTPUT_DIR, "\n")
cat("============================================\n")
