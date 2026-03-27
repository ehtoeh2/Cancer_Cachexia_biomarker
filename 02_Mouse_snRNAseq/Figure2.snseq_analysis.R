# ============================================================================
# Figure 2: Single-nucleus RNA-seq Analysis for Cancer Cachexia
# ============================================================================
# This script performs:
#   - QC & Doublet Removal      : per-sample filtering + scDblFinder
#   - Integration & Clustering  : SCTransform → PCA → Harmony → UMAP
#   - Cell Type Annotation      : marker-based cluster assignment
#   - Figure 2A : UMAP (annotated cell types)
#   - Figure 2B : Stacked bar plot (cell type proportions)
#   - Figure 2C : Heatmap (secretome DEGs in IIb / IIa/IIx / FAPs)
#   - Figure 2D : Dot plot (Serpina3n, Lox, Ctsl across cell types)
#   - Figure 2E : Feature plot (Tnfrsf12a, split by condition)
#   - Supp Fig 1B : Dot plot (canonical marker genes)
#
# Author : Dae-Hwan Kim, Jebeom Ko
# Paper  : "Diagnostic Platform for Cancer Cachexia: Multi-Omics Based Novel Biomarker Discovery"
#
# Public datasets used (snRNA-seq, Mus musculus, Gastrocnemius/Tibialis anterior):
#   - GSE272085     : C57BL/6, male, 8-12 weeks, Lewis lung carcinoma (LLC model)
#   - Zenodo.11090497 : Mixed background (C57BL/6 + FVB), male, KIC model
# ============================================================================

# ----------------------------------------------------------------------------
# 0. Load Required Packages
# ----------------------------------------------------------------------------
library(Seurat)
library(harmony)
library(scDblFinder)
library(glmGamPoi) # SCTransform v2 backend
library(openxlsx)
library(writexl)
library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)
library(patchwork)
library(RColorBrewer)
library(viridis)
library(clusterProfiler)
library(org.Mm.eg.db)
library(msigdbr)
library(fgsea)
library(enrichplot)
library(ComplexHeatmap)
library(circlize)
library(grid)

set.seed(1234)

options(future.globals.maxSize = 1e15) # allow large objects for SCTransform

# ----------------------------------------------------------------------------
# 1. Data Paths Configuration
# ----------------------------------------------------------------------------
# NOTE: Modify these paths according to your local environment.
#
# Datasets used (two independent snRNA-seq cachexia studies):
#   - GSE272085 (JCSM study)      : Con_filtered_feature_bc_matrix / cac_filtered_feature_bc_matrix
#   - Zenodo.11090497 (KIC model) : ctrl_feature_bc_matrix / ccx_feature_bc_matrix

DATA_DIR <- "/Users/daehwankim/Desktop/Seqencing_Practicing/Single seq analysis practice/Single Cell Cachexia(통합)"
OUTPUT_DIR <- "/Users/daehwankim/Desktop/KIST_folder/Progress/Cachexia model/Serpina3_Paper/260119"

# Paths for pre-processed secretome gene list
ALLGC_PATH <- file.path(OUTPUT_DIR, "allGC_secretion[260119].xlsx")

# Create output directory if needed
if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)

# ----------------------------------------------------------------------------
# 2. Load Raw 10X Data
# ----------------------------------------------------------------------------
control.data1 <- Read10X(data.dir = file.path(DATA_DIR, "JCSM/Con_filtered_feature_bc_matrix"))
control.data2 <- Read10X(data.dir = file.path(DATA_DIR, "Cell reports (Myod)/ctrl_feature_bc_matrix"))
cachexia.data1 <- Read10X(data.dir = file.path(DATA_DIR, "JCSM/cac_filtered_feature_bc_matrix"))
cachexia.data2 <- Read10X(data.dir = file.path(DATA_DIR, "Cell reports (Myod)/ccx_feature_bc_matrix"))

# ----------------------------------------------------------------------------
# 3. Create Seurat Objects
# ----------------------------------------------------------------------------
Control1 <- CreateSeuratObject(control.data1, project = "Control")
Control2 <- CreateSeuratObject(control.data2, project = "Control")
Cachexia1 <- CreateSeuratObject(cachexia.data1, project = "Cachexia")
Cachexia2 <- CreateSeuratObject(cachexia.data2, project = "Cachexia")

cat("Initial cell counts:\n")
cat("  Control1  :", ncol(Control1), "\n")
cat("  Control2  :", ncol(Control2), "\n")
cat("  Cachexia1 :", ncol(Cachexia1), "\n")
cat("  Cachexia2 :", ncol(Cachexia2), "\n")

# ----------------------------------------------------------------------------
# 4. QC: Compute Mitochondrial & Ribosomal Percentages
# ----------------------------------------------------------------------------
seurat_list <- list(Control1, Control2, Cachexia1, Cachexia2)

seurat_list <- lapply(seurat_list, function(obj) {
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^mt-")
  obj[["percent.ribo"]] <- PercentageFeatureSet(obj, pattern = "^Rps|^Rpl")
  obj
})

# Optional: visualize QC metrics before filtering
lapply(seurat_list, function(obj) {
  VlnPlot(obj,
    features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"),
    ncol = 4
  )
})

# ----------------------------------------------------------------------------
# 5. QC: Cell Filtering
# ----------------------------------------------------------------------------
seurat_list <- lapply(seurat_list, function(obj) {
  subset(obj,
    subset = nFeature_RNA > 200 &
      nFeature_RNA < 5000 &
      percent.mt < 5
  )
})

cat("\nCell counts after QC filtering:\n")
cat("  Control1  :", ncol(seurat_list[[1]]), "\n")
cat("  Control2  :", ncol(seurat_list[[2]]), "\n")
cat("  Cachexia1 :", ncol(seurat_list[[3]]), "\n")
cat("  Cachexia2 :", ncol(seurat_list[[4]]), "\n")

# ----------------------------------------------------------------------------
# 6. Doublet Removal with scDblFinder
# ----------------------------------------------------------------------------
run_scDbl <- function(sobj, dbr = NULL, seed = 1234) {
  # Convert to SingleCellExperiment, run doublet detection, return singlets only
  set.seed(seed)
  sce <- as.SingleCellExperiment(sobj)
  sce <- scDblFinder(sce, dbr = dbr)
  sobj$scDblFinder.class <- colData(sce)$scDblFinder.class
  sobj$scDblFinder.score <- colData(sce)$scDblFinder.score
  subset(sobj, subset = scDblFinder.class == "singlet")
}

seurat_list <- lapply(seurat_list, run_scDbl)

cat("\nCell counts after doublet removal:\n")
cat("  Control1  :", ncol(seurat_list[[1]]), "\n")
cat("  Control2  :", ncol(seurat_list[[2]]), "\n")
cat("  Cachexia1 :", ncol(seurat_list[[3]]), "\n")
cat("  Cachexia2 :", ncol(seurat_list[[4]]), "\n")

# ----------------------------------------------------------------------------
# 7. Merge Samples into One Seurat Object
# ----------------------------------------------------------------------------
Control_merged <- merge(seurat_list[[1]],
  y = seurat_list[[2]],
  add.cell.ids = c("Control1", "Control2")
)
Cachexia_merged <- merge(seurat_list[[3]],
  y = seurat_list[[4]],
  add.cell.ids = c("Cachexia1", "Cachexia2")
)

combined <- merge(Control_merged,
  y = Cachexia_merged,
  add.cell.ids = c("Control", "Cachexia")
)

# Free memory
rm(
  seurat_list, Control_merged, Cachexia_merged,
  control.data1, control.data2, cachexia.data1, cachexia.data2
)
gc()

# Set condition factor levels
combined$orig.ident <- factor(combined$orig.ident, levels = c("Control", "Cachexia"))

# ----------------------------------------------------------------------------
# 8. Normalization: SCTransform (v2)
# ----------------------------------------------------------------------------
DefaultAssay(combined) <- "RNA"

combined <- SCTransform(
  combined,
  vars.to.regress  = "percent.mt",
  vst.flavor       = "v2",
  conserve.memory  = TRUE,
  do.correct.umi   = FALSE,
  verbose          = TRUE
)

# ----------------------------------------------------------------------------
# 9. Dimensionality Reduction & Batch Correction (Harmony)
# ----------------------------------------------------------------------------
combined <- RunPCA(combined, npcs = 50)

# Visualize elbow plot to select appropriate number of PCs
ElbowPlot(combined, ndims = 50)

# Harmony integration (batch = orig.ident)
combined <- RunHarmony(combined, group.by.vars = "orig.ident", dims = 1:10)
combined <- FindNeighbors(combined, reduction = "harmony", dims = 1:10)
combined <- RunUMAP(combined, reduction = "harmony", dims = 1:10)
combined <- FindClusters(combined, resolution = 0.4)

# Quick check: unlabeled cluster UMAP
DimPlot(combined,
  reduction = "umap", label = TRUE, label.size = 4,
  pt.size = 0.1
) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  ) +
  labs(x = "UMAP 1", y = "UMAP 2")

# ----------------------------------------------------------------------------
# 10. Identify Cluster Marker Genes
# ----------------------------------------------------------------------------
combined <- PrepSCTFindMarkers(combined)

combined_marker_genes <- FindAllMarkers(
  combined,
  only.pos = TRUE,
  min.pct = 0.01,
  logfc.threshold = 1
)

write.xlsx(
  combined_marker_genes,
  file.path(OUTPUT_DIR, "cluster_allmarker.xlsx")
)

# ----------------------------------------------------------------------------
# 11. Cell Type Annotation
# ----------------------------------------------------------------------------
# Canonical marker genes used for annotation
marker_genes <- c(
  "Myh4", # IIb
  "Myh1", "Myh2", # IIa/IIx
  "Pdgfra", "Cd34", # FAPs
  "Pax7", "Calcr", # MuSCs
  "Pecam1", # Endothelial
  "Pdgfrb", # Pericyte
  "Chrne", # NMJ
  "Mrc1", # Macrophage
  "Col22a1" # MTJ
)

# Supplementary Figure 1B: canonical marker dot plot
DefaultAssay(combined) <- "SCT"
DotPlot(combined,
  features = marker_genes,
  cols = c("#FFFFCC", "blue3")
) + RotatedAxis()

# Cluster-to-cell-type mapping (update as needed based on marker expression)
Idents(combined) <- "seurat_clusters"

Annotation <- c(
  "0"  = "IIb",
  "1"  = "IIa/IIx",
  "2"  = "IIb",
  "3"  = "FAPs",
  "4"  = "IIa/IIx",
  "5"  = "FAPs",
  "6"  = "IIb",
  "7"  = "IIb",
  "8"  = "IIa/IIx",
  "9"  = "IIb",
  "10" = "Endothelial",
  "11" = "MTJ",
  "12" = "Pericyte",
  "13" = "Macrophage",
  "14" = "MuSCs",
  "15" = "NMJ",
  "16" = "IIb"
)

combined$Annotation <- unname(Annotation[as.character(Idents(combined))])
cat("Unannotated cells:", sum(is.na(combined$Annotation)), "\n")

# Set display order and convert to factor
cell_type_order <- c(
  "IIb", "IIa/IIx", "FAPs", "MuSCs",
  "Endothelial", "Pericyte", "NMJ", "Macrophage", "MTJ"
)
combined$Annotation <- factor(combined$Annotation, levels = cell_type_order)

# Switch active identity
Idents(combined) <- "Annotation"

# Save annotated object
saveRDS(combined, file.path(OUTPUT_DIR, "combined.rds"))
# combined <- readRDS(file.path(OUTPUT_DIR, "combined.rds"))

# Shared color palette for all figures
my_colors <- c(
  "IIb"         = "#3B6E91",
  "IIa/IIx"     = "#8AB6D6",
  "FAPs"        = "#6FA07F",
  "MuSCs"       = "#E6C229",
  "Endothelial" = "#C97586",
  "Pericyte"    = "#D1A6A8",
  "NMJ"         = "#8F80BA",
  "Macrophage"  = "#D98C5F",
  "MTJ"         = "#8894A8"
)

# ----------------------------------------------------------------------------
# 12. Figure 2A: UMAP (annotated, colored by cell type)
# ----------------------------------------------------------------------------
fig2A <- DimPlot(
  combined,
  reduction  = "umap",
  label      = FALSE,
  pt.size    = 0.4,
  order      = FALSE,
  cols       = my_colors
) +
  theme(
    legend.position = "right",
    axis.text       = element_blank(),
    axis.ticks      = element_blank(),
    axis.line       = element_blank()
  ) +
  labs(x = "UMAP 1", y = "UMAP 2")

print(fig2A)
ggsave(file.path(OUTPUT_DIR, "Figure2A_UMAP.pdf"), fig2A, width = 7, height = 5)

# ----------------------------------------------------------------------------
# 13. Figure 2B: Stacked Bar Plot (cell type proportions per condition)
# ----------------------------------------------------------------------------
prop_df <- combined@meta.data %>%
  mutate(
    group   = as.character(orig.ident),
    cluster = as.character(Annotation)
  ) %>%
  filter(!is.na(cluster)) %>%
  dplyr::count(group, cluster, name = "n") %>%
  group_by(group) %>%
  mutate(frac = n / sum(n)) %>%
  ungroup() %>%
  mutate(
    cluster = factor(cluster, levels = cell_type_order),
    group   = factor(group, levels = c("Control", "Cachexia"))
  )

fig2B <- ggplot(prop_df, aes(x = group, y = frac, fill = cluster)) +
  geom_col(width = 0.7, color = "black", linewidth = 0.3) +
  geom_text(
    aes(label = ifelse(frac > 0.05, as.character(cluster), "")),
    position  = position_stack(vjust = 0.5),
    color     = "black",
    size      = 2,
    fontface  = "bold"
  ) +
  scale_fill_manual(values = my_colors) +
  scale_y_continuous(labels = percent_format(accuracy = 1), expand = c(0, 0)) +
  labs(x = NULL, y = "Cell proportion (%)") +
  theme_classic(base_size = 14) +
  theme(
    axis.text = element_text(color = "black"),
    axis.line = element_line(linewidth = 0.3, color = "black"),
    axis.ticks = element_line(color = "black"),
    legend.position = "none",
    plot.margin = margin(10, 10, 10, 10)
  )

print(fig2B)
ggsave(file.path(OUTPUT_DIR, "Figure2B_CellProportion.pdf"), fig2B, width = 4, height = 6)

# ----------------------------------------------------------------------------
# 14. Figure 2C: Heatmap (secretome DEGs in IIb / IIa/IIx / FAPs)
# ----------------------------------------------------------------------------

# -- 14-1. Subset muscle & FAP clusters
sub3 <- subset(combined, subset = Annotation %in% c("IIb", "IIa/IIx", "FAPs"))
sub3$Group <- factor(as.character(sub3$orig.ident), levels = c("Control", "Cachexia"))
sub3$group_id <- paste(sub3$Annotation, sub3$Group, sep = "-")

# -- 14-2. Find upregulated DEGs in Cachexia vs. Control
DefaultAssay(sub3) <- "SCT"
sub3 <- PrepSCTFindMarkers(sub3)

markers_sub3 <- FindMarkers(
  sub3,
  assay = "SCT",
  ident.1 = "Cachexia",
  ident.2 = "Control",
  only.pos = TRUE,
  min.pct = 0.01,
  logfc.threshold = 1
)

markers_sub3_sig <- markers_sub3[markers_sub3$p_val_adj < 0.05, ]
markers_sub3_sig$genes <- rownames(markers_sub3_sig)

cat(
  "\nUpregulated DEGs (Cachexia vs. Control) in IIb/IIa/IIx/FAPs:",
  nrow(markers_sub3_sig), "\n"
)

# -- 14-3. Intersect with secretome candidate list
allGC <- read.xlsx(ALLGC_PATH)
gene_candidates <- allGC[allGC$gene %in% markers_sub3_sig$genes, ]

GC2_genes <- intersect(
  gene_candidates$gene[gene_candidates$cluster_label == "GC2"],
  rownames(sub3)
)
GC3_genes <- intersect(
  gene_candidates$gene[gene_candidates$cluster_label == "GC3"],
  rownames(sub3)
)
all_genes <- unique(c(GC2_genes, GC3_genes))

cat("Secretome candidate genes (GC2):", length(GC2_genes), "\n")
cat("Secretome candidate genes (GC3):", length(GC3_genes), "\n")

# -- 14-4. Compute average expression per group_id
avg_exp <- AverageExpression(
  sub3,
  features  = all_genes,
  group.by  = "group_id",
  assays    = "SCT",
  slot      = "data"
)$SCT

# Order columns: Control before Cachexia within each cell type
col_order <- c(
  "IIb-Control", "IIb-Cachexia",
  "IIa/IIx-Control", "IIa/IIx-Cachexia",
  "FAPs-Control", "FAPs-Cachexia"
)
valid_cols <- col_order[col_order %in% colnames(avg_exp)]
avg_exp <- avg_exp[, valid_cols, drop = FALSE]

# -- 14-5. Z-score scale genes (row-wise)
scaled_mat <- t(scale(t(avg_exp)))
scaled_mat[is.na(scaled_mat)] <- 0

# -- 14-6. Order genes by avg_log2FC (descending)
rank_tbl <- markers_sub3_sig
rank_tbl$gene <- rank_tbl$genes
ordered_genes <- rank_tbl$gene[order(rank_tbl$avg_log2FC, decreasing = TRUE)]
ordered_genes <- intersect(ordered_genes, rownames(scaled_mat))

scaled_mat <- scaled_mat[ordered_genes, , drop = FALSE]

# -- 14-7. Gene group annotation (GC2 / GC3)
gene_group_vec <- rep(NA_character_, nrow(scaled_mat))
names(gene_group_vec) <- rownames(scaled_mat)
gene_group_vec[rownames(scaled_mat) %in% GC2_genes] <- "GC2"
gene_group_vec[rownames(scaled_mat) %in% GC3_genes] <- "GC3"
gene_group_vec <- factor(gene_group_vec, levels = c("GC2", "GC3"))

# -- 14-8. Column split by cell type
col_split <- factor(sub("-.*", "", colnames(scaled_mat)),
  levels = c("IIa/IIx", "IIb", "FAPs")
)

# -- 14-9. Color function
col_fun <- colorRamp2(
  c(-2, 0, 2),
  c("#2166AC", "white", "#B2182B")
)

# -- 14-10. Draw heatmap
cell_w <- unit(13, "mm")
cell_h <- unit(7, "mm")

ht_fig2C <- Heatmap(
  scaled_mat,
  name             = "Z-score",
  col              = col_fun,
  row_split        = gene_group_vec,
  column_split     = col_split,
  cluster_rows     = FALSE,
  cluster_columns  = FALSE,
  column_gap       = unit(c(0, 4), "mm"),
  width            = cell_w * ncol(scaled_mat),
  height           = cell_h * nrow(scaled_mat),
  row_names_gp     = gpar(fontsize = 10, fontface = "bold.italic"),
  column_names_gp  = gpar(fontsize = 10, fontface = "bold")
)

pdf(file.path(OUTPUT_DIR, "Figure2C_Heatmap.pdf"), width = 8, height = 10)
draw(ht_fig2C)
dev.off()

# ----------------------------------------------------------------------------
# 15. Figure 2D: Dot Plot (Serpina3n, Lox, Ctsl across cell types)
# ----------------------------------------------------------------------------
DefaultAssay(combined) <- "SCT"

fig2D <- DotPlot(
  sub3,
  features  = c("Serpina3n", "Lox", "Ctsl"),
  dot.scale = 9,
  cols      = c("lightgrey", "#d7301f")
) +
  RotatedAxis() +
  theme(axis.title = element_blank())

print(fig2D)
ggsave(file.path(OUTPUT_DIR, "Figure2D_DotPlot.pdf"), fig2D, width = 5, height = 4)

# ----------------------------------------------------------------------------
# 16. Figure 2E: Feature Plot (Tnfrsf12a, split by condition)
# ----------------------------------------------------------------------------
fig2E <- FeaturePlot(
  combined,
  features     = "Tnfrsf12a",
  split.by     = "orig.ident",
  keep.scale   = "feature",
  order        = TRUE,
  pt.size      = 0.1,
  min.cutoff   = "q05",
  max.cutoff   = "q95",
  cols         = c("grey90", "#d7301f")
)

print(fig2E)
ggsave(file.path(OUTPUT_DIR, "Figure2E_FeaturePlot_Tnfrsf12a.pdf"), fig2E,
  width = 10, height = 5
)

# ----------------------------------------------------------------------------
# 17. Save Session Info
# ----------------------------------------------------------------------------
saveRDS(combined, file.path(OUTPUT_DIR, "combined_final.rds"))
writeLines(
  capture.output(sessionInfo()),
  file.path(OUTPUT_DIR, "session_info_Figure2.txt")
)

cat("\n============================================\n")
cat("Figure 2 analysis completed!\n")
cat("Output files saved to:", OUTPUT_DIR, "\n")
cat("============================================\n")
