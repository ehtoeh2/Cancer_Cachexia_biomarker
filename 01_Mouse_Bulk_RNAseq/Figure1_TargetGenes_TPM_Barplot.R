# ============================================================================
# Figure 1: Integrated Bulk RNA-seq Analysis
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
# 6. Extract TPM values for Target Genes
# ----------------------------------------------------------------------------
target_genes <- c("Trim63", "Fbxo32", "Eda2r")

# Abundance matrix contains TPM values (from tximport of kallisto quantifiers)
tpm_mat <- txi.kallisto$abundance

# Check if target genes are present
genes_present <- intersect(target_genes, rownames(tpm_mat))
cat("Genes found in kallisto data:", paste(genes_present, collapse = ", "), "\n")

if(length(genes_present) == 0) {
  stop("None of the target genes were found in the dataset.")
}

tpm_subset <- tpm_mat[genes_present, , drop = FALSE]

# Reshape data for plotting
tpm_df <- as.data.frame(tpm_subset) %>%
  rownames_to_column("Gene") %>%
  pivot_longer(cols = -Gene, names_to = "Sample", values_to = "TPM")

# Add condition information
sample_conditions <- data.frame(Sample = rownames(sampleTable), condition = sampleTable$condition)
tpm_df <- tpm_df %>% left_join(sample_conditions, by = "Sample")

# Ensure factor levels for ordering in plot
tpm_df$condition <- factor(tpm_df$condition, levels = c("Control", "Cachexia"))
tpm_df$Gene <- factor(tpm_df$Gene, levels = target_genes[target_genes %in% genes_present])

# ----------------------------------------------------------------------------
# 7. Generate Bar Plot
# ----------------------------------------------------------------------------
# Calculate summary statistics (mean and Standard Error of the Mean)
summary_df <- tpm_df %>%
  group_by(Gene, condition) %>%
  summarise(
    mean_TPM = mean(TPM),
    se_TPM = sd(TPM) / sqrt(n()),
    .groups = 'drop'
  )

# Define colors (matching original script aesthetics)
fill_colors <- c("Control" = "lightgray", "Cachexia" = "skyblue")
point_colors <- c("Control" = "#373838", "Cachexia" = "#20416C")

# Create Plot
tpm_plot <- ggplot(summary_df, aes(x = condition, y = mean_TPM, fill = condition)) +
  geom_bar(stat = "identity", position = position_dodge(), color = "black", alpha = 0.8, width = 0.7) +
  geom_errorbar(aes(ymin = mean_TPM - se_TPM, ymax = mean_TPM + se_TPM), 
                width = 0.2, position = position_dodge(0.9)) +
  # Add individual sample points (jittered)
  geom_jitter(data = tpm_df, aes(x = condition, y = TPM, color = condition), 
              width = 0.2, alpha = 0.8, size = 1.5) +
  scale_fill_manual(values = fill_colors) +
  scale_color_manual(values = point_colors) +
  facet_wrap(~ Gene, scales = "free_y") +
  labs(title = "TPM Expression of Target Genes", x = "", y = "TPM") +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    axis.text.x = element_text(size = 12, face = "bold", color = "black", angle = 45, hjust = 1, vjust = 1),
    axis.text.y = element_text(size = 10, color = "black"),
    axis.title.y = element_text(size = 12, face = "bold"),
    strip.text = element_text(size = 12, face = "bold.italic"),
    strip.background = element_rect(fill = "white", color = "black"),
    legend.position = "none",
    panel.grid.major.x = element_blank()
  )

view(tpm_df)
write.xlsx(tpm_df, "/Users/daehwankim/Desktop/KIST_folder/2025/Progress/Cachexia model/Serpina3_Paper/Nature communications/tpm.xlsx")

# ----------------------------------------------------------------------------
# 8. Save Plot and Session Info
# ----------------------------------------------------------------------------
print(tpm_plot)

# Define Figure Output Directory
FIGURE_DIR <- "/Users/daehwankim/Desktop/KIST_folder/2025/Progress/Cachexia model/Serpina3_Paper/Nature communications/Figureset_natcom"
if (!dir.exists(FIGURE_DIR)) dir.create(FIGURE_DIR, recursive = TRUE)

# Save PDF
ggsave(file.path(FIGURE_DIR, "Figure1_TargetGenes_TPM_Barplot.pdf"), tpm_plot, width = 8, height = 4)
# Save PNG
ggsave(file.path(FIGURE_DIR, "Figure1_TargetGenes_TPM_Barplot.png"), tpm_plot, width = 8, height = 4, dpi = 300)

cat("\n============================================\n")
cat("TPM Barplot generation completed!\n")
cat("Output files saved to:", OUTPUT_DIR, "\n")
cat("============================================\n")
