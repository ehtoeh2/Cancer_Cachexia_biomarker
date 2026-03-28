# ============================================================================
# Figure 4: In vivo cancer cachexia model validates sex-independent 
# circulating biomarkers detectability
# ============================================================================
# This script performs:
#   - Figure 4H: Volcano plot

# Author: Dae-Hwan Kim, Jebeom Ko
# Paper: "Diagnostic Platform for Cancer Cachexia: Multi-Omics Based Novel Biomarker Discovery"
#
# Dataset used (Mouse serum LC-MS/MS proteomics):
#   - In-house generated mass spectrometry data (C57BL/6J, male, 8 weeks, LLC model)
#     Quadriceps muscle, serum; Cachexia vs. Control comparison
# ============================================================================

# ----------------------------------------------------------------------------
# 1. Figure 4H: Proteomics Volcano plot
# ----------------------------------------------------------------------------

# Required input files:
#   - Mouse serum proteomics data (Mass spectrometry results)
#   - Per-sample proteomics abundance data

Mouse_protein <- read.xlsx(file.path(DATA_DIR, "Mouse_serum_proteomics.xlsx"))

# See Figure 4B volcano plot
