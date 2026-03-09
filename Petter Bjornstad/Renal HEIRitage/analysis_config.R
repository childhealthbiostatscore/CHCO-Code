# ============================================================================
# Shared analysis configuration
# Source this file in any script that needs folders or celltype_groups.
# When adding a new analysis, update ONLY this file.
# ============================================================================

# ---------------------------------------------------------------------------
# Analysis folders: named vector mapping S3 folder names -> file suffixes
# ---------------------------------------------------------------------------
folders <- c(
  # Unadjusted
  "T2D_GLP_Y_vs_T2D_GLP_N" = "t2d_glpyn",
  "T2D_GLP_N_vs_HC" = "t2d_glpn_hc",
  # Adjusted for age + sex.x
  "T2D_GLP_Y_vs_T2D_GLP_N_adj_age_sex" = "t2d_glpyn_adj_age_sex",
  "T2D_GLP_N_vs_HC_adj_age_sex" = "t2d_glpn_hc_adj_age_sex"
)

# ---------------------------------------------------------------------------
# Cell type groupings (all resolution levels)
# ---------------------------------------------------------------------------
celltype_groups <- list(
  PT = c("PT-S1/S2", "PT-S3", "aPT"),
  aPT = "aPT",
  `PT-S1/S2` = "PT-S1/S2",
  `PT-S3` = "PT-S3",
  TAL = c("C-TAL-1", "C-TAL-2", "aTAL", "dTAL"),
  `C-TAL-1` = "C-TAL-1",
  `C-TAL-2` = "C-TAL-2",
  aTAL = "aTAL",
  dTAL = "dTAL",
  EC = c("EC-AVR", "EC-GC", "EC-PTC", "EC-AEA", "EC-LYM", "EC/VSMC", "EC-A"),
  `EC-AVR` = "EC-AVR",
  `EC-GC`  = "EC-GC",
  `EC-PTC` = "EC-PTC",
  `EC-AEA` = "EC-AEA",
  `EC-LYM` = "EC-LYM",
  `EC/VSMC` = "EC/VSMC",
  Immune = c("MAC", "MON", "cDC", "pDC", "CD4+ T", "CD8+ T", "B", "NK", "cycT"),
  Immune_Myeloid = c("MAC", "MON", "cDC", "pDC"),
  Immune_Lymphoid = c("CD4+ T", "CD8+ T", "B", "NK", "cycT"),
  PC = c("CCD-PC", "CNT-PC", "dCCD-PC", "M-PC", "tPC-IC"),
  IC = c("IC-A", "IC-B", "aIC"),
  DTL_ATL = c("DTL", "aDTL", "ATL"),
  DCT_CNT = c("DCT", "dDCT", "CNT"),
  VSMC_P_FIB = c("VSMC/P", "FIB"),
  POD = "POD",
  MC = "MC",
  PEC = "PEC",
  Schwann = "SchwannCells"
)

# High-resolution cell type subset (for DEG summary plots)
celltype_groups_high <- list(
  aPT = "aPT",
  `PT-S1/S2` = "PT-S1/S2",
  `PT-S3` = "PT-S3",
  `C-TAL-1` = "C-TAL-1",
  `C-TAL-2` = "C-TAL-2",
  aTAL = "aTAL",
  dTAL = "dTAL",
  `EC-AVR` = "EC-AVR",
  `EC-GC`  = "EC-GC",
  `EC-PTC` = "EC-PTC",
  `EC-AEA` = "EC-AEA",
  `EC-LYM` = "EC-LYM",
  `EC/VSMC` = "EC/VSMC",
  Immune_Myeloid = c("MAC", "MON", "cDC", "pDC"),
  Immune_Lymphoid = c("CD4+ T", "CD8+ T", "B", "NK", "cycT")
)

# Low-resolution cell type subset (used in reversal/pathway analyses)
celltype_groups_low <- list(
  PT = c("PT-S1/S2", "PT-S3", "aPT"),
  TAL = c("C-TAL-1", "C-TAL-2", "aTAL", "dTAL"),
  EC = c("EC-AVR", "EC-GC", "EC-PTC", "EC-AEA", "EC-LYM", "EC/VSMC", "EC-A"),
  Immune = c("MAC", "MON", "cDC", "pDC", "CD4+ T", "CD8+ T", "B", "NK", "cycT"),
  PC = c("CCD-PC", "CNT-PC", "dCCD-PC", "M-PC", "tPC-IC"),
  IC = c("IC-A", "IC-B", "aIC"),
  DTL_ATL = c("DTL", "aDTL", "ATL"),
  DCT_CNT = c("DCT", "dDCT", "CNT"),
  VSMC_P_FIB = c("VSMC/P", "FIB"),
  POD = "POD",
  MC = "MC",
  PEC = "PEC",
  Schwann = "SchwannCells"
)
