# T1D adiposity analysis

# create clinical dataset
library(tidyverse)
library(arsenal)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(ggpubr)
library(Hmisc)
library(readxl)
library(purrr)
library(aws.s3)
library(jsonlite)

user <- Sys.info()[["user"]]

if (user == "choiyej") {
  root_path <- "/Users/choiyej/Library/CloudStorage/OneDrive-UW/Bjornstad/Biostatistics Core Shared Drive/"
  git_path <- "/Users/choiyej/GitHub/CHCO-Code/Petter Bjornstad/"
  keys <- fromJSON("/Users/choiyej/Library/CloudStorage/OneDrive-TheUniversityofColoradoDenver/Bjornstad Pyle Lab/keys.json")
} else if (user == "yejichoi") {
  root_path <- "/mmfs1/gscratch/togo/yejichoi/"
  git_path <- "/mmfs1/gscratch/togo/yejichoi/CHCO-Code/Petter Bjornstad"
  keys <- fromJSON("/mmfs1/home/yejichoi/keys.json")
} else if (user == "pylell") {
  root_path <- "/Users/pylell/Library/CloudStorage/OneDrive-SharedLibraries-UW/Bjornstad/Biostatistics Core Shared Drive"
  git_path <- "/Users/pylell/Documents/GitHub/CHCO-Code/Petter Bjornstad"
} else if (user == "hhampson") {
  root_path <- "/Volumes/Peds Endo"
} else {
  stop("Unknown user: please specify root path for this user.")
}

## Create an S3 client

Sys.setenv(
  "AWS_ACCESS_KEY_ID" = keys$MY_ACCESS_KEY,
  "AWS_SECRET_ACCESS_KEY" = keys$MY_SECRET_KEY,
  "AWS_DEFAULT_REGION" = "",
  "AWS_REGION" = "",
  "AWS_S3_ENDPOINT" = "s3.kopah.uw.edu"
)

# merge with Seurat
t1d_hc_dat <- s3readRDS(object = "data_clean/t1d_hc_clinical_data.csv", 
                        bucket = "t1d.adiposity",
                        region = "")

pb90_attempt <- s3readRDS(object = "data_clean/pb90_attempt_integrated_processed.rds",
                          bucket = "scrna",
                          region = "")

pb90_attempt_meta_subset <- pb90_attempt@meta.data %>%
  dplyr::select(1:7, record_id, kit_id, visit, source, celltype_rpca, celltype_harmony, KPMP_celltype) %>%
  dplyr::mutate(kit_id = gsub("kl", "KL", gsub("KI", "KL", kit_id)),
                visit = case_when(visit == "screening" ~ "baseline",
                                  T ~ visit)) %>%
  left_join(t1d_hc_dat) %>%
  dplyr::mutate(rownames = paste0(source, "_", barcode)) %>%
  filter(!is.na(uuid)) %>%
  filter(record_id != "CRC-55") # remove HC with IgAN

# create KPMP_celltype_general
celltype_groups <- list(
  PT = c("PT-S1/S2", "PT-S3", "aPT"),
  TAL = c("C-TAL-1", "C-TAL-2", "aTAL", "dTAL"),
  PC = c("CCD-PC", "CNT-PC", "dCCD-PC", "M-PC", "tPC-IC"),
  IC = c("IC-A", "IC-B", "aIC"),
  DTL_ATL = c("DTL", "aDTL", "ATL"),   # grouped thin limbs
  DCT_CNT = c("DCT", "dDCT", "CNT"),   # grouped distal tubule/connecting tubule
  EC = c("EC-AVR", "EC-GC", "EC-PTC", "EC-AEA", "EC-LYM", "EC/VSMC"), 
  Immune = c("MAC", "MON", "cDC", "pDC", "CD4+ T", "CD8+ T", "B", "NK", "cycT"),
  # Immune_Myeloid = c("MAC", "MON", "cDC", "pDC"),
  # Immune_Lymphoid = c("CD4+ T", "CD8+ T", "B", "NK", "cycT"),
  VSMC_P_FIB = c("VSMC/P", "FIB"),
  POD = "POD",
  MC = "MC",                         # mesangial cells
  PEC = "PEC",                       # parietal epithelial cells
  Schwann = "SchwannCells",
  Other = c("non-specific")          # catchall
)

map_celltype_to_general <- function(celltype, celltype_groups) {
  for (group_name in names(celltype_groups)) {
    if (celltype %in% celltype_groups[[group_name]]) {
      return(group_name)
    }
  }
  return("Other")  # For any celltype not found in groups
}

pb90_attempt_meta_subset$KPMP_celltype_general <- sapply(
  pb90_attempt_meta_subset$KPMP_celltype, 
  map_celltype_to_general, 
  celltype_groups = celltype_groups)

rownames(pb90_attempt_meta_subset) <- pb90_attempt_meta_subset$rownames

table((pb90_attempt_meta_subset %>%
  distinct(record_id, visit, .keep_all = T))$study) # confirmed 16 ATTEMPT BL T1D, 5 PANDA T1D, 28 CROC T1D, 12 CROC HC

pb90_attempt_subset <- pb90_attempt[, colnames(pb90_attempt) %in% rownames(pb90_attempt_meta_subset)]
pb90_attempt_subset@meta.data <- pb90_attempt_meta_subset

s3saveRDS(pb90_attempt_subset, 
          object = "data_clean/t1d_hc_scrna_w_clinical.rds", 
          bucket = "t1d.adiposity",
          region = "",
          multipart = T)

table((pb90_attempt_subset@meta.data %>%
         distinct(record_id, visit, .keep_all = T))$dxa_obesity,
      (pb90_attempt_subset@meta.data %>%
         distinct(record_id, visit, .keep_all = T))$group)

table((pb90_attempt_subset@meta.data %>%
         distinct(record_id, visit, .keep_all = T))$bmi_obesity,
      (pb90_attempt_subset@meta.data %>%
         distinct(record_id, visit, .keep_all = T))$group)

table((pb90_attempt_subset@meta.data %>%
         distinct(record_id, visit, .keep_all = T))$group)
