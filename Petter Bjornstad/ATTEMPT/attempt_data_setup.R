library(Seurat)
library(future)
library(tidyverse)
library(colorspace)
library(patchwork)
library(ggdendro)
library(cowplot)
library(ggpubr)
library(rstatix)
library(arsenal)
library(Biobase)
library(ReactomeGSA)
library(GSEABase)
library(msigdbr)
library(kableExtra)
library(knitr)
library(slingshot); library(SingleCellExperiment)
library(REDCapR)
library(edgeR)
library(data.table)
library(MAST)
library(limma)      # For linear modeling of microarray and RNA-seq data
library(muscat)

user <- Sys.info()[["user"]]

if (user == "choiyej") { # local version
  root_path <- "/Users/choiyej/Library/CloudStorage/OneDrive-UW/Bjornstad/Biostatistics Core Shared Drive"
  git_path <- "/Users/choiyej/GitHub/CHCO-Code/Petter Bjornstad"
  keys <- fromJSON("/Users/choiyej/Library/CloudStorage/OneDrive-TheUniversityofColoradoDenver/Bjornstad Pyle Lab/keys.json")
} else if (user == "yejichoi") { # hyak version
  root_path <- "/mmfs1/gscratch/togo/yejichoi"
  git_path <- "/mmfs1/gscratch/togo/yejichoi/CHCO-Code/Petter Bjornstad"
  keys <- fromJSON("/mmfs1/home/yejichoi/keys.json")
} else if (user == "pylell") {
  root_path <- "/Users/pylell/Library/CloudStorage/OneDrive-SharedLibraries-UW/Bjornstad/Biostatistics Core Shared Drive"
  git_path <- "/Users/pylell/Documents/GitHub/CHCO-Code/Petter Bjornstad"
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

# clinical data clean & save
# Load dictionary function and file
source(file.path(git_path, "/Resources/YC/R Functions/label_harmonized_function.R"))
source(file.path(git_path, "/CROCODILE/crocodile_functions.R"))

token_dat <- read.csv(file.path(root_path, "Data Harmonization/api_tokens.csv"))
attempt_token <- token_dat$Token[token_dat$Study == "ATTEMPT"]
attempt_dat <- redcap_read(redcap_uri = "https://redcap.ucdenver.edu/api/",
                           token = attempt_token)
attempt_grp <- read.csv(file.path(root_path, "ATTEMPT/Data Raw/ATTEMPT_unblinding.csv"))

race_names <- c(
  "American Indian or Alaskan Native",
  "Asian",
  "Hawaiian or Pacific Islander",
  "Black or African American",
  "White",
  "Unknown",
  "Other"
)
ethnicity_names <- c(
  "Hispanic",
  "Non-Hispanic",
  "Unknown")

dat <- attempt_dat$data %>%
  group_by(subject_id) %>%
  fill(everything(), .direction = "down") %>%
  ungroup() %>% rowwise() %>%
  dplyr::mutate(visit = case_when(redcap_event_name == "visit_2_arm_1" ~ "PRE",
                                  redcap_event_name == "visit_3_arm_1" ~ "POST",
                                  T ~ "Screening"),
                age = case_when(is.na(bx_date) ~ age_consent,
                                T ~ as.integer((bx_date - dob) / dyears(1))),
                race = case_when(sum(c_across(starts_with("race___")), na.rm = TRUE) == 0 ~ "Unknown",
                                 sum(c_across(starts_with("race___")), na.rm = TRUE) > 1 ~ "More than one race",
                                 race___1 == 1 ~ race_names[1],
                                 race___2 == 1 ~ race_names[2],
                                 race___3 == 1 ~ race_names[3],
                                 race___4 == 1 ~ race_names[4],
                                 race___5 == 1 ~ race_names[5],
                                 race___6 == 1 ~ race_names[6],
                                 race___7 == 1 ~ race_names[7]),
                ethnicity = case_when(sum(c_across(starts_with("ethnicity___")), na.rm = TRUE) == 0 ~ "Unknown",
                                      sum(c_across(starts_with("ethnicity___")), na.rm = TRUE) > 1 ~ "More than one ethnicity",
                                      ethnicity___1 == 1 ~ ethnicity_names[1],
                                      ethnicity___2 == 1 ~ ethnicity_names[2],
                                      ethnicity___3 == 1 ~ ethnicity_names[3]),
                diabetes_duration = case_when(visit == "Screening" ~ as.integer((consent_date - t1d_date) / dyears(1)),
                                              T ~ as.integer((bx_date - t1d_date) / dyears(1)))) %>%
  filter(visit != "Screening") %>%
  left_join(attempt_grp)

# manually adding participants missing in REDCap for now...
dat[nrow(dat) + 1, "subject_id"] <- 30186
dat[nrow(dat), "visit"] <- "PRE"
dat[nrow(dat), "treatment"] <- "Placebo"
dat[nrow(dat) + 1, "subject_id"] <- 30186
dat[nrow(dat), "visit"] <- "POST"
dat[nrow(dat), "treatment"] <- "Placebo"

dat$visit <- factor(dat$visit, levels = c("PRE", "POST"))

saveRDS(dat, file = file.path(root_path, "ATTEMPT/Data Clean/attempt_clinical_data.RDS"))

# scRNA
plan(multisession, workers = 16)
options(future.globals.maxSize=2e9)
so_attempt <- s3readRDS("Kidney transcriptomics/Single cell RNA seq/PB_attempt_harmony_rpca_Sept2024.RDS", "scrna", region = "")
so_attempt_meta <- so_attempt@meta.data %>%
  mutate(subject_id = Subject.ID,
         visit = case_when(Visit == "BL" ~ "PRE", 
                           Visit == "4M" ~ "POST"))
so_attempt_meta <- left_join(so_attempt_meta, dat, by = c("subject_id", "visit"))
rownames(so_attempt_meta) <- so_attempt_meta$barcode
so_attempt <- AddMetaData(so_attempt, so_attempt_meta)
saveRDS(so_attempt, file = file.path(root_path, "ATTEMPT/Data Clean/so_attempt.RDS"))

# PT cells
so_attempt$celltype_pt <- ifelse(grepl("PT-", so_attempt$celltype),
                                 "PT", as.character(so_attempt$celltype))
so_attempt_pt <- subset(so_attempt, celltype_pt == "PT" & celltype != "PT_lowQuality")
options(future.globals.maxSize = 3000 * 1024^3)
so_attempt_pt <- NormalizeData(so_attempt_pt)
so_attempt_pt <- ScaleData(so_attempt_pt)
ElbowPlot(so_attempt_pt)
so_attempt_pt <- RunPCA(so_attempt_pt, ncomponents = 10, features = VariableFeatures(object = so_attempt_pt))
so_attempt_pt <- FindNeighbors(so_attempt_pt)
so_attempt_pt <- FindClusters(so_attempt_pt)
so_attempt_pt <- RunUMAP(so_attempt_pt, dims = 1:30, reduction.key = "UMAP_")
gc()
so_attempt_pt$visit <- factor(so_attempt_pt$visit, levels = c("PRE", "POST"))
so_attempt_pt$treatment <- factor(so_attempt_pt$treatment, levels = c("Placebo", "Dapagliflozin 5mg"),
                                  labels = c("Placebo", "Dapagliflozin"))

saveRDS(so_attempt_pt, file = file.path(root_path, "ATTEMPT/Data Clean/attempt_so_pt.RDS"))


# TAL cells
so_attempt$celltype_tal <- ifelse(grepl("TAL-", so_attempt$celltype),
                                 "TAL", as.character(so_attempt$celltype))
so_attempt_tal <- subset(so_attempt, celltype_tal == "TAL" & celltype != "TAL_highUMI")
options(future.globals.maxSize = 3000 * 1024^3)
plan("multicore", workers = 4) 
so_attempt_tal <- NormalizeData(so_attempt_tal)
so_attempt_tal <- ScaleData(so_attempt_tal)
ElbowPlot(so_attempt_tal)
so_attempt_tal <- RunPCA(so_attempt_tal, ncomponents = 10, features = VariableFeatures(object = so_attempt_tal))
so_attempt_tal <- FindNeighbors(so_attempt_tal)
so_attempt_tal <- FindClusters(so_attempt_tal)
so_attempt_tal <- RunUMAP(so_attempt_tal, dims = 1:30, reduction.key = "UMAP_")
gc()
so_attempt_tal$visit <- factor(so_attempt_tal$visit, levels = c("PRE", "POST"))
so_attempt_tal$treatment <- factor(so_attempt_tal$treatment, levels = c("Placebo", "Dapagliflozin 5mg"),
                                  labels = c("Placebo", "Dapagliflozin"))

saveRDS(so_attempt_tal, file = file.path(root_path, "ATTEMPT/Data Clean/attempt_so_tal.RDS"))

# POD cells
so_attempt$celltype_pod <- ifelse(grepl("POD", so_attempt$celltype),
                                  "POD", as.character(so_attempt$celltype))
so_attempt_pod <- subset(so_attempt, celltype_pod == "POD")
options(future.globals.maxSize = 3000 * 1024^3)
plan("multicore", workers = 4) 
so_attempt_pod <- NormalizeData(so_attempt_pod)
so_attempt_pod <- ScaleData(so_attempt_pod)
ElbowPlot(so_attempt_pod)
so_attempt_pod <- RunPCA(so_attempt_pod, ncomponents = 10, features = VariableFeatures(object = so_attempt_pod))
so_attempt_pod <- FindNeighbors(so_attempt_pod)
so_attempt_pod <- FindClusters(so_attempt_pod)
so_attempt_pod <- RunUMAP(so_attempt_pod, dims = 1:30, reduction.key = "UMAP_")
gc()
so_attempt_pod$visit <- factor(so_attempt_pod$visit, levels = c("PRE", "POST"))
so_attempt_pod$treatment <- factor(so_attempt_pod$treatment, levels = c("Placebo", "Dapagliflozin 5mg"),
                                   labels = c("Placebo", "Dapagliflozin"))

saveRDS(so_attempt_pod, file = file.path(root_path, "ATTEMPT/Data Clean/attempt_so_pod.RDS"))

# EC cells
so_attempt$celltype_ec <- ifelse(grepl("EC-", so_attempt$celltype),
                                  "EC", as.character(so_attempt$celltype))
so_attempt_ec <- subset(so_attempt, celltype_ec == "EC")
options(future.globals.maxSize = 3000 * 1024^3)
plan("multicore", workers = 4) 
so_attempt_ec <- NormalizeData(so_attempt_ec)
so_attempt_ec <- ScaleData(so_attempt_ec)
ElbowPlot(so_attempt_ec)
so_attempt_ec <- RunPCA(so_attempt_ec, ncomponents = 10, features = VariableFeatures(object = so_attempt_ec))
so_attempt_ec <- FindNeighbors(so_attempt_ec)
so_attempt_ec <- FindClusters(so_attempt_ec)
so_attempt_ec <- RunUMAP(so_attempt_ec, dims = 1:30, reduction.key = "UMAP_")
gc()
so_attempt_ec$visit <- factor(so_attempt_ec$visit, levels = c("PRE", "POST"))
so_attempt_ec$treatment <- factor(so_attempt_ec$treatment, levels = c("Placebo", "Dapagliflozin 5mg"),
                                   labels = c("Placebo", "Dapagliflozin"))

saveRDS(so_attempt_ec, file = file.path(root_path, "ATTEMPT/Data Clean/attempt_so_ec.RDS"))

# PEC cells
so_attempt$celltype_pec <- ifelse(grepl("PEC", so_attempt$celltype),
                                 "PEC", as.character(so_attempt$celltype))
so_attempt_pec <- subset(so_attempt, celltype_pec == "PEC")
options(future.globals.maxSize = 3000 * 1024^3)
plan("multicore", workers = 4) 
so_attempt_pec <- NormalizeData(so_attempt_pec)
so_attempt_pec <- ScaleData(so_attempt_pec)
ElbowPlot(so_attempt_pec)
so_attempt_pec <- RunPCA(so_attempt_pec, ncomponents = 10, features = VariableFeatures(object = so_attempt_pec))
so_attempt_pec <- FindNeighbors(so_attempt_pec)
so_attempt_pec <- FindClusters(so_attempt_pec)
so_attempt_pec <- RunUMAP(so_attempt_pec, dims = 1:30, reduction.key = "UMAP_")
gc()
so_attempt_pec$visit <- factor(so_attempt_pec$visit, levels = c("PRE", "POST"))
so_attempt_pec$treatment <- factor(so_attempt_pec$treatment, levels = c("Placebo", "Dapagliflozin 5mg"),
                                  labels = c("Placebo", "Dapagliflozin"))

saveRDS(so_attempt_pec, file = file.path(root_path, "ATTEMPT/Data Clean/attempt_so_pec.RDS"))

# immune cells
so_attempt$celltype_immune <- ifelse(grepl("MAC|MON|B|T|NKT/NKC", so_attempt$celltype),
                                  "IMMUNE", as.character(so_attempt$celltype))
so_attempt_immune <- subset(so_attempt, celltype_immune == "IMMUNE")
so_attempt_immune <- NormalizeData(so_attempt_immune)
so_attempt_immune <- ScaleData(so_attempt_immune)
ElbowPlot(so_attempt_immune)
so_attempt_immune <- RunPCA(so_attempt_immune, ncomponents = 10, features = VariableFeatures(object = so_attempt_immune))
so_attempt_immune <- FindNeighbors(so_attempt_immune)
so_attempt_immune <- FindClusters(so_attempt_immune)
so_attempt_immune <- RunUMAP(so_attempt_immune, dims = 1:30, reduction.key = "UMAP_")
gc()
so_attempt_immune$visit <- factor(so_attempt_immune$visit, levels = c("PRE", "POST"))
so_attempt_immune$treatment <- factor(so_attempt_immune$treatment, levels = c("Placebo", "Dapagliflozin 5mg"),
                                   labels = c("Placebo", "Dapagliflozin"))

saveRDS(so_attempt_immune, file = file.path(root_path, "ATTEMPT/Data Clean/attempt_so_immune.RDS"))

# DTL cells
so_attempt$celltype_dtl <- ifelse(grepl("DTL-", so_attempt$celltype),
                                 "DTL", as.character(so_attempt$celltype))
so_attempt_dtl <- subset(so_attempt, celltype_dtl == "DTL")
options(future.globals.maxSize = 3000 * 1024^3)
so_attempt_dtl <- NormalizeData(so_attempt_dtl)
so_attempt_dtl <- ScaleData(so_attempt_dtl)
ElbowPlot(so_attempt_dtl)
so_attempt_dtl <- RunPCA(so_attempt_dtl, ncomponents = 10, features = VariableFeatures(object = so_attempt_dtl))
so_attempt_dtl <- FindNeighbors(so_attempt_dtl)
so_attempt_dtl <- FindClusters(so_attempt_dtl)
so_attempt_dtl <- RunUMAP(so_attempt_dtl, dims = 1:30, reduction.key = "UMAP_")
gc()
so_attempt_dtl$visit <- factor(so_attempt_dtl$visit, levels = c("PRE", "POST"))
so_attempt_dtl$treatment <- factor(so_attempt_dtl$treatment, levels = c("Placebo", "Dapagliflozin 5mg"),
                                  labels = c("Placebo", "Dapagliflozin"))

saveRDS(so_attempt_dtl, file = file.path(root_path, "ATTEMPT/Data Clean/attempt_so_dtl.RDS"))

# ATL cells
so_attempt$celltype_atl <- ifelse(grepl("ATL", so_attempt$celltype),
                                  "ATL", as.character(so_attempt$celltype))
so_attempt_atl <- subset(so_attempt, celltype_atl == "ATL")
options(future.globals.maxSize = 3000 * 1024^3)
so_attempt_atl <- NormalizeData(so_attempt_atl)
so_attempt_atl <- ScaleData(so_attempt_atl)
ElbowPlot(so_attempt_atl)
so_attempt_atl <- RunPCA(so_attempt_atl, ncomponents = 10, features = VariableFeatures(object = so_attempt_atl))
so_attempt_atl <- FindNeighbors(so_attempt_atl)
so_attempt_atl <- FindClusters(so_attempt_atl)
so_attempt_atl <- RunUMAP(so_attempt_atl, dims = 1:30, reduction.key = "UMAP_")
gc()
so_attempt_atl$visit <- factor(so_attempt_atl$visit, levels = c("PRE", "POST"))
so_attempt_atl$treatment <- factor(so_attempt_atl$treatment, levels = c("Placebo", "Dapagliflozin 5mg"),
                                   labels = c("Placebo", "Dapagliflozin"))

saveRDS(so_attempt_atl, file = file.path(root_path, "ATTEMPT/Data Clean/attempt_so_atl.RDS"))

# DCT cells
so_attempt$celltype_dct <- ifelse(grepl("DCT", so_attempt$celltype),
                                  "DCT", as.character(so_attempt$celltype))
so_attempt_dct <- subset(so_attempt, celltype_dct == "DCT")
options(future.globals.maxSize = 3000 * 1024^3)
so_attempt_dct <- NormalizeData(so_attempt_dct)
so_attempt_dct <- ScaleData(so_attempt_dct)
ElbowPlot(so_attempt_dct)
so_attempt_dct <- RunPCA(so_attempt_dct, ncomponents = 10, features = VariableFeatures(object = so_attempt_dct))
so_attempt_dct <- FindNeighbors(so_attempt_dct)
so_attempt_dct <- FindClusters(so_attempt_dct)
so_attempt_dct <- RunUMAP(so_attempt_dct, dims = 1:30, reduction.key = "UMAP_")
gc()
so_attempt_dct$visit <- factor(so_attempt_dct$visit, levels = c("PRE", "POST"))
so_attempt_dct$treatment <- factor(so_attempt_dct$treatment, levels = c("Placebo", "Dapagliflozin 5mg"),
                                   labels = c("Placebo", "Dapagliflozin"))

saveRDS(so_attempt_dct, file = file.path(root_path, "ATTEMPT/Data Clean/attempt_so_dct.RDS"))

# CNT cells
so_attempt$celltype_cnt <- ifelse(grepl("CNT", so_attempt$celltype),
                                  "CNT", as.character(so_attempt$celltype))
so_attempt_cnt <- subset(so_attempt, celltype_cnt == "CNT")
options(future.globals.maxSize = 3000 * 1024^3)
so_attempt_cnt <- NormalizeData(so_attempt_cnt)
so_attempt_cnt <- ScaleData(so_attempt_cnt)
ElbowPlot(so_attempt_cnt)
so_attempt_cnt <- RunPCA(so_attempt_cnt, ncomponents = 10, features = VariableFeatures(object = so_attempt_cnt))
so_attempt_cnt <- FindNeighbors(so_attempt_cnt)
so_attempt_cnt <- FindClusters(so_attempt_cnt)
so_attempt_cnt <- RunUMAP(so_attempt_cnt, dims = 1:30, reduction.key = "UMAP_")
gc()
so_attempt_cnt$visit <- factor(so_attempt_cnt$visit, levels = c("PRE", "POST"))
so_attempt_cnt$treatment <- factor(so_attempt_cnt$treatment, levels = c("Placebo", "Dapagliflozin 5mg"),
                                   labels = c("Placebo", "Dapagliflozin"))

saveRDS(so_attempt_cnt, file = file.path(root_path, "ATTEMPT/Data Clean/attempt_so_cnt.RDS"))

# PC cells
so_attempt$celltype_pc <- ifelse(grepl("PC-", so_attempt$celltype),
                                  "PC", as.character(so_attempt$celltype))
so_attempt_pc <- subset(so_attempt, celltype_pc == "PC")
options(future.globals.maxSize = 3000 * 1024^3)
so_attempt_pc <- NormalizeData(so_attempt_pc)
so_attempt_pc <- ScaleData(so_attempt_pc)
ElbowPlot(so_attempt_pc)
so_attempt_pc <- RunPCA(so_attempt_pc, ncomponents = 10, features = VariableFeatures(object = so_attempt_pc))
so_attempt_pc <- FindNeighbors(so_attempt_pc)
so_attempt_pc <- FindClusters(so_attempt_pc)
so_attempt_pc <- RunUMAP(so_attempt_pc, dims = 1:30, reduction.key = "UMAP_")
gc()
so_attempt_pc$visit <- factor(so_attempt_pc$visit, levels = c("PRE", "POST"))
so_attempt_pc$treatment <- factor(so_attempt_pc$treatment, levels = c("Placebo", "Dapagliflozin 5mg"),
                                   labels = c("Placebo", "Dapagliflozin"))

saveRDS(so_attempt_pc, file = file.path(root_path, "ATTEMPT/Data Clean/attempt_so_pc.RDS"))

# IC cells
so_attempt$celltype_ic <- ifelse(grepl("IC-", so_attempt$celltype),
                                 "IC", as.character(so_attempt$celltype))
so_attempt_ic <- subset(so_attempt, celltype_ic == "IC")
options(future.globals.maxSize = 3000 * 1024^3)
so_attempt_ic <- NormalizeData(so_attempt_ic)
so_attempt_ic <- ScaleData(so_attempt_ic)
ElbowPlot(so_attempt_ic)
so_attempt_ic <- RunPCA(so_attempt_ic, ncomponents = 10, features = VariableFeatures(object = so_attempt_ic))
so_attempt_ic <- FindNeighbors(so_attempt_ic)
so_attempt_ic <- FindClusters(so_attempt_ic)
so_attempt_ic <- RunUMAP(so_attempt_ic, dims = 1:30, reduction.key = "UMAP_")
gc()
so_attempt_ic$visit <- factor(so_attempt_ic$visit, levels = c("PRE", "POST"))
so_attempt_ic$treatment <- factor(so_attempt_ic$treatment, levels = c("Placebo", "Dapagliflozin 5mg"),
                                  labels = c("Placebo", "Dapagliflozin"))

saveRDS(so_attempt_ic, file = file.path(root_path, "ATTEMPT/Data Clean/attempt_so_ic.RDS"))

# MC cells
so_attempt$celltype_mc <- ifelse(grepl("MC-", so_attempt$celltype),
                                 "MC", as.character(so_attempt$celltype))
so_attempt_mc <- subset(so_attempt, celltype_mc == "MC")
options(future.globals.maxSize = 3000 * 1024^3)
so_attempt_mc <- NormalizeData(so_attempt_mc)
so_attempt_mc <- ScaleData(so_attempt_mc)
ElbowPlot(so_attempt_mc)
so_attempt_mc <- RunPCA(so_attempt_mc, ncomponents = 10, features = VariableFeatures(object = so_attempt_mc))
so_attempt_mc <- FindNeighbors(so_attempt_mc)
so_attempt_mc <- FindClusters(so_attempt_mc)
so_attempt_mc <- RunUMAP(so_attempt_mc, dims = 1:30, reduction.key = "UMAP_")
gc()
so_attempt_mc$visit <- factor(so_attempt_mc$visit, levels = c("PRE", "POST"))
so_attempt_mc$treatment <- factor(so_attempt_mc$treatment, levels = c("Placebo", "Dapagliflozin 5mg"),
                                  labels = c("Placebo", "Dapagliflozin"))

saveRDS(so_attempt_mc, file = file.path(root_path, "ATTEMPT/Data Clean/attempt_so_mc.RDS"))

# VSMC/P cells
so_attempt$celltype_vsmc_p <- ifelse(grepl("VSMC/P", so_attempt$celltype),
                                 "VSMC/P", as.character(so_attempt$celltype))
so_attempt_vsmc_p <- subset(so_attempt, celltype_vsmc_p == "VSMC/P")
options(future.globals.maxSize = 3000 * 1024^3)
so_attempt_vsmc_p <- NormalizeData(so_attempt_vsmc_p)
so_attempt_vsmc_p <- ScaleData(so_attempt_vsmc_p)
ElbowPlot(so_attempt_vsmc_p)
so_attempt_vsmc_p <- RunPCA(so_attempt_vsmc_p, ncomponents = 10, features = VariableFeatures(object = so_attempt_vsmc_p))
so_attempt_vsmc_p <- FindNeighbors(so_attempt_vsmc_p)
so_attempt_vsmc_p <- FindClusters(so_attempt_vsmc_p)
so_attempt_vsmc_p <- RunUMAP(so_attempt_vsmc_p, dims = 1:30, reduction.key = "UMAP_")
gc()
so_attempt_vsmc_p$visit <- factor(so_attempt_vsmc_p$visit, levels = c("PRE", "POST"))
so_attempt_vsmc_p$treatment <- factor(so_attempt_vsmc_p$treatment, levels = c("Placebo", "Dapagliflozin 5mg"),
                                  labels = c("Placebo", "Dapagliflozin"))

saveRDS(so_attempt_vsmc_p, file = file.path(root_path, "ATTEMPT/Data Clean/attempt_so_vsmc_p.RDS"))

# FIB cells
so_attempt$celltype_fib <- ifelse(grepl("FIB", so_attempt$celltype),
                                 "FIB", as.character(so_attempt$celltype))
so_attempt_fib <- subset(so_attempt, celltype_fib == "FIB")
options(future.globals.maxSize = 3000 * 1024^3)
so_attempt_fib <- NormalizeData(so_attempt_fib)
so_attempt_fib <- ScaleData(so_attempt_fib)
ElbowPlot(so_attempt_fib)
so_attempt_fib <- RunPCA(so_attempt_fib, ncomponents = 10, features = VariableFeatures(object = so_attempt_fib))
so_attempt_fib <- FindNeighbors(so_attempt_fib)
so_attempt_fib <- FindClusters(so_attempt_fib)
so_attempt_fib <- RunUMAP(so_attempt_fib, dims = 1:30, reduction.key = "UMAP_")
gc()
so_attempt_fib$visit <- factor(so_attempt_fib$visit, levels = c("PRE", "POST"))
so_attempt_fib$treatment <- factor(so_attempt_fib$treatment, levels = c("Placebo", "Dapagliflozin 5mg"),
                                  labels = c("Placebo", "Dapagliflozin"))

saveRDS(so_attempt_fib, file = file.path(root_path, "ATTEMPT/Data Clean/attempt_so_fib.RDS"))

# Fetch all data upfront
gene_names <- rownames(so_attempt_pt)
all_data <- FetchData(so_attempt_pt, vars = c("subject_id", "treatment", "visit", "ident", gene_names)) %>%
  group_by(subject_id, treatment, visit) %>%
  dplyr::summarise(across(where(is.numeric), ~ mean(.x > 0, na.rm = TRUE)), .groups = 'drop')

all_data <- all_data %>%
  group_by(subject_id) %>%
  filter(n() > 1) %>%
  ungroup() %>%
  mutate(treatment = case_when(treatment == "Placebo" ~ "Placebo",
                               T ~ "Dapagliflozin"))
gc()

all_data$treatment <- factor(all_data$treatment)
all_data$treatment <- relevel(all_data$treatment, ref = "Placebo")
saveRDS(all_data, file = file.path(root_path, "ATTEMPT/Data Clean/attempt_pb.RDS"))