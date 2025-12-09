library(tidyverse)
library(readxl)
library(Seurat)
library(pedbp)

# specify user for paths
user <- Sys.info()[["user"]]
if (user == "laurapyle") {
  data_path <- '/Users/lpyle/Library/CloudStorage/OneDrive-UW/Bjornstad'
} else if (user == "lpyle") {
  data_path <- "/Users/lpyle/Library/CloudStorage/OneDrive-UW/Bjornstad"
} else if (user == "pylell") {
  data_path <- "/Users/lpyle/Library/CloudStorage/OneDrive-UW/Bjornstad"
} else {
  stop("Unknown user: please specify root path for this user.")
}

setwd(data_path)

# Import top proteins for MIC, MAC, etc. at baseline
top_mic_df <- read_excel("./Biostatistics Core Shared Drive/TODAY subaward/Results/Linear and Cox models/TODAY somalogic Cox models scaled baseline adjusted.xlsx", sheet = "MIC CPH")
de_genes_mic <- top_mic_df[top_mic_df$p.value <= 0.05, c("EntrezGeneID", "estimate")]
de_genes_mic <- setNames(de_genes_mic$estimate, de_genes_mic$EntrezGeneID)
top_mic <- top_mic_df %>%
  filter(adj.p.value <= 0.05) %>%
  slice_max(abs(log(estimate)), n = 5) %>%
  pull(AptName)
top_mac_df <- read_excel("./Biostatistics Core Shared Drive/TODAY subaward/Results/Linear and Cox models/TODAY somalogic Cox models scaled baseline adjusted.xlsx", sheet = "MAC CPH")
de_genes_mac <- top_mac_df[top_mac_df$p.value <= 0.05, c("EntrezGeneID", "estimate")]
de_genes_mac <- setNames(de_genes_mac$estimate, de_genes_mac$EntrezGeneID)
top_mac <- top_mac_df %>%
  filter(adj.p.value <= 0.05) %>%
  slice_max(abs(log(estimate)), n = 5) %>%
  pull(AptName)
top_mic.or.mac_df <- read_excel("./Biostatistics Core Shared Drive/TODAY subaward/Results/Linear and Cox models/TODAY somalogic Cox models scaled baseline adjusted.xlsx", sheet = "MIC.OR.MAC CPH")
de_genes_mic.or.mac <- top_mic.or.mac_df[top_mic.or.mac_df$p.value <= 0.05, c("EntrezGeneID", "estimate")]
de_genes_mic.or.mac <- setNames(de_genes_mic.or.mac$estimate, de_genes_mic.or.mac$EntrezGeneID)
top_mic.or.mac <- top_mic.or.mac_df %>%
  filter(adj.p.value <= 0.05) %>%
  slice_max(abs(log(estimate)), n = 5) %>%
  pull(AptName)
top_hyp_df <- read_excel("./Biostatistics Core Shared Drive/TODAY subaward/Results/Linear and Cox models/TODAY somalogic Cox models scaled baseline adjusted.xlsx", sheet = "HYP CPH")
de_genes_hyp <- top_hyp_df[top_hyp_df$p.value <= 0.05, c("EntrezGeneID", "estimate")]
de_genes_hyp <- setNames(de_genes_hyp$estimate, de_genes_hyp$EntrezGeneID)
top_hyp <- top_hyp_df %>%
  filter(adj.p.value <= 0.05) %>%
  slice_max(abs(log(estimate)), n = 5) %>%
  pull(AptName)
top_rapid_df <- read_excel("./Biostatistics Core Shared Drive/TODAY subaward/Results/Linear and Cox models/TODAY somalogic Cox models scaled baseline adjusted.xlsx", sheet = "RAPID CPH")
de_genes_rapid <- top_rapid_df[top_rapid_df$p.value <= 0.05, c("EntrezGeneID", "estimate")]
de_genes_rapid <- setNames(de_genes_rapid$estimate, de_genes_rapid$EntrezGeneID)
top_rapid <- top_rapid_df %>%
  filter(adj.p.value <= 0.05) %>%
  slice_max(abs(log(estimate)), n = 5) %>%
  pull(AptName)
top_neuro_df <- read_excel("./Biostatistics Core Shared Drive/TODAY subaward/Results/Linear and Cox models/TODAY somalogic Cox models scaled baseline adjusted.xlsx", sheet = "NEURO CPH")
de_genes_neuro <- top_neuro_df[top_neuro_df$p.value <= 0.05, c("EntrezGeneID", "estimate")]
de_genes_neuro <- setNames(de_genes_neuro$estimate, de_genes_neuro$EntrezGeneID)
top_neuro <- top_neuro_df %>%
  filter(adj.p.value <= 0.05) %>%
  slice_max(abs(log(estimate)), n = 5) %>%
  pull(AptName)
top_retino_df <- read_excel("./Biostatistics Core Shared Drive/TODAY subaward/Results/Linear and Cox models/TODAY somalogic Cox models scaled baseline adjusted.xlsx", sheet = "RETINO CPH")
de_genes_retino <- top_retino_df[top_retino_df$p.value <= 0.05, c("EntrezGeneID", "estimate")]
de_genes_retino <- setNames(de_genes_retino$estimate, de_genes_retino$EntrezGeneID)
top_retino <- top_retino_df %>%
  filter(adj.p.value <= 0.05) %>%
  slice_max(abs(log(estimate)), n = 5) %>%
  pull(AptName)
top_glyc_df <- read_excel("./Biostatistics Core Shared Drive/TODAY subaward/Results/Linear and Cox models/TODAY somalogic Cox models scaled baseline adjusted.xlsx", sheet = "GLYC CPH")
de_genes_glyc <- top_glyc_df[top_glyc_df$p.value <= 0.05, c("EntrezGeneID", "estimate")]
de_genes_glyc <- setNames(de_genes_glyc$estimate, de_genes_glyc$EntrezGeneID)
top_glyc <- top_glyc_df %>%
  filter(adj.p.value <= 0.05) %>%
  slice_max(abs(log(estimate)), n = 5) %>%
  pull(AptName)
top_glyc_a1c_df <- read_excel("./Biostatistics Core Shared Drive/TODAY subaward/Results/Linear and Cox models/TODAY somalogic Cox models scaled baseline adjusted.xlsx", sheet = "GLYC with A1c CPH")
de_genes_glyc_a1c <- top_glyc_a1c_df[top_glyc_a1c_df$p.value <= 0.05, c("EntrezGeneID", "estimate")]
de_genes_glyc_a1c <- setNames(de_genes_glyc_a1c$estimate, de_genes_glyc_a1c$EntrezGeneID)
top_glyc_a1c <- top_glyc_a1c_df %>%
  filter(adj.p.value <= 0.05) %>%
  slice_max(abs(log(estimate)), n = 5) %>%
  pull(AptName)
# Import 10 year results
top_mac_df_10 <- read_excel("./Biostatistics Core Shared Drive/TODAY subaward/Results/Linear and Cox models/TODAY somalogic limma yr10 adjusted.xlsx", sheet = "MAC_moderated_FDR")
de_genes_mac_10 <- top_mac_df_10[top_mac_df_10$P.Value <= 0.05, c("EntrezGeneID", "logFC")]
de_genes_mac_10 <- setNames(de_genes_mac_10$logFC, de_genes_mac_10$EntrezGeneID)
top_mic_df_10 <- read_excel("./Biostatistics Core Shared Drive/TODAY subaward/Results/Linear and Cox models/TODAY somalogic limma yr10 adjusted.xlsx", sheet = "MIC_moderated_FDR")
de_genes_mic_10 <- top_mic_df_10[top_mic_df_10$P.Value <= 0.05, c("EntrezGeneID", "logFC")]
de_genes_mic_10 <- setNames(de_genes_mic_10$logFC, de_genes_mic_10$EntrezGeneID)
top_hyp_df_10 <- read_excel("./Biostatistics Core Shared Drive/TODAY subaward/Results/Linear and Cox models/TODAY somalogic limma yr10 adjusted.xlsx", sheet = "hyperfilt_moderated_FDR")
de_genes_hyp_10 <- top_hyp_df_10[top_hyp_df_10$P.Value <= 0.05, c("EntrezGeneID", "logFC")]
de_genes_hyp_10 <- setNames(de_genes_hyp_10$logFC, de_genes_hyp_10$EntrezGeneID)
top_rapid_df_10 <- read_excel("./Biostatistics Core Shared Drive/TODAY subaward/Results/Linear and Cox models/TODAY somalogic limma yr10 adjusted.xlsx", sheet = "rapid_moderated_FDR")
de_genes_rapid_10 <- top_rapid_df_10[top_rapid_df_10$P.Value <= 0.05, c("EntrezGeneID", "logFC")]
de_genes_rapid_10 <- setNames(de_genes_rapid_10$logFC, de_genes_rapid_10$EntrezGeneID)
top_neuro_df_10 <- read_excel("./Biostatistics Core Shared Drive/TODAY subaward/Results/Linear and Cox models/TODAY somalogic limma yr10 adjusted.xlsx", sheet = "neuro_moderated_FDR")
de_genes_neuro_10 <- top_neuro_df_10[top_neuro_df_10$P.Value <= 0.05, c("EntrezGeneID", "logFC")]
de_genes_neuro_10 <- setNames(de_genes_neuro_10$logFC, de_genes_neuro_10$EntrezGeneID)
top_glyc_df_10 <- read_excel("./Biostatistics Core Shared Drive/TODAY subaward/Results/Linear and Cox models/TODAY somalogic limma yr10 adjusted.xlsx", sheet = "glyc_moderated_FDR")
de_genes_glyc_10 <- top_glyc_df_10[top_glyc_df_10$P.Value <= 0.05, c("EntrezGeneID", "logFC")]
de_genes_glyc_10 <- setNames(de_genes_glyc_10$logFC, de_genes_glyc_10$EntrezGeneID)
top_glyc_a1c_df_10 <- read_excel("./Biostatistics Core Shared Drive/TODAY subaward/Results/Linear and Cox models/TODAY somalogic limma yr10 adjusted.xlsx", sheet = "glyc_with_a1c_moderated_FDR")
de_genes_glyc_a1c_10 <- top_glyc_a1c_df_10[top_glyc_a1c_df_10$P.Value <= 0.05, c("EntrezGeneID", "logFC")]
de_genes_glyc_a1c_10 <- setNames(de_genes_glyc_a1c_10$logFC, de_genes_glyc_a1c_10$EntrezGeneID)
# Import and clean data
df <- read.csv("./Biostatistics Core Shared Drive/Data Harmonization/Data Clean/soma_harmonized_dataset.csv", na.strings = "")
df <- df %>%
  filter(
    study %in% c("IMPROVE", "RENAL-HEIR","RENAL-HEIRitage"),
    !grepl("IT2D", co_enroll_id), participation_status == "Participated"
  ) %>%
  group_by(record_id, visit) %>%
  summarise(across(where(is.character), ~ last(na.omit(.x))),
    across(where(is.factor), ~ last(na.omit(.x))),
    across(where(is.numeric), ~ mean(.x, na.rm = T)),
    .groups = "drop"
  ) %>%
  mutate_all(~ ifelse(is.nan(.), NA, .))
# Import proteomics data for RH and IMPROVE
load("./Biostatistics Core Shared Drive/Data Harmonization/Combined SomaScan/analytes.Rdata")
soma <- df %>% select(record_id, visit, contains("seq."))
# Transform
soma[, 3:ncol(soma)] <- lapply(soma[, 3:ncol(soma)], log)
# BP percentiles and HTN
df$bp_age <- df$age * 12
df$bp_age[df$bp_age >= 19 * 12] <- 19 * 12 - 1e-8
df$male <- ifelse(df$sex == "Male", 1, 0)
bps <- p_bp(
  q_sbp = df$sbp, q_dbp = df$dbp, age = df$bp_age,
  male = df$male , height = df$height
)
df$sbp_perc <- bps$sbp_percentile
df$sbp_perc[df$age >= 19] <- NA
df$dbp_perc <- bps$dbp_percentile
df$dbp_perc[df$age >= 19] <- NA
df$htn <- df$sbp >= 130 | df$dbp >= 80 | df$sbp_perc >= 0.95 | df$dbp_perc >= 0.95
df$htn[is.na(df$htn)] <- F
df$htn[is.na(df$sbp) & is.na(df$dbp) & is.na(df$sbp_perc) & is.na(df$dbp_perc)] <- NA
df$htn <- factor(df$htn, levels = c(F, T), labels = c("HTN-", "HTN+"))
# Hyperfiltration
df$hyp <- factor(df$eGFR_fas_cr >= 135, levels = c(F, T), labels = c("eGFR < 135", "eGFR >= 135"))
# UACR
df$elevated_uacr <- factor(df$acr_u >= 30, labels = c("UACR < 30", "UACR >= 30"))
# Merge
df <- df %>% select(-contains("seq."))
df <- full_join(df, soma, by = c("record_id", "visit"))
# only baseline
df <- df %>% filter(visit == "baseline")

# Import Olink data and combine
olink_map <- read.csv("./Biostatistics Core Shared Drive/Olink Data/Data_Clean/olink_id_map.csv")
load("./Biostatistics Core Shared Drive/Olink Data/improve_olink_plasma.Rdata")
load("./Biostatistics Core Shared Drive/Olink Data/improve_olink_urine.Rdata")
load("./Biostatistics Core Shared Drive/Olink Data/rh_olink_plasma.Rdata")
load("./Biostatistics Core Shared Drive/Olink Data/rh_olink_urine.Rdata")
olink_plasma <- rbind(improve_olink_plasma, rh_olink_plasma)
olink_urine <- rbind(improve_olink_urine, rh_olink_urine)
# Add visit column, format IDs
olink_plasma$visit <- sapply(str_split(olink_plasma$record_id, "_"), "[", 3)
olink_plasma <- olink_plasma %>%
  mutate(
    visit = case_when(
      visit == "BL" ~ "baseline",
      visit == "3M" ~ "3_months_post_surgery",
      visit == "12M" ~ "12_months_post_surgery",
      is.na(visit) ~ "baseline"
    ),
    record_id = sub("_BL|_3M|_12M", "", record_id)
  )
olink_urine$visit <- sapply(str_split(olink_urine$record_id, "_"), "[", 3)
olink_urine <- olink_urine %>%
  mutate(
    visit = case_when(
      visit == "BL" ~ "baseline",
      visit == "3M" ~ "3_months_post_surgery",
      visit == "12M" ~ "12_months_post_surgery",
      is.na(visit) ~ "baseline"
    ),
    record_id = sub("_BL|_3M|_12M", "", record_id)
  )
# Baseline visit only 
df <- df %>% filter(visit == "baseline")
# Combine
#plasma <- left_join(df, olink_plasma, by = c("record_id", "visit"))
#urine <- left_join(df, olink_urine, by = c("record_id", "visit"))
# Limit df to those with all data
#ids <- intersect(soma$record_id, olink_plasma$record_id)
#df <- df %>% filter(record_id %in% ids)
# Save
# save(df, plasma, urine, analytes, olink_map,
#   top_mac, top_mic, top_mic.or.mac, top_hyp, top_rapid, top_glyc,
#   top_mac_df, top_mic_df, top_mic.or.mac_df, top_hyp_df, top_rapid_df, top_glyc_df,
#   de_genes_mac, de_genes_mic, de_genes_mic.or.mac, de_genes_hyp, de_genes_rapid, de_genes_glyc,
#   de_genes_mac_10, de_genes_mic_10, de_genes_hyp_10, de_genes_rapid_10, de_genes_glyc_10,
#   file = "./Proteomics and DKD/Data_Cleaned/analysis_dataset.RData"
# )
# Save without plasma or urine - Olink data
save(df, analytes, olink_map,
  top_mac, top_mic, top_mic.or.mac, top_hyp, top_rapid, top_glyc,
  top_mac_df, top_mic_df, top_mic.or.mac_df, top_hyp_df, top_rapid_df, top_glyc_df,
  de_genes_mac, de_genes_mic, de_genes_mic.or.mac, de_genes_hyp, de_genes_rapid, de_genes_glyc,
  de_genes_mac_10, de_genes_mic_10, de_genes_hyp_10, de_genes_rapid_10, de_genes_glyc_10,
  file = "./Biostatistics Core Shared Drive/Proteomics and DKD/Data_Cleaned/analysis_dataset.RData"
)

# # Read in scRNA object
# so <- readRDS("/Volumes/Shared/Shared Projects/Laura/Peds Endo/Petter Bjornstad/scRNA/Data_Clean/seurat_data_no_computations.RDS")
# # Limit to RH/IMPROVE baseline visit
# so <- so[, !grepl("_12M", so$michigan_id)]
# # Limit to those with SOMA and Olink
# so$michigan_id <- sub("_BL", "", so$michigan_id)
# so <- so[, so$michigan_id %in% ids]
# # Combined groups
# so$diabetes <- sub("i", "", so$T2D_HC_Phil)
# so$SGLT2i <- factor(so$T2D_HC_Phil,
#   levels = c("HC", "T2D", "T2Di"),
#   labels = c("SGLT2i-", "SGLT2i-", "SGLT2i+")
# )
# so$elevated_uacr <- df$elevated_uacr[match(so$michigan_id, df$record_id)]
# so$htn <- df$htn[match(so$michigan_id, df$record_id)]
# so$hyp <- df$hyp[match(so$michigan_id, df$record_id)]
# # Normalize and scale
# so <- NormalizeData(so)
# so <- ScaleData(so, features = rownames(so))
# # PCA
# so <- RunPCA(so, features = VariableFeatures(object = so))
# # Cluster cells
# so <- FindNeighbors(so)
# so <- FindClusters(so)
# # Perform UMAP
# so <- RunUMAP(so, dims = 1:20)
# # General cell types as identifiers
# so$generaltype <- sub("_.*", "", so$LR_clusters)
# Idents(so) <- so$generaltype
# # Save
# save(so,
#   file = "/Volumes/Shared/Shared Projects/Laura/Peds Endo/Petter Bjornstad/Proteomics and DKD/Data_Cleaned/analysis_seurat_object.RData"
# )
