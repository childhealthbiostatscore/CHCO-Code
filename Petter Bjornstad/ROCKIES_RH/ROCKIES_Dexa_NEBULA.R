#DEXA Analysis
library(scran)
library(future)
library(future.apply)
library(tidyverse)
library(colorspace)
library(patchwork)
library(ggdendro)
library(cowplot)
library(ggpubr)
library(rstatix)
library(arsenal)
library(Biobase)
library(msigdbr)
library(kableExtra)
library(knitr)
library(REDCapR)
library(data.table)
library(emmeans)
library(NMF)
library(pheatmap)
library(UpSetR)
library(enrichR)
library(WriteXLS)
library(SAVER)
library(readxl)
library(limma)
library(edgeR)
library(BiocGenerics)
library(GSEABase)
library(slingshot)
library(SingleCellExperiment)
library(MAST)
library(muscat)
library(scater)
library(Seurat)
library(jsonlite)
library(dplyr)
library(glmmTMB)
library(reshape2)
library(broom.mixed)
library(nebula)
library(doParallel)




harmonized_data <- read.csv("C:/Users/netio/OneDrive - UW/Laura Pyle's files - Biostatistics Core Shared Drive/Data Harmonization/Data Clean/harmonized_dataset.csv", na = '')


dat <- harmonized_data %>%
  arrange(date_of_screen) %>% 
  dplyr::summarise(across(where(negate(is.numeric)), ~ ifelse(all(is.na(.x)), NA_character_, first(na.omit(.x)))),
                   across(where(is.numeric), ~ ifelse(all(is.na(.x)), NA_real_, first(na.omit(.x)))),
                   .by = c(record_id, visit))


dat2 <- dat %>% dplyr::select(record_id, visit, epic_sglti2_1, starts_with('dexa_')) %>% 
  filter(!is.na(epic_sglti2_1))


dat2 <- dat2 %>% filter(!is.na(dexa_body_fat))






dex_var <- c("dexa_ag_ratio",  "dexa_body_fat",            "dexa_bone_mineral_density", "dexa_est_vat",             
             "dexa_fat_kg",               "dexa_lean_kg",              "dexa_lean_mass",            "dexa_trunk_kg",            
            "dexa_trunk_mass"  )




gene_list <- c(
      "LEPR", "ADIPOR1", "ADIPOR2",
      "PPARG", "PPARA", "PPARD",
      "SREBF1", "SREBF2",
      "STAT3", "JAK2", "SOCS3",
      "TNF", "IL6", "IL1B", "CCL2",
      "NFKB1", "RELA", "IKBKB",
      "TLR2", "TLR4", "MYD88",
      "NLRP3", "CASP1", "IL18",
      "FASN", "ACACA", "SCD",
      "CPT1A", "CPT2", "ACOX1",
      "LDLR", "SCARB1", "ABCA1",
      "APOE", "APOB", "LPL",
      "DGAT1", "DGAT2", "PLIN1",
      "TGFB1", "TGFB2", "SMAD2", "SMAD3",
      "COL1A1", "COL3A1", "COL4A1",
      "FN1", "CTGF", "ACTA2",
      "MMP2", "MMP9", "TIMP1",
      "NOX1", "NOX2", "NOX4",
      "SOD1", "SOD2", "CAT", "GPX1",
      "NRF2", "KEAP1", "NQO1",
      "HMOX1", "GCLC", "GCLM",
      "HMOX1", "GCLC", "GCLM",
      "IRS1", "IRS2", "INSR",
      "GLUT4", "GLUT2",
      "SOCS1", "SOCS3", "PTPN1",
      "PRKAA1", "PRKAA2",
      "HSPA5", "EIF2AK3", "ERN1", "ATF6", # UPR sensors
      "XBP1", "ATF4", "DDIT3",
      "EDEM1", "DNAJB9", "HYOU1"
    )



#cell-type specific analyses 


 
"PT" = c(
                        cell_type,
                        # Fatty acid uptake and metabolism
                        "CD36", "FABP1", "FABP3", "SLC27A2", "SLC27A4",
                        # Glucose/lipid interaction
                        "PPARGC1A", "FOXO1", "SIRT1", "AMPK",
                        # PT-specific transporters affected by lipids
                        "SLC5A2", "SLC9A3", "SLC34A1"
                      )
              "POD" = c(
                        # Podocyte lipotoxicity
                        "NPHS1", "NPHS2", "PODXL", "SYNPO", "CD2AP",
                        # Cholesterol efflux
                        "ABCA1", "ABCG1", "APOE",
                      "CHOP", "BIP", "XBP1s"

)
"Endothelial" = c(
  # Endothelial dysfunction
  "VCAM1", "ICAM1", "SELE", "SELP",
  # NO synthesis
  "NOS3", "GUCY1A3", "GUCY1B3",
  # Angiogenesis
  "VEGFA", "VEGFR2", "ANGPT2"
)
"Immune" = c(
  # Macrophage polarization
  "CD68", "CD86", "CD163", "CD206",
  # Cytokines
  "IL10", "IL12", "TGFB1",
  # Chemokines
  "CCL2", "CCL3", "CCL5", "CXCL10"
)
"TAL" = c(
  # Metabolic adaptation
  "PPARGC1A", "NRF1", "TFAM",
  # Transport regulation
  "SLC12A1", "KCNJ1", "ATP1A1")

"DCT" = c(
  # Calcium handling
  "SLC12A3", "TRPM6", "CALB1",
  # Metabolic genes
  "HIF1A", "EPAS1"
)







#LC vs. T2D (no SLGT2)














#T2D SGLT2 Analysis 










