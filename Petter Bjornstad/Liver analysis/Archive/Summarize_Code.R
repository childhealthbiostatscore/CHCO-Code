#Summary Code
#Directories
dir.results <- c("/Users/hhampson/Library/CloudStorage/OneDrive-SharedLibraries-UW/Laura Pyle - Biostatistics Core Shared Drive/Liver project/Results")

#Packages
library(tidyverse)

#Load significant results
hep1 <- read.csv(fs::path(dir.results,"Hep_1_Steatosis.csv")) %>% 
  filter(fdr<0.05) %>% 
  dplyr::select(-X)
hep2 <- read.csv(fs::path(dir.results,"Hep_2_Steatosis.csv"))%>% 
  filter(fdr<0.05)%>% 
  dplyr::select(-X)
hep3 <- read.csv(fs::path(dir.results,"Hep_3_Steatosis.csv"))%>% 
  filter(fdr<0.05)%>% 
  dplyr::select(-X)
hep4 <- read.csv(fs::path(dir.results,"Hep_4_Steatosis.csv"))%>% 
  filter(fdr<0.05)%>% 
  dplyr::select(-X)
hep5 <- read.csv(fs::path(dir.results,"Hep_5_Steatosis.csv"))%>% 
  filter(fdr<0.05)%>% 
  dplyr::select(-X)

#Merge together
dataframes <- list(hep1, hep2, hep3, hep4, hep5)

# Merge all dataframes by Gene (keep only genes that are present in all)
merged_df <- Reduce(function(x, y) merge(x, y, by = "Gene"), dataframes)

# Add column names to the merged dataframe to identify which dataset the beta belongs to
colnames(merged_df) <- c("Gene", paste("Beta", 1:5, sep = "_"))

# Create a new dataframe that stores the sign of betas across all datasets (Positive, Negative, NA)
sign_df <- merged_df
for (i in 2:ncol(merged_df)) {
  sign_df[[i]] <- ifelse(sign_df[[i]] > 0, "Positive", ifelse(sign_df[[i]] < 0, "Negative", "NA"))
}
