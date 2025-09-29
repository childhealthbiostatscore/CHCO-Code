### ROCKIES Module Scores 

library(tidyverse) 
library(stringr)
library(dplyr)
library(ggplot2)
library(lme4)
library(lmerTest)
library(broom.mixed)




module_scores <- data.table::fread('C:/Users/netio/Documents/UofW/Rockies/Rockies_updates_9.26.25/module_scores/HALLMARK_GO_pathways_modulescores.txt')






# List of module score columns (excluding metadata columns)
module_columns <- c("TCA_score1", "OxPhos_score1", "Response_to_Insulin1", 
                    "Insulin_Receptor_Signaling1", "Neg_Reg_Insulin_Signaling1", 
                    "Pos_Reg_Insulin_Signaling1", "IGF_Receptor_Binding1", 
                    "OxPhos_Hallmark1", "Fatty_Acid_Metabolism1", "Glycolysis1",
                    "mTORC1_Signaling1", "Adipogenesis1", "Peroxisome1", 
                    "Inflammatory_Response1", "IFN_Gamma_Response1", 
                    "IFN_Alpha_Response1", "IL6_JAK_STAT31", "TNFa_NF_kB1",
                    "Complement1", "Allograft_Rejection1", "Hypoxia1", 
                    "ROS_Pathway1", "PI3K_AKT_mTOR1")

# Function to run lmer for one module
run_lmer_module <- function(data, module_name, group_var = "group", 
                            person_var = "record_id") {
  
  # Create formula
  formula_str <- paste0(module_name, " ~ ", group_var, " + (1|", person_var, ")")
  
  # Fit model
  tryCatch({
    model <- lmer(as.formula(formula_str), data = data)
    
    # Extract results
    fixed_effects <- summary(model)$coefficients
    anova_result <- anova(model)
    
    # Return tidy results
    result <- data.frame(
      module = module_name,
      term = rownames(fixed_effects),
      estimate = fixed_effects[, "Estimate"],
      std_error = fixed_effects[, "Std. Error"],
      t_value = fixed_effects[, "t value"],
      p_value = fixed_effects[, "Pr(>|t|)"],
      stringsAsFactors = FALSE
    )
    
    return(result)
    
  }, error = function(e) {
    message(paste("Error fitting model for", module_name, ":", e$message))
    return(NULL)
  })
}

# Run models for all modules
all_results <- lapply(module_columns, function(module) {
  run_lmer_module(module_scores, module)
}) %>%
  bind_rows()

# Filter to group comparison rows only (not intercept)
group_results <- all_results %>%
  filter(term != "(Intercept)") %>%
  arrange(p_value)

# Add FDR correction
group_results$p_adj_fdr <- p.adjust(group_results$p_value, method = "fdr")

# View results
print(group_results)


write.table(group_results, '/Users/netio/Documents/UofW/Rockies/Rockies_updates_9.26.25/module_scores/LMER_AllCells_T2DvsLC_Results.txt', row.names=F, quote=F, sep='\t')





library(ggplot2)
library(ggrepel)

# Volcano plot



png('/Users/netio/Documents/UofW/Rockies/Rockies_updates_9.26.25/module_scores/LMER_AllCells_T2DvsLC_Results_Volcanoplot.png')
ggplot(group_results, aes(x = estimate, y = -log10(p_value))) +
  geom_point(aes(color = p_value < 0.05), alpha = 0.6, size = 3) +
  geom_text_repel(data = filter(group_results, p_value < 0.05),
                  aes(label = module), size = 3) +
  scale_color_manual(values = c("grey", "red")) +
  geom_hline(yintercept = -log10(0.05))+
  geom_vline(xintercept = 0)+
  theme_classic() +
  labs(x = "Effect Size (Estimate)", y = "-log10(p-value)")
dev.off()




############################# PT Cells Only 



module_scores <- data.table::fread('C:/Users/netio/Documents/UofW/Rockies/Rockies_updates_9.26.25/module_scores/HALLMARK_GO_pathways_modulescores.txt')
module_scores <- module_scores %>% filter(celltype2 == 'PT')





# List of module score columns (excluding metadata columns)
module_columns <- c("TCA_score1", "OxPhos_score1", "Response_to_Insulin1", 
                    "Insulin_Receptor_Signaling1", "Neg_Reg_Insulin_Signaling1", 
                    "Pos_Reg_Insulin_Signaling1", "IGF_Receptor_Binding1", 
                    "OxPhos_Hallmark1", "Fatty_Acid_Metabolism1", "Glycolysis1",
                    "mTORC1_Signaling1", "Adipogenesis1", "Peroxisome1", 
                    "Inflammatory_Response1", "IFN_Gamma_Response1", 
                    "IFN_Alpha_Response1", "IL6_JAK_STAT31", "TNFa_NF_kB1",
                    "Complement1", "Allograft_Rejection1", "Hypoxia1", 
                    "ROS_Pathway1", "PI3K_AKT_mTOR1")

# Function to run lmer for one module
run_lmer_module <- function(data, module_name, group_var = "group", 
                            person_var = "record_id") {
  
  # Create formula
  formula_str <- paste0(module_name, " ~ ", group_var, " + (1|", person_var, ")")
  
  # Fit model
  tryCatch({
    model <- lmer(as.formula(formula_str), data = data)
    
    # Extract results
    fixed_effects <- summary(model)$coefficients
    anova_result <- anova(model)
    
    # Return tidy results
    result <- data.frame(
      module = module_name,
      term = rownames(fixed_effects),
      estimate = fixed_effects[, "Estimate"],
      std_error = fixed_effects[, "Std. Error"],
      t_value = fixed_effects[, "t value"],
      p_value = fixed_effects[, "Pr(>|t|)"],
      stringsAsFactors = FALSE
    )
    
    return(result)
    
  }, error = function(e) {
    message(paste("Error fitting model for", module_name, ":", e$message))
    return(NULL)
  })
}

# Run models for all modules
all_results <- lapply(module_columns, function(module) {
  run_lmer_module(module_scores, module)
}) %>%
  bind_rows()

# Filter to group comparison rows only (not intercept)
group_results <- all_results %>%
  filter(term != "(Intercept)") %>%
  arrange(p_value)

# Add FDR correction
group_results$p_adj_fdr <- p.adjust(group_results$p_value, method = "fdr")

# View results
print(group_results)


write.table(group_results, '/Users/netio/Documents/UofW/Rockies/Rockies_updates_9.26.25/module_scores/LMER_AllPTCells_T2DvsLC_Results.txt', row.names=F, quote=F, sep='\t')





library(ggplot2)
library(ggrepel)

# Volcano plot



png('/Users/netio/Documents/UofW/Rockies/Rockies_updates_9.26.25/module_scores/LMER_AllPTCells_T2DvsLC_Results_Volcanoplot.png')
ggplot(group_results, aes(x = estimate, y = -log10(p_value))) +
  geom_point(aes(color = p_value < 0.05), alpha = 0.6, size = 3) +
  geom_text_repel(data = filter(group_results, p_value < 0.05),
                  aes(label = module), size = 3) +
  scale_color_manual(values = c("grey", "red")) +
  geom_hline(yintercept = -log10(0.05))+
  geom_vline(xintercept = 0)+
  theme_classic() +
  labs(x = "Effect Size (Estimate)", y = "-log10(p-value)")
dev.off()





#################### PT Cell Subtypes 

remove(list=ls())


celltypes <- c('PT-S1/S2', 'PT-S3', 'aPT')


for(celltype in celltypes){
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
}















