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





#################### Cell Subtypes 

remove(list=ls())


celltypes <- c('All', 'PT', 'TAL', 'EC', 'PT-S1/S2', 'PT-S3', 'aPT')


for(celltype in celltypes){
  
  
  celltype2 <- str_replace_all(celltype,"/","_")
  celltype2 <- str_replace_all(celltype2,"-","_")
  
  
  
  
  
  
  
  
  module_scores <- data.table::fread('C:/Users/netio/Documents/UofW/Rockies/Rockies_updates_9.26.25/module_scores/HALLMARK_GO_pathways_modulescores.txt')
  
  if(celltype == 'All'){
    module_scores <- module_scores
  }else if(celltype %in% c('PT', 'TAL', 'EC')){
    module_scores <- module_scores %>% 
      filter(celltype2 == celltype)
  }else{
    module_scores <- module_scores %>% filter(KPMP_celltype == celltype) 
  }
  
  
  
  
  
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
  
  
  write.table(group_results, 
              paste0('/Users/netio/Documents/UofW/Rockies/Rockies_updates_9.26.25/new_analyses/LMER_', celltype2, 'Cells_T2DvsLC_Results.txt'), 
              row.names=F, quote=F, sep='\t')
  
  
  
  
  
  library(ggplot2)
  library(ggrepel)
  
  # Volcano plot
  
  
  
  png(paste0('/Users/netio/Documents/UofW/Rockies/Rockies_updates_9.26.25/new_analyses/LMER_', celltype2, 'Cells_T2DvsLC_Results_Volcanoplot.png'))
  
  tmp_graph <- ggplot(group_results, aes(x = estimate, y = -log10(p_value))) +
    geom_point(aes(color = p_value < 0.05), alpha = 0.6, size = 3) +
    geom_text_repel(data = filter(group_results, p_value < 0.05),
                    aes(label = module), size = 3) +
    scale_color_manual(values = c("grey", "red")) +
    geom_hline(yintercept = -log10(0.05))+
    geom_vline(xintercept = 0)+
    theme_classic() +
    labs(x = "Effect Size (Estimate)", y = "-log10(p-value)")
  print(tmp_graph)
  dev.off()
  
  
  
  
  
  
  
}







########## Quantile Regressions

remove(list=ls())

module_scores <- data.table::fread('C:/Users/netio/Documents/UofW/Rockies/Rockies_updates_9.26.25/module_scores/HALLMARK_GO_pathways_modulescores.txt')

celltypes <- c('All', 'PT', 'TAL', 'EC', 'PT-S1/S2', 'PT-S3', 'aPT')



for(celltype in celltypes){
  
  module_scores <- data.table::fread('C:/Users/netio/Documents/UofW/Rockies/Rockies_updates_9.26.25/module_scores/HALLMARK_GO_pathways_modulescores.txt')
  
  
  celltype2 <- str_replace_all(celltype,"/","_")
  celltype2 <- str_replace_all(celltype2,"-","_")
  
  if(celltype == 'All'){
    module_scores <- module_scores
  }else if(celltype %in% c('PT', 'TAL', 'EC')){
    module_scores <- module_scores %>% 
      filter(celltype2 == 'PT')
  }else{
    module_scores <- module_scores %>%
      filter(KPMP_celltype == celltype)
  }
  

library(quantreg)

# Test at different quantiles (median, upper quartile)
person_means <- module_scores %>%
  group_by(record_id, group) %>%
  summarize(mean_TCA = mean(TCA_score1), .groups = 'drop')

# Median regression
rq_model <- rq(mean_TCA ~ group, tau = 0.5, data = person_means)
summary(rq_model)

# Upper quartile (metabolically active cells)
rq_model_75 <- rq(mean_TCA ~ group, tau = 0.75, data = person_means)
summary(rq_model_75)



library(quantreg)
library(dplyr)
library(ggplot2)
library(tidyr)

# Define all pathway columns (all score columns)
all_pathways <- c("TCA_score1", "OxPhos_score1", "Response_to_Insulin1", 
                  "Insulin_Receptor_Signaling1", "Neg_Reg_Insulin_Signaling1", 
                  "Pos_Reg_Insulin_Signaling1", "IGF_Receptor_Binding1", 
                  "OxPhos_Hallmark1", "Fatty_Acid_Metabolism1", "Glycolysis1",
                  "mTORC1_Signaling1", "Adipogenesis1", "Peroxisome1", 
                  "Inflammatory_Response1", "IFN_Gamma_Response1", 
                  "IFN_Alpha_Response1", "IL6_JAK_STAT31", "TNFa_NF_kB1",
                  "Complement1", "Allograft_Rejection1", "Hypoxia1", 
                  "ROS_Pathway1", "PI3K_AKT_mTOR1")

# Function to test multiple quantiles for one pathway
test_multiple_quantiles <- function(data, formula_str, quantiles = c(0.25, 0.5, 0.75)) {
  
  results <- lapply(quantiles, function(tau) {
    # Fit model
    model <- rq(as.formula(formula_str), tau = tau, data = data)
    
    # Get summary with bootstrap SE
    model_summary <- summary(model, se = "boot", R = 1000)
    
    # Extract coefficients (skip intercept)
    coef_df <- as.data.frame(model_summary$coefficients)
    coef_df$term <- rownames(coef_df)
    coef_df$quantile <- tau
    
    # Rename columns for clarity
    names(coef_df) <- c("estimate", "std_error", "t_value", "p_value", "term", "quantile")
    
    return(coef_df)
  })
  
  bind_rows(results) %>%
    filter(term != "(Intercept)")
}

# Function to test all pathways
test_all_pathways_qr <- function(data, pathways, quantiles = c(0.1, 0.25, 0.5, 0.75, 0.9)) {
  
  all_results <- list()
  
  for (pathway in pathways) {
    cat("Testing", pathway, "\n")
    
    # Aggregate to person level
    person_data <- data %>%
      group_by(record_id, group) %>%
      summarize(score = mean(.data[[pathway]], na.rm = TRUE), .groups = 'drop')
    
    # Skip if not enough data
    if (nrow(person_data) < 10) {
      cat("  Skipping - insufficient data\n")
      next
    }
    
    # Test quantiles
    tryCatch({
      pathway_results <- test_multiple_quantiles(
        person_data,
        "score ~ group",
        quantiles = quantiles
      )
      
      pathway_results$pathway <- pathway
      all_results[[pathway]] <- pathway_results
    }, error = function(e) {
      cat("  Error:", e$message, "\n")
    })
  }
  
  bind_rows(all_results)
}

# Run quantile regression for ALL pathways
qr_all_results <- test_all_pathways_qr(
  module_scores, 
  all_pathways,
  quantiles = c(0.1, 0.25, 0.5, 0.75, 0.9)
)

# View significant results (NO FDR CORRECTION)
qr_significant <- qr_all_results %>%
  filter(p_value < 0.05) %>%
  arrange(quantile, p_value)

cat("\n=== SIGNIFICANT RESULTS (p < 0.05) ===\n")
print(qr_significant)

# Summary by pathway
qr_summary <- qr_all_results %>%
  group_by(pathway) %>%
  summarize(
    sig_at_any_quantile = any(p_value < 0.05),
    sig_quantiles = paste(quantile[p_value < 0.05], collapse = ", "),
    n_sig_quantiles = sum(p_value < 0.05),
    strongest_effect_quantile = quantile[which.min(p_value)],
    strongest_effect_pval = min(p_value),
    effect_direction = sign(estimate[which.min(p_value)]),
    mean_effect = mean(estimate),
    .groups = 'drop'
  ) %>%
  arrange(desc(sig_at_any_quantile), strongest_effect_pval)

cat("\n=== PATHWAY SUMMARY ===\n")
print(qr_summary)

# Export results
write.csv(qr_all_results, 
          paste0("/Users/netio/Documents/UofW/Rockies/Rockies_updates_9.26.25/new_analyses/quantile_regression_", celltype2, "_cells_all_results.csv"), 
          row.names = FALSE)
write.csv(qr_summary, 
          paste0("/Users/netio/Documents/UofW/Rockies/Rockies_updates_9.26.25/new_analyses/quantile_regression_", celltype2, "_cells_summary.csv"), 
          row.names = FALSE)

# ===== VISUALIZATIONS =====

# 1. Heatmap of -log10(p-values) with corrplot colors
p1 <- ggplot(qr_all_results, aes(x = as.factor(quantile), y = pathway, 
                                 fill = -log10(p_value))) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "white", mid = "yellow", high = "red",
                       midpoint = -log10(0.05),
                       name = "-log10(p)") +
  geom_text(aes(label = ifelse(p_value < 0.05, "*", "")), size = 4) +
  labs(title = "Quantile Regression P-values Across All Pathways",
       subtitle = "* indicates p < 0.05",
       x = "Quantile (τ)",
       y = "Pathway") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 8),
        axis.text.x = element_text(size = 10))

print(p1)
ggsave(paste0("qr_", celltype2, "_cells_heatmap_pvalues.pdf"), 
       p1, width = 10, height = 12)

# 2. Effect sizes across quantiles for significant pathways
sig_pathways <- qr_summary %>% 
  filter(sig_at_any_quantile) %>% 
  pull(pathway)

if (length(sig_pathways) > 0) {
  qr_sig_data <- qr_all_results %>% 
    filter(pathway %in% sig_pathways)
  
  p2 <- ggplot(qr_sig_data, aes(x = quantile, y = estimate, color = pathway)) +
    geom_line(size = 1) +
    geom_ribbon(aes(ymin = estimate - 1.96*std_error, 
                    ymax = estimate + 1.96*std_error,
                    fill = pathway), 
                alpha = 0.2, color = NA) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    facet_wrap(~pathway, scales = "free_y", ncol = 3) +
    labs(title = "Group Effect Across Distribution Quantiles",
         subtitle = "Only pathways with p < 0.05 at any quantile",
         x = "Quantile (τ)",
         y = "Effect Size (Estimate)") +
    theme_minimal() +
    theme(legend.position = "none")
  
  print(p2)
  ggsave(paste0("/Users/netio/Documents/UofW/Rockies/Rockies_updates_9.26.25/new_analyses/qr_", celltype2, "_cells_effect_sizes_significant.pdf"), 
         p2, width = 14, height = 10)
}

# 3. Heatmap of effect sizes - RED = DOWN (negative), BLUE = UP (positive)
p3 <- ggplot(qr_all_results, aes(x = as.factor(quantile), y = pathway, 
                                 fill = estimate)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "red", mid = "white", high = "blue",
                       midpoint = 0,
                       name = "Effect Size") +
  geom_text(aes(label = ifelse(p_value < 0.05, "*", "")), size = 4) +
  labs(title = "Effect Sizes Across All Pathways",
       subtitle = "Red = Down, Blue = Up, * indicates p < 0.05",
       x = "Quantile (τ)",
       y = "Pathway") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 8),
        axis.text.x = element_text(size = 10))

print(p3)
ggsave("qr_heatmap_effects.pdf", p3, width = 10, height = 12)

# 4. Summary bar plot
if (sum(qr_summary$sig_at_any_quantile) > 0) {
  p4 <- ggplot(qr_summary %>% filter(sig_at_any_quantile), 
               aes(x = reorder(pathway, n_sig_quantiles), y = n_sig_quantiles)) +
    geom_col(aes(fill = mean_effect)) +
    coord_flip() +
    scale_fill_gradient2(low = "red", high = "blue", mid = "white", midpoint = 0,
                         name = "Mean Effect\n(Red=Down)") +
    labs(title = "Number of Significant Quantiles per Pathway",
         subtitle = "No FDR correction applied",
         x = "Pathway",
         y = "Number of Quantiles with p < 0.05") +
    theme_minimal()
  
  print(p4)
  ggsave(paste0("/Users/netio/Documents/UofW/Rockies/Rockies_updates_9.26.25/new_analyses/qr_", celltype2, "_cells_summary_barplot.pdf"), p4, width = 10, height = 8)
}

# 5. Additional: Combined heatmap showing both effect and significance
qr_all_results$sig_label <- ifelse(qr_all_results$p_value < 0.05, 
                                   ifelse(qr_all_results$p_value < 0.01, "**", "*"), 
                                   "")

p5 <- ggplot(qr_all_results, aes(x = as.factor(quantile), y = pathway, 
                                 fill = estimate)) +
  geom_tile(color = "white", size = 0.5) +
  scale_fill_gradient2(low = "red", mid = "white", high = "blue",
                       midpoint = 0,
                       limits = c(-max(abs(qr_all_results$estimate)), 
                                  max(abs(qr_all_results$estimate))),
                       name = "Effect Size\n(Red=Down\nBlue=Up)") +
  geom_text(aes(label = sig_label), size = 5, fontface = "bold") +
  labs(title = "Quantile Regression: Effect Sizes and Significance",
       subtitle = "* p < 0.05, ** p < 0.01",
       x = "Quantile (τ)",
       y = "Pathway") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 9),
        axis.text.x = element_text(size = 10, face = "bold"),
        panel.grid = element_blank())

print(p5)
ggsave(paste0("/Users/netio/Documents/UofW/Rockies/Rockies_updates_9.26.25/new_analyses/qr_", celltype2, "_cells_combined_heatmap.pdf"), 
       p5, width = 11, height = 13)

cat("\n=== ANALYSIS COMPLETE ===\n")
cat("Files saved:\n")
cat("  - quantile_regression_all_results.csv\n")
cat("  - quantile_regression_summary.csv\n")
cat("  - qr_heatmap_pvalues.pdf\n")
cat("  - qr_heatmap_effects.pdf (RED=down, BLUE=up)\n")
cat("  - qr_effect_sizes_significant.pdf\n")
cat("  - qr_summary_barplot.pdf\n")
cat("  - qr_combined_heatmap.pdf (NEW!)\n")



}








#### K-S Testing

remove(list=ls())


celltypes <- c('All', 'PT', 'TAL', 'EC', 'PT-S1/S2', 'PT-S3', 'aPT')



for(celltype in celltypes){
  
  module_scores <- data.table::fread('C:/Users/netio/Documents/UofW/Rockies/Rockies_updates_9.26.25/module_scores/HALLMARK_GO_pathways_modulescores.txt')
  
  # Filter by cell type if not "All"
  if(celltype == 'All'){
    module_scores <- module_scores
  }else if(celltype %in% c('PT', 'TAL', 'EC')){
    module_scores <- module_scores %>% filter(celltype2 == celltype)
  }else{
    module_scores <- module_scores %>% filter(KPMP_celltype == celltype)
  }
  
  
  
  celltype2 <- str_replace_all(celltype,"/","_")
  celltype2 <- str_replace_all(celltype2,"-","_")
  
  
  
  
  
  all_pathways <- c("TCA_score1", "OxPhos_score1", "Response_to_Insulin1", 
                    "Insulin_Receptor_Signaling1", "Neg_Reg_Insulin_Signaling1", 
                    "Pos_Reg_Insulin_Signaling1", "IGF_Receptor_Binding1", 
                    "OxPhos_Hallmark1", "Fatty_Acid_Metabolism1", "Glycolysis1",
                    "mTORC1_Signaling1", "Adipogenesis1", "Peroxisome1", 
                    "Inflammatory_Response1", "IFN_Gamma_Response1", 
                    "IFN_Alpha_Response1", "IL6_JAK_STAT31", "TNFa_NF_kB1",
                    "Complement1", "Allograft_Rejection1", "Hypoxia1", 
                    "ROS_Pathway1", "PI3K_AKT_mTOR1")
  
  # Test for distribution differences (non-parametric)
  ks_test_all_pathways <- function(data, pathways) {
    
    results <- lapply(pathways, function(pathway) {
      # Aggregate to person level
      person_data <- data %>%
        group_by(record_id, group) %>%
        summarize(score = mean(.data[[pathway]], na.rm = TRUE), .groups = 'drop')
      
      # Split by group
      groups <- unique(person_data$group)
      group1_scores <- person_data %>% filter(group == groups[1]) %>% pull(score)
      group2_scores <- person_data %>% filter(group == groups[2]) %>% pull(score)
      
      # KS test
      ks_result <- ks.test(group1_scores, group2_scores)
      
      data.frame(
        pathway = pathway,
        D_statistic = ks_result$statistic,
        p_value = ks_result$p.value
      )
    })
    
    bind_rows(results) %>%
      mutate(p_adj_fdr = p.adjust(p_value, method = "fdr"))
  }
  
  ks_results <- ks_test_all_pathways(module_scores, all_pathways)
  
  # Significant distribution differences
  ks_significant <- ks_results %>%
    arrange(p_adj_fdr) %>%
    mutate(
      pathway_index = row_number(),
      # Create label only if p < 0.1
      pathway_label = ifelse(p_value < 0.1, 
                             gsub("_score1|1$", "", pathway),  # Clean up pathway names
                             "")
    )
  
  print(ks_significant)
  
  # Create plot
  p <- ggplot(ks_significant, aes(x = pathway_index, y = -log10(p_value))) +
    geom_point(aes(color = p_value < 0.1), size = 2) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
    geom_hline(yintercept = -log10(0.1), linetype = "dashed", color = "orange") +
    geom_text(aes(label = pathway_label), 
              hjust = -0.1, vjust = 0.5, size = 3, angle = 45) +
    scale_color_manual(values = c("FALSE" = "gray60", "TRUE" = "darkblue"),
                       name = "p < 0.1") +
    labs(title = paste0("KS Test Results - ", celltype),
         x = "Pathways",
         y = "-log10(p-value)") +
    theme_minimal() +
    theme(
      axis.text.x = element_blank(),  # Remove x-axis labels
      axis.ticks.x = element_blank(), # Remove x-axis ticks
      panel.grid.major.x = element_blank(),
      plot.title = element_text(hjust = 0.5, face = "bold")
    )
  
  print(p)
  
  # Save plot as PNG
  ggsave(paste0("/Users/netio/Documents/UofW/Rockies/Rockies_updates_9.26.25/new_analyses/KS_test_", celltype2, ".png"), 
         plot = p, width = 10, height = 6, dpi = 300)
}














########### Test interaction between group and gene expression

library(lme4)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(patchwork)
library(stringr)

# Load and prepare data
harmonized_data <- read.csv("C:/Users/netio/OneDrive - UW/Laura Pyle's files - Biostatistics Core Shared Drive/Data Harmonization/Data Clean/harmonized_dataset.csv", na = '')

dat <- harmonized_data %>% 
  dplyr::select(-dob) %>% 
  arrange(date_of_screen) %>% 
  dplyr::summarise(
    across(where(negate(is.numeric)), ~ ifelse(all(is.na(.x)), NA_character_, last(na.omit(.x)))),
    across(where(is.numeric), ~ ifelse(all(is.na(.x)), NA_real_, mean(na.omit(.x), na.rm=T))),
    .by = c(record_id, visit)
  )

dat2 <- dat %>% 
  mutate(avg_c_k2_f = avg_c_k2 / avg_c_f) %>%
  dplyr::select(mrn, avg_c_k2, avg_c_f, avg_c_k2_f)

# Define pathways, cell types, and PET outcomes
all_pathways <- c("TCA_score1", "OxPhos_score1", "Response_to_Insulin1", 
                  "Insulin_Receptor_Signaling1", "Neg_Reg_Insulin_Signaling1", 
                  "Pos_Reg_Insulin_Signaling1", "IGF_Receptor_Binding1", 
                  "OxPhos_Hallmark1", "Fatty_Acid_Metabolism1", "Glycolysis1",
                  "mTORC1_Signaling1", "Adipogenesis1", "Peroxisome1", 
                  "Inflammatory_Response1", "IFN_Gamma_Response1", 
                  "IFN_Alpha_Response1", "IL6_JAK_STAT31", "TNFa_NF_kB1",
                  "Complement1", "Allograft_Rejection1", "Hypoxia1", 
                  "ROS_Pathway1", "PI3K_AKT_mTOR1")

celltypes <- c('All', 'PT', 'TAL', 'EC', 'PT-S1/S2', 'PT-S3', 'aPT')

pet_outcomes <- c("avg_c_k2", "avg_c_f", "avg_c_k2_f")

# Function to test interaction for one pathway and one outcome
test_interaction <- function(data, pathway, outcome) {
  
  # Filter complete cases
  interaction_data <- data %>%
    filter(!is.na(.data[[outcome]]), !is.na(.data[[pathway]])) %>%
    select(record_id, group, all_of(pathway), all_of(outcome))
  
  # Skip if insufficient data
  if(nrow(interaction_data) < 10) {
    return(NULL)
  }
  
  # Fit interaction model
  formula_str <- paste0(outcome, " ~ ", pathway, " * group")
  model <- lm(as.formula(formula_str), data = interaction_data)
  model_summary <- summary(model)
  
  # Extract interaction p-value
  coef_table <- coef(model_summary)
  interaction_term <- grep(":", rownames(coef_table), value = TRUE)
  
  if(length(interaction_term) > 0) {
    interaction_p <- coef_table[interaction_term, "Pr(>|t|)"]
  } else {
    interaction_p <- NA
  }
  
  return(list(
    data = interaction_data,
    model = model,
    interaction_p = interaction_p,
    pathway = pathway,
    outcome = outcome
  ))
}

# Function to create IMPROVED plot for one pathway-outcome combination
create_pathway_plot <- function(result) {
  
  if(is.null(result)) return(NULL)
  
  # Clean pathway name - remove suffixes and make readable
  pathway_clean <- gsub("_score1|1$", "", result$pathway)
  pathway_clean <- gsub("_", " ", pathway_clean)
  
  # Determine significance marker
  sig_marker <- case_when(
    is.na(result$interaction_p) ~ "",
    result$interaction_p < 0.001 ~ "***",
    result$interaction_p < 0.01 ~ "**",
    result$interaction_p < 0.05 ~ "*",
    TRUE ~ ""  # Don't show "ns" to reduce clutter
  )
  
  # Create cleaner title
  if(!is.na(result$interaction_p) && result$interaction_p < 0.1) {
    title_text <- paste0(pathway_clean, " ", sig_marker, 
                         sprintf("\np=%.3f", result$interaction_p))
  } else {
    title_text <- pathway_clean
  }
  
  # Clean outcome name for y-axis
  outcome_clean <- case_when(
    result$outcome == "avg_c_k2" ~ "Cortical K2",
    result$outcome == "avg_c_f" ~ "Cortical F",
    result$outcome == "avg_c_k2_f" ~ "Cortical K2/F",
    TRUE ~ result$outcome
  )
  
  # Create IMPROVED plot
  p <- ggplot(result$data, aes(x = .data[[result$pathway]], y = .data[[result$outcome]], color = group)) +
    # Use semi-transparent points
    geom_point(alpha = 0.5, size = 1.5, shape = 16) +
    # Add regression lines WITHOUT confidence bands (less clutter)
    geom_smooth(method = "lm", se = FALSE, linewidth = 1.2) +
    # Better color scheme
    scale_color_manual(
      values = c("Lean_Control" = "#3498db", "Type_2_Diabetes" = "#e74c3c"),
      labels = c("Lean_Control" = "Lean", "Type_2_Diabetes" = "T2D")
    ) +
    labs(
      title = title_text,
      x = NULL,  # Remove x-axis label to reduce clutter
      y = NULL   # Will add common y-axis label to overall plot
    ) +
    theme_minimal(base_size = 10) +
    theme(
      plot.title = element_text(size = 9, face = "bold", hjust = 0.5),
      axis.text = element_text(size = 7),
      legend.position = "none",
      panel.grid.minor = element_blank(),  # Remove minor grid lines
      panel.grid.major = element_line(size = 0.3, color = "grey90"),  # Lighter grid
      plot.background = element_rect(fill = "white", color = NA),
      panel.border = element_rect(color = "grey80", fill = NA, size = 0.5)
    )
  
  return(p)
}

# Main loop through cell types and PET outcomes
for(celltype in celltypes){
  
  # Load module scores
  module_scores <- data.table::fread('C:/Users/netio/Documents/UofW/Rockies/Rockies_updates_9.26.25/module_scores/HALLMARK_GO_pathways_modulescores.txt')
  
  
  # Filter by cell type
  if(celltype == 'All'){
    module_scores <- module_scores
  } else if(celltype %in% c('PT', 'TAL', 'EC')){
    module_scores <- module_scores %>% filter(celltype2 == celltype)
  } else {
    module_scores <- module_scores %>% filter(KPMP_celltype == celltype)
  }
  
  
  
  
  module_scores <- module_scores %>%
    group_by(record_id, mrn, group) %>%
    summarize(
      # Calculate mean for all numeric module score columns
      across(c(TCA_score1, OxPhos_score1, 
               # GO terms
               Response_to_Insulin1, Insulin_Receptor_Signaling1, 
               Neg_Reg_Insulin_Signaling1, Pos_Reg_Insulin_Signaling1, 
               IGF_Receptor_Binding1,
               # HALLMARK pathways - Metabolism
               OxPhos_Hallmark1, Fatty_Acid_Metabolism1, Glycolysis1, 
               mTORC1_Signaling1, Adipogenesis1, Peroxisome1,
               # HALLMARK pathways - Immune/Inflammatory
               Inflammatory_Response1, IFN_Gamma_Response1, IFN_Alpha_Response1,
               IL6_JAK_STAT31, TNFa_NF_kB1, Complement1, Allograft_Rejection1,
               # HALLMARK pathways - Cross-cutting
               Hypoxia1, ROS_Pathway1, PI3K_AKT_mTOR1), 
             list(mean = ~mean(.x, na.rm = TRUE)),
             .names = "{.col}"),
      .groups = 'drop'
    )
  
  
  
  
  
  # Combine with PET data
  combined_df <- module_scores %>% left_join(dat2, by='mrn')
  
  # Clean cell type name for filename
  celltype2 <- str_replace_all(celltype, "/", "_")
  celltype2 <- str_replace_all(celltype2, "-", "_")
  
  cat("Processing cell type:", celltype, "\n")
  
  # Loop through each PET outcome
  for(outcome in pet_outcomes) {
    
    cat("  Testing outcome:", outcome, "\n")
    
    # Test interactions for all pathways with this outcome
    results_list <- lapply(all_pathways, function(pathway) {
      test_interaction(combined_df, pathway, outcome)
    })
    
    # Create plots for all pathways
    plot_list <- lapply(results_list, create_pathway_plot)
    
    # Remove NULL plots
    plot_list <- plot_list[!sapply(plot_list, is.null)]
    
    if(length(plot_list) == 0) {
      cat("    No plots generated for", outcome, "\n")
      next
    }
    
    # Calculate layout - use fixed 5 columns for consistency
    n_plots <- length(plot_list)
    n_cols <- 5
    n_rows <- ceiling(n_plots / n_cols)
    
    # Clean outcome name for title
    outcome_clean <- case_when(
      outcome == "avg_c_k2" ~ "Cortical K2",
      outcome == "avg_c_f" ~ "Cortical F",
      outcome == "avg_c_k2_f" ~ "Cortical K2/F",
      TRUE ~ outcome
    )
    
    # Create combined plot with patchwork
    combined_plot <- wrap_plots(plot_list, ncol = n_cols) +
      plot_annotation(
        title = paste0("Group × Gene Expression Interactions: ", celltype),
        subtitle = paste0("Outcome: ", outcome_clean, " | *p<0.05  **p<0.01  ***p<0.001"),
        theme = theme(
          plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
          plot.subtitle = element_text(size = 12, hjust = 0.5, color = "grey30")
        )
      )
    
    # Add a single legend at the bottom
    combined_plot <- combined_plot & 
      theme(legend.position = "bottom") &
      guides(color = guide_legend(
        title = "Group",
        override.aes = list(size = 3, alpha = 1)
      ))
    
    # Clean outcome name for filename
    outcome_file <- str_replace_all(outcome, "_", "")
    
    # Save plot with better dimensions
    ggsave(
      paste0("C:/Users/netio/Documents/UofW/Rockies/Rockies_updates_9.26.25/new_analyses/SlopeDifferences_", celltype2, "_", outcome_file, ".png"),
      plot = combined_plot, 
      width = 16,  # Fixed width for consistency
      height = 3 * n_rows, 
      dpi = 300,
      limitsize = FALSE,
      bg = "white"
    )
    
    cat("    Saved plot with", n_plots, "pathways\n")
    
    # Create and save summary table for this outcome
    summary_table <- data.frame(
      pathway = sapply(results_list, function(x) if(!is.null(x)) x$pathway else NA),
      outcome = outcome,
      interaction_p = sapply(results_list, function(x) if(!is.null(x)) x$interaction_p else NA)
    ) %>%
      filter(!is.na(pathway)) %>%
      arrange(interaction_p) %>%
      mutate(
        significant = interaction_p < 0.05,
        p_adj_fdr = p.adjust(interaction_p, method = "fdr")
      )
    
    write.csv(summary_table, 
              paste0("C:/Users/netio/Documents/UofW/Rockies/Rockies_updates_9.26.25/new_analyses/SlopeDifferences_", celltype2, "_", outcome_file, "_summary.csv"),
              row.names = FALSE)
  }
}












############ Metabolic per person graphing 

library(dplyr)
library(ggplot2)
library(pheatmap)
library(viridis)
library(tidyr)
library(stringr)

# Define cell types and metabolic pathways
celltypes <- c('All', 'PT', 'TAL', 'EC', 'PT-S1/S2', 'PT-S3', 'aPT')
metabolic_pathways <- c("TCA_score1", "OxPhos_score1", "OxPhos_Hallmark1", 
                        "Fatty_Acid_Metabolism1", "Glycolysis1", "mTORC1_Signaling1")

# Loop through cell types
for(celltype in celltypes){
  
  cat("\nProcessing cell type:", celltype, "\n")
  
  # Load module scores
  module_scores <- data.table::fread('C:/Users/netio/Documents/UofW/Rockies/Rockies_updates_9.26.25/module_scores/HALLMARK_GO_pathways_modulescores.txt')
  
  # Filter by cell type
  if(celltype == 'All'){
    module_scores <- module_scores
  } else if(celltype %in% c('PT', 'TAL', 'EC')){
    module_scores <- module_scores %>% filter(celltype2 == celltype)
  } else {
    module_scores <- module_scores %>% filter(KPMP_celltype == celltype)
  }
  
  # Check if enough data
  if(nrow(module_scores) == 0) {
    cat("  No data for", celltype, "\n")
    next
  }
  
  # Aggregate to participant level (average across all cells per person)
  participant_metabolic <- module_scores %>%
    group_by(record_id, group) %>%
    summarize(across(all_of(metabolic_pathways), mean, na.rm = TRUE), .groups = 'drop') %>%
    arrange(group, record_id)
  
  # Check if enough participants
  if(nrow(participant_metabolic) < 5) {
    cat("  Insufficient participants for", celltype, "\n")
    next
  }
  
  # Clean cell type name for filename
  celltype2 <- str_replace_all(celltype, "/", "_")
  celltype2 <- str_replace_all(celltype2, "-", "_")
  
  # Create matrix for heatmap (participants x pathways)
  heatmap_matrix <- participant_metabolic %>%
    select(all_of(metabolic_pathways)) %>%
    as.matrix()
  
  # Add row names (participant IDs with group)
  rownames(heatmap_matrix) <- paste0(participant_metabolic$record_id, " (", 
                                     gsub("_", " ", participant_metabolic$group), ")")
  
  # Clean column names
  colnames(heatmap_matrix) <- gsub("_score1|1$", "", colnames(heatmap_matrix))
  colnames(heatmap_matrix) <- gsub("_", " ", colnames(heatmap_matrix))
  
  # Create annotation for groups
  annotation_row <- data.frame(
    Group = participant_metabolic$group,
    row.names = rownames(heatmap_matrix)
  )
  
  # Define colors for groups
  ann_colors <- list(
    Group = c("Lean_Control" = "#3498db", "Type_2_Diabetes" = "#e74c3c")
  )
  
  # Create pheatmap (RED = LOW, BLUE = HIGH)
  tryCatch({
    pheatmap(
      heatmap_matrix,
      color = colorRampPalette(c("#e74c3c", "white", "#3498db"))(100),
      scale = "column",
      cluster_rows = TRUE,
      cluster_cols = TRUE,
      annotation_row = annotation_row,
      annotation_colors = ann_colors,
      show_rownames = TRUE,
      show_colnames = TRUE,
      fontsize_row = 8,
      fontsize_col = 10,
      main = paste0("Metabolic Pathway Scores by Participant: ", celltype),
      filename = paste0("C:/Users/netio/Documents/UofW/Rockies/Rockies_updates_9.26.25/new_analyses/Participant_Metabolic_Heatmap_", celltype2, "_pheatmap.png"),
      width = 8,
      height = 12,
      dpi = 300
    )
  }, error = function(e) {
    cat("  Error creating pheatmap for", celltype, ":", e$message, "\n")
  })
  
  # Prepare data for ggplot
  heatmap_long <- participant_metabolic %>%
    mutate(participant_label = paste0(record_id, "\n(", str_replace(group, "_", " "), ")")) %>%
    pivot_longer(cols = all_of(metabolic_pathways), 
                 names_to = "pathway", 
                 values_to = "score") %>%
    mutate(pathway_clean = gsub("_score1|1$", "", pathway),
           pathway_clean = gsub("_", " ", pathway_clean))
  
  # Z-score normalization within each pathway
  heatmap_long <- heatmap_long %>%
    group_by(pathway) %>%
    mutate(score_scaled = scale(score)[,1]) %>%
    ungroup()
  
  # Order participants by group and average score
  participant_order <- heatmap_long %>%
    group_by(record_id, group, participant_label) %>%
    summarize(avg_score = mean(score_scaled, na.rm = TRUE), .groups = 'drop') %>%
    arrange(group, avg_score) %>%
    pull(participant_label)
  
  heatmap_long$participant_label <- factor(heatmap_long$participant_label, 
                                           levels = participant_order)
  
  # Create ggplot heatmap (RED = LOW, BLUE = HIGH)
  p <- ggplot(heatmap_long, aes(x = pathway_clean, y = participant_label, fill = score_scaled)) +
    geom_tile(color = "white", linewidth = 0.5) +
    scale_fill_gradient2(
      low = "#e74c3c",
      mid = "white",
      high = "#3498db",
      midpoint = 0,
      name = "Z-score\n(Red=Low\nBlue=High)",
      limits = c(-3, 3),
      oob = scales::squish
    ) +
    labs(
      title = paste0("Metabolic Pathway Expression by Participant: ", celltype),
      subtitle = "Averaged across all cells per participant | Z-score normalized by pathway",
      x = "Metabolic Pathway",
      y = "Participant ID (Group)"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
      axis.text.y = element_text(size = 7),
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 10, hjust = 0.5, color = "grey30"),
      legend.position = "right",
      panel.grid = element_blank()
    )
  
  ggsave(
    paste0("C:/Users/netio/Documents/UofW/Rockies/Rockies_updates_9.26.25/new_analyses/Participant_Metabolic_Heatmap_", celltype2, "_ggplot.png"),
    plot = p,
    width = 10,
    height = 14,
    dpi = 300
  )
  
  cat("  Saved heatmaps for", celltype, "\n")
  
  # Summary statistics
  cat("\nSummary for", celltype, ":\n")
  summary_stats <- participant_metabolic %>%
    group_by(group) %>%
    summarize(
      n_participants = n(),
      across(all_of(metabolic_pathways), 
             list(mean = ~mean(., na.rm=TRUE), sd = ~sd(., na.rm=TRUE)),
             .names = "{.col}_{.fn}")
    )
  print(summary_stats)
  
  # Save summary table
  write.csv(
    summary_stats,
    paste0("C:/Users/netio/Documents/UofW/Rockies/Rockies_updates_9.26.25/new_analyses/Participant_Metabolic_Summary_", celltype2, ".csv"),
    row.names = FALSE
  )
}















########### Quadratic Equations




library(lme4)
library(lmerTest)
library(dplyr)
library(ggplot2)
library(broom.mixed)
library(stringr)

# Load data
module_scores <- data.table::fread('C:/Users/netio/Documents/UofW/Rockies/Rockies_updates_9.26.25/module_scores/HALLMARK_GO_pathways_modulescores.txt')

harmonized_data <- read.csv("C:/Users/netio/OneDrive - UW/Laura Pyle's files - Biostatistics Core Shared Drive/Data Harmonization/Data Clean/harmonized_dataset.csv", na = '')

dat <- harmonized_data %>% 
  dplyr::select(-dob) %>% 
  arrange(date_of_screen) %>% 
  dplyr::summarise(
    across(where(negate(is.numeric)), ~ ifelse(all(is.na(.x)), NA_character_, last(na.omit(.x)))),
    across(where(is.numeric), ~ ifelse(all(is.na(.x)), NA_real_, mean(na.omit(.x), na.rm=T))),
    .by = c(record_id, visit)
  )

dat2 <- dat %>% 
  mutate(avg_c_k2_f = avg_c_k2 / avg_c_f) %>%
  dplyr::select(mrn, avg_c_k2, avg_c_f, avg_c_k2_f)

# Combine data
combined_df <- module_scores %>% left_join(dat2, by='mrn')

# Define pathways and outcomes to test
metabolic_pathways <- c("TCA_score1", "OxPhos_score1", "OxPhos_Hallmark1", 
                        "Fatty_Acid_Metabolism1", "Glycolysis1", "mTORC1_Signaling1")
pet_outcomes <- c("avg_c_k2", "avg_c_f", "avg_c_k2_f")

# Function to test quadratic relationship
test_quadratic_lmer <- function(data, pathway, outcome) {
  
  # Filter complete cases
  test_data <- data %>%
    filter(!is.na(.data[[outcome]]), !is.na(.data[[pathway]])) %>%
    select(record_id, group, all_of(pathway), all_of(outcome))
  
  # Need at least 20 observations for quadratic models
  if(nrow(test_data) < 20) {
    return(NULL)
  }
  
  # Center the predictor (important for quadratic terms)
  test_data <- test_data %>%
    mutate(pathway_centered = scale(.data[[pathway]], center = TRUE, scale = FALSE)[,1],
           pathway_squared = pathway_centered^2)
  
  # Fit models
  tryCatch({
    # Model 1: Linear only
    formula_linear <- paste0(outcome, " ~ pathway_centered + group + (1|record_id)")
    model_linear <- lmer(as.formula(formula_linear), data = test_data)
    
    # Model 2: Quadratic (parabolic)
    formula_quad <- paste0(outcome, " ~ pathway_centered + pathway_squared + group + (1|record_id)")
    model_quad <- lmer(as.formula(formula_quad), data = test_data)
    
    # Model 3: Quadratic with group interaction
    formula_quad_int <- paste0(outcome, " ~ (pathway_centered + pathway_squared) * group + (1|record_id)")
    model_quad_int <- lmer(as.formula(formula_quad_int), data = test_data)
    
    # Extract quadratic term significance
    quad_summary <- summary(model_quad)
    quad_coef <- coef(quad_summary)
    quad_p <- quad_coef["pathway_squared", "Pr(>|t|)"]
    quad_estimate <- quad_coef["pathway_squared", "Estimate"]
    
    # Extract interaction term significance (if exists)
    int_summary <- summary(model_quad_int)
    int_coef <- coef(int_summary)
    
    # Check for quadratic*group interaction
    int_term <- grep("pathway_squared.*group|group.*pathway_squared", rownames(int_coef), value = TRUE)
    if(length(int_term) > 0) {
      int_p <- int_coef[int_term[1], "Pr(>|t|)"]
      int_estimate <- int_coef[int_term[1], "Estimate"]
    } else {
      int_p <- NA
      int_estimate <- NA
    }
    
    # Calculate AIC/BIC
    model_comparison <- data.frame(
      pathway = pathway,
      outcome = outcome,
      linear_AIC = AIC(model_linear),
      quad_AIC = AIC(model_quad),
      quad_int_AIC = AIC(model_quad_int),
      linear_BIC = BIC(model_linear),
      quad_BIC = BIC(model_quad),
      quad_int_BIC = BIC(model_quad_int),
      quad_term_estimate = quad_estimate,
      quad_term_p = quad_p,
      quad_int_estimate = int_estimate,
      quad_int_p = int_p,
      best_model = ifelse(AIC(model_quad_int) == min(c(AIC(model_linear), AIC(model_quad), AIC(model_quad_int))), 
                          "Quadratic_Interaction",
                          ifelse(AIC(model_quad) < AIC(model_linear), "Quadratic", "Linear"))
    )
    
    return(list(
      comparison = model_comparison,
      model_linear = model_linear,
      model_quad = model_quad,
      model_quad_int = model_quad_int,
      data = test_data
    ))
    
  }, error = function(e) {
    cat("    Error fitting models:", e$message, "\n")
    return(NULL)
  })
}

# Test all pathway-outcome combinations
celltypes <- c('All', 'PT', 'TAL', 'EC', 'PT-S1/S2', 'PT-S3', 'aPT')

for(celltype in celltypes) {
  
  cat("\n=== Processing cell type:", celltype, "===\n")
  
  # Filter by cell type
  if(celltype == 'All'){
    filtered_data <- combined_df
  } else if(celltype %in% c('PT', 'TAL', 'EC')){
    filtered_data <- combined_df %>% filter(celltype2 == celltype)
  } else {
    filtered_data <- combined_df %>% filter(KPMP_celltype == celltype)
  }
  
  # Check if enough data
  if(nrow(filtered_data) < 20) {
    cat("  Insufficient data for", celltype, "\n")
    next
  }
  
  celltype_results <- list()
  
  for(outcome in pet_outcomes) {
    for(pathway in metabolic_pathways) {
      
      cat("  Testing:", pathway, "vs", outcome, "\n")
      
      result <- test_quadratic_lmer(filtered_data, pathway, outcome)
      
      if(!is.null(result)) {
        result$comparison$celltype <- celltype
        celltype_results[[paste(outcome, pathway, sep = "_")]] <- result
      }
    }
  }
  
  # Compile summary table for this cell type
  if(length(celltype_results) > 0) {
    summary_table <- lapply(celltype_results, function(x) x$comparison) %>% bind_rows()
    
    # Identify significant quadratic relationships
    summary_table <- summary_table %>%
      mutate(
        quadratic_significant = quad_term_p < 0.05,
        interaction_significant = !is.na(quad_int_p) & quad_int_p < 0.05,
        delta_AIC_quad_vs_linear = linear_AIC - quad_AIC,  # Positive = quadratic is better
        delta_AIC_int_vs_quad = quad_AIC - quad_int_AIC     # Positive = interaction is better
      ) %>%
      arrange(quad_term_p)
    
    # Clean celltype name for filename
    celltype2 <- str_replace_all(celltype, "/", "_")
    celltype2 <- str_replace_all(celltype2, "-", "_")
    
    # Save results for this cell type
    write.csv(summary_table, 
              paste0("C:/Users/netio/Documents/UofW/Rockies/Rockies_updates_9.26.25/new_analyses/Quadratic_LMER_Results_", celltype2, ".csv"),
              row.names = FALSE)
    
    # Print significant quadratic relationships
    cat("\n=== SIGNIFICANT QUADRATIC RELATIONSHIPS (p < 0.05) for", celltype, "===\n")
    sig_quad <- summary_table %>%
      filter(quadratic_significant) %>%
      select(pathway, outcome, quad_term_estimate, quad_term_p, best_model, delta_AIC_quad_vs_linear)
    
    if(nrow(sig_quad) > 0) {
      print(sig_quad)
      
      # Visualize significant quadratic relationships
      for(i in 1:nrow(sig_quad)) {
        
        pw <- sig_quad$pathway[i]
        oc <- sig_quad$outcome[i]
        
        # Get the result
        result <- celltype_results[[paste(oc, pw, sep = "_")]]
        
        if(is.null(result)) next
        
        # Create visualization
        plot_data <- result$data
        
        # Generate predictions from best model
        pred_grid <- expand.grid(
          pathway_centered = seq(min(plot_data$pathway_centered), 
                                 max(plot_data$pathway_centered), 
                                 length.out = 100),
          group = unique(plot_data$group)
        )
        pred_grid$pathway_squared <- pred_grid$pathway_centered^2
        
        if(sig_quad$best_model[i] == "Quadratic_Interaction") {
          pred_grid$pred <- predict(result$model_quad_int, newdata = pred_grid, re.form = NA)
        } else {
          # For non-interaction model, just use one group for prediction
          pred_grid_single <- pred_grid %>% filter(group == unique(plot_data$group)[1])
          pred_grid_single$pred <- predict(result$model_quad, newdata = pred_grid_single, re.form = NA)
          # Replicate for both groups
          pred_grid <- pred_grid_single %>%
            bind_rows(pred_grid_single %>% mutate(group = unique(plot_data$group)[2]))
        }
        
        # Create plot
        pathway_clean <- gsub("_score1|1$", "", pw)
        pathway_clean <- gsub("_", " ", pathway_clean)
        outcome_clean <- case_when(
          oc == "avg_c_k2" ~ "Cortical K2",
          oc == "avg_c_f" ~ "Cortical F",
          oc == "avg_c_k2_f" ~ "Cortical K2/F",
          TRUE ~ oc
        )
        
        p <- ggplot(plot_data, aes(x = pathway_centered, y = .data[[oc]])) +
          geom_point(aes(color = group), alpha = 0.5, size = 2) +
          geom_line(data = pred_grid, aes(y = pred, color = group, group = group), 
                    linewidth = 1.5) +
          scale_color_manual(values = c("Lean_Control" = "#3498db", "Type_2_Diabetes" = "#e74c3c"),
                             labels = c("Lean", "T2D")) +
          labs(
            title = paste0("Quadratic Relationship: ", pathway_clean, " vs ", outcome_clean),
            subtitle = paste0(celltype, " | p=", sprintf("%.4f", sig_quad$quad_term_p[i]), 
                              " | ", sig_quad$best_model[i]),
            x = paste0(pathway_clean, " (centered)"),
            y = outcome_clean,
            color = "Group"
          ) +
          theme_minimal(base_size = 12) +
          theme(
            plot.title = element_text(face = "bold", hjust = 0.5),
            plot.subtitle = element_text(hjust = 0.5, color = "grey30"),
            legend.position = "bottom"
          )
        
        pathway_file <- str_replace_all(pathway_clean, " ", "_")
        outcome2 <- str_replace_all(oc, "_", "")
        
        ggsave(
          paste0("C:/Users/netio/Documents/UofW/Rockies/Rockies_updates_9.26.25/new_analyses/Quadratic_", 
                 celltype2, "_", pathway_file, "_", outcome2, ".png"),
          plot = p,
          width = 8,
          height = 6,
          dpi = 300
        )
      }
    } else {
      cat("  No significant quadratic relationships found.\n")
    }
  } else {
    cat("  No valid models for", celltype, "\n")
  }
}

cat("\n=== Analysis complete! ===\n")
cat("Summary tables saved per cell type\n")
cat("Plots created for all significant quadratic relationships\n")






################# Correlations between module scores 


library(dplyr)
library(ggplot2)
library(corrplot)
library(lme4)
library(lmerTest)
library(stringr)

# Create output directory if it doesn't exist
output_dir <- "C:/Users/netio/Documents/UofW/Rockies/Rockies_updates_9.26.25/new_analyses/Module_Correlations/"
if(!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Load module scores
module_scores <- data.table::fread('C:/Users/netio/Documents/UofW/Rockies/Rockies_updates_9.26.25/module_scores/HALLMARK_GO_pathways_modulescores.txt')

# Define all pathways
all_pathways <- c("TCA_score1", "OxPhos_score1", "Response_to_Insulin1", 
                  "Insulin_Receptor_Signaling1", "Neg_Reg_Insulin_Signaling1", 
                  "Pos_Reg_Insulin_Signaling1", "IGF_Receptor_Binding1", 
                  "OxPhos_Hallmark1", "Fatty_Acid_Metabolism1", "Glycolysis1",
                  "mTORC1_Signaling1", "Adipogenesis1", "Peroxisome1", 
                  "Inflammatory_Response1", "IFN_Gamma_Response1", 
                  "IFN_Alpha_Response1", "IL6_JAK_STAT31", "TNFa_NF_kB1",
                  "Complement1", "Allograft_Rejection1", "Hypoxia1", 
                  "ROS_Pathway1", "PI3K_AKT_mTOR1")

# Define cell types
celltypes <- c('All', 'PT', 'TAL', 'EC', 'PT-S1/S2', 'PT-S3', 'aPT')

# Function to calculate mixed-effects correlations (cell-level, accounting for person)
calculate_mixed_correlations <- function(data, pathways) {
  
  n_pathways <- length(pathways)
  cor_matrix <- matrix(NA, n_pathways, n_pathways)
  p_matrix <- matrix(NA, n_pathways, n_pathways)
  rownames(cor_matrix) <- colnames(cor_matrix) <- pathways
  rownames(p_matrix) <- colnames(p_matrix) <- pathways
  
  cat("    Computing", n_pathways * (n_pathways - 1) / 2, "pairwise correlations...\n")
  
  counter <- 0
  total <- n_pathways * (n_pathways - 1) / 2
  
  for(i in 1:n_pathways) {
    for(j in i:n_pathways) {
      
      if(i == j) {
        cor_matrix[i, j] <- 1
        p_matrix[i, j] <- 0
        next
      }
      
      counter <- counter + 1
      if(counter %% 50 == 0) {
        cat("      Progress:", counter, "/", total, "\n")
      }
      
      tryCatch({
        # Prepare data - remove missing values
        test_data <- data[!is.na(data[[pathways[i]]]) & !is.na(data[[pathways[j]]]), ]
        
        # Need at least 20 observations
        if(nrow(test_data) < 20) {
          cor_matrix[i, j] <- cor_matrix[j, i] <- NA
          p_matrix[i, j] <- p_matrix[j, i] <- NA
          next
        }
        
        # Standardize variables (within this subset)
        test_data$x_std <- scale(test_data[[pathways[i]]])[,1]
        test_data$y_std <- scale(test_data[[pathways[j]]])[,1]
        
        # Fit mixed model: y ~ x + (1|person)
        model <- lmer(y_std ~ x_std + (1|record_id), data = test_data, REML = FALSE)
        
        # Extract fixed effect (standardized regression = correlation accounting for clustering)
        model_summary <- summary(model)
        cor_val <- fixef(model)["x_std"]
        p_val <- coef(model_summary)[2, "Pr(>|t|)"]
        
        cor_matrix[i, j] <- cor_matrix[j, i] <- cor_val
        p_matrix[i, j] <- p_matrix[j, i] <- p_val
        
      }, error = function(e) {
        cor_matrix[i, j] <- cor_matrix[j, i] <- NA
        p_matrix[i, j] <- p_matrix[j, i] <- NA
      })
    }
  }
  
  return(list(
    cor_matrix = cor_matrix,
    p_matrix = p_matrix,
    n = nrow(data)
  ))
}

# Function to create correlation plot
create_correlation_plot <- function(cor_matrix, p_matrix, celltype, group_filter = NULL) {
  
  # Clean pathway names
  rownames(cor_matrix) <- gsub("_score1|1$", "", rownames(cor_matrix))
  colnames(cor_matrix) <- gsub("_score1|1$", "", colnames(cor_matrix))
  rownames(cor_matrix) <- gsub("_", " ", rownames(cor_matrix))
  colnames(cor_matrix) <- gsub("_", " ", colnames(cor_matrix))
  
  if(!is.null(p_matrix)) {
    rownames(p_matrix) <- rownames(cor_matrix)
    colnames(p_matrix) <- colnames(cor_matrix)
  }
  
  title_suffix <- ifelse(!is.null(group_filter), 
                         paste0(" - ", gsub("_", " ", group_filter)),
                         " - All Participants")
  
  # Create corrplot
  corrplot(
    cor_matrix,
    method = "color",
    type = "upper",
    order = "hclust",
    p.mat = p_matrix,
    sig.level = 0.05,
    insig = "blank",
    tl.col = "black",
    tl.srt = 45,
    tl.cex = 0.7,
    cl.cex = 0.8,
    title = paste0("Module Score Correlations (Mixed-Effects): ", celltype, title_suffix),
    mar = c(0, 0, 2, 0)
  )
}

# Loop through cell types
for(celltype in celltypes) {
  
  cat("\n=== Processing cell type:", celltype, "===\n")
  
  # Filter by cell type
  if(celltype == 'All'){
    filtered_data <- module_scores
  } else if(celltype %in% c('PT', 'TAL', 'EC')){
    filtered_data <- module_scores %>% filter(celltype2 == celltype)
  } else {
    filtered_data <- module_scores %>% filter(KPMP_celltype == celltype)
  }
  
  # Check if enough data
  if(nrow(filtered_data) < 50) {
    cat("  Insufficient data for", celltype, "\n")
    next
  }
  
  # Clean celltype name for filename
  celltype2 <- str_replace_all(celltype, "/", "_")
  celltype2 <- str_replace_all(celltype2, "-", "_")
  
  # Calculate correlations - All participants (using all cells, accounting for person)
  cat("  Calculating mixed-effects correlations for all participants...\n")
  cor_result_all <- calculate_mixed_correlations(filtered_data, all_pathways)
  
  # Save correlation matrix
  cor_matrix_clean <- cor_result_all$cor_matrix
  rownames(cor_matrix_clean) <- gsub("_score1|1$", "", rownames(cor_matrix_clean))
  colnames(cor_matrix_clean) <- gsub("_score1|1$", "", colnames(cor_matrix_clean))
  
  write.csv(
    cor_matrix_clean,
    paste0(output_dir, "Correlation_Matrix_MixedEffects_", celltype2, "_All.csv")
  )
  
  # Save p-value matrix
  p_matrix_clean <- cor_result_all$p_matrix
  rownames(p_matrix_clean) <- rownames(cor_matrix_clean)
  colnames(p_matrix_clean) <- colnames(cor_matrix_clean)
  
  write.csv(
    p_matrix_clean,
    paste0(output_dir, "Correlation_PValues_MixedEffects_", celltype2, "_All.csv")
  )
  
  # Create plot - All participants
  png(
    paste0(output_dir, "Correlation_Plot_MixedEffects_", celltype2, "_All.png"),
    width = 14, height = 14, units = "in", res = 300
  )
  create_correlation_plot(cor_result_all$cor_matrix, cor_result_all$p_matrix, celltype, group_filter = NULL)
  dev.off()
  
  # Calculate correlations by group
  for(group in c("Lean_Control", "Type_2_Diabetes")) {
    
    cat("  Calculating correlations for", group, "...\n")
    
    group_data <- filtered_data %>% filter(group == !!group)
    
    if(nrow(group_data) < 50) {
      cat("    Insufficient data for", group, "\n")
      next
    }
    
    cor_result_group <- calculate_mixed_correlations(group_data, all_pathways)
    
    group_clean <- gsub("_", "", group)
    
    # Create plot by group
    png(
      paste0(output_dir, "Correlation_Plot_MixedEffects_", celltype2, "_", group_clean, ".png"),
      width = 14, height = 14, units = "in", res = 300
    )
    create_correlation_plot(cor_result_group$cor_matrix, cor_result_group$p_matrix, celltype, group_filter = group)
    dev.off()
    
    # Save group-specific matrices
    write.csv(
      cor_result_group$cor_matrix,
      paste0(output_dir, "Correlation_Matrix_MixedEffects_", celltype2, "_", group_clean, ".csv")
    )
    
    write.csv(
      cor_result_group$p_matrix,
      paste0(output_dir, "Correlation_PValues_MixedEffects_", celltype2, "_", group_clean, ".csv")
    )
  }
  
  # Create difference heatmap (if both groups have data)
  cat("  Creating difference heatmap...\n")
  
  lean_data <- filtered_data %>% filter(group == "Lean_Control")
  t2d_data <- filtered_data %>% filter(group == "Type_2_Diabetes")
  
  if(nrow(lean_data) >= 50 && nrow(t2d_data) >= 50) {
    
    cor_lean <- calculate_mixed_correlations(lean_data, all_pathways)
    cor_t2d <- calculate_mixed_correlations(t2d_data, all_pathways)
    
    # Calculate difference
    cor_diff <- cor_t2d$cor_matrix - cor_lean$cor_matrix
    
    # Clean names
    rownames(cor_diff) <- gsub("_score1|1$", "", rownames(cor_diff))
    colnames(cor_diff) <- gsub("_score1|1$", "", colnames(cor_diff))
    rownames(cor_diff) <- gsub("_", " ", rownames(cor_diff))
    colnames(cor_diff) <- gsub("_", " ", colnames(cor_diff))
    
    # Save difference matrix
    write.csv(
      cor_diff,
      paste0(output_dir, "Correlation_Difference_MixedEffects_T2D_vs_Lean_", celltype2, ".csv")
    )
    
    # Plot difference
    png(
      paste0(output_dir, "Correlation_Difference_MixedEffects_", celltype2, ".png"),
      width = 14, height = 14, units = "in", res = 300
    )
    
    corrplot(
      cor_diff,
      method = "color",
      type = "upper",
      is.corr = FALSE,
      col = colorRampPalette(c("#3498db", "white", "#e74c3c"))(200),
      tl.col = "black",
      tl.srt = 45,
      tl.cex = 0.7,
      cl.cex = 0.8,
      title = paste0("Correlation Difference (T2D - Lean): ", celltype),
      mar = c(0, 0, 2, 0),
      addCoef.col = "black",
      number.cex = 0.5
    )
    
    dev.off()
  }
  
  cat("  Completed", celltype, "\n")
}

cat("\n=== Analysis complete! ===\n")
cat("All files saved to:", output_dir, "\n")

















