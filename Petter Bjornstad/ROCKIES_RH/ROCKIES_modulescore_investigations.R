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


module_scores <- data.table::fread('C:/Users/netio/Documents/UofW/Rockies/Rockies_updates_9.26.25/module_scores/HALLMARK_GO_pathways_modulescores.txt')



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
  filter(p_adj_fdr < 0.05) %>%
  arrange(p_adj_fdr)

print(ks_significant)



# Test interaction between group and gene expression
library(lme4)

interaction_data <- combined_df %>%
  filter(!is.na(avg_c_k2f)) %>%
  select(record_id, group, OxPhos_score1, avg_c_k2f)

# Linear model with interaction
model <- lm(avg_c_k2f ~ OxPhos_score1 * group, data = interaction_data)
summary(model)

# Visualize
ggplot(interaction_data, aes(x = OxPhos_score1, y = avg_c_k2f, color = group)) +
  geom_point(size = 4) +
  geom_smooth(method = "lm", se = TRUE) +
  stat_cor(method = "spearman", aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~"))) +
  scale_color_manual(values = c("Lean_Control" = "blue", "Type_2_Diabetes" = "red")) +
  labs(title = "Do Lean and T2D Have Different Slopes?",
       subtitle = "Testing for interaction between group and gene expression",
       x = "OxPhos Module Score",
       y = "Cortical K2/F (voxel)") +
  theme_minimal()


















