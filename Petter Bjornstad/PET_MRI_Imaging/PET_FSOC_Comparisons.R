#PET_FSOC_Comparisons.R




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

















#Step 1: Hypermetabolism in T2D vs. HC (PET Analysis)



#harmonized_data <- read.csv("C:/Users/netio/Documents/Harmonized_data/harmonized_dataset.csv", na = '')




harmonized_data <- read.csv("C:/Users/netio/OneDrive - UW/Laura Pyle's files - Biostatistics Core Shared Drive/Data Harmonization/Data Clean/harmonized_dataset.csv", na = '')

#date_of_screen
#screen_date

dat <- harmonized_data %>% dplyr::select(-dob) %>% 
  arrange(date_of_screen) %>% 
  dplyr::summarise(across(where(negate(is.numeric)), ~ ifelse(all(is.na(.x)), NA_character_, first(na.omit(.x)))),
                   across(where(is.numeric), ~ ifelse(all(is.na(.x)), NA_real_, first(na.omit(.x)))),
                   .by = c(record_id, visit))


PET_avg <- function(data){
  tmp_df <- data %>% dplyr::select(lc_k2_wo_cyst_vw, rc_k2_wo_cyst_vw, lm_k2_wo_cyst_vw, rm_k2_wo_cyst_vw,
                                   lc_f, rc_f, lm_f, rm_f)
  avg_c_k2 <- tmp_df %>%
    dplyr::select(lc_k2_wo_cyst_vw, rc_k2_wo_cyst_vw) %>% rowMeans(na.rm=T)
  
  avg_m_k2 <- tmp_df %>% 
    dplyr::select(lm_k2_wo_cyst_vw, rm_k2_wo_cyst_vw) %>% rowMeans(na.rm=T)
  
  avg_c_f <- tmp_df %>% 
    dplyr::select(lc_f, rc_f) %>% rowMeans(na.rm=T)
  
  avg_m_f <- tmp_df %>% 
    dplyr::select(lm_f, rm_f) %>% rowMeans(na.rm=T)
  
  avg_c_k2_f <- avg_c_k2 / avg_c_f
  
  avg_m_k2_f <- avg_m_k2/ avg_m_f
  
  results <- bind_cols(avg_c_k2, avg_m_k2, avg_c_f, avg_m_f, 
                       avg_c_k2_f, avg_m_k2_f) %>% as.data.frame()
  names(results) <- c('avg_c_k2_vw', 'avg_m_k2_vw', 'avg_c_f_vw', 'avg_m_f_vw', 
                      'avg_c_k2_f_vw', 'avg_m_k2_f_vw')
  
  return(results)
  
}


tmp_results_vw <- PET_avg(dat)


dat_results <- dat




PET_avg <- function(data){
  tmp_df <- data %>% dplyr::select(lc_k2, rc_k2, lm_k2, rm_k2,
                                   lc_f, rc_f, lm_f, rm_f)
  avg_c_k2 <- tmp_df %>%
    dplyr::select(lc_k2, rc_k2) %>% rowMeans(na.rm=T)
  
  avg_m_k2 <- tmp_df %>% 
    dplyr::select(lm_k2, rm_k2) %>% rowMeans(na.rm=T)
  
  avg_c_f <- tmp_df %>% 
    dplyr::select(lc_f, rc_f) %>% rowMeans(na.rm=T)
  
  avg_m_f <- tmp_df %>% 
    dplyr::select(lm_f, rm_f) %>% rowMeans(na.rm=T)
  
  avg_c_k2_f <- avg_c_k2 / avg_c_f
  
  avg_m_k2_f <- avg_m_k2 / avg_m_f
  
  results <- bind_cols(avg_c_k2, avg_m_k2, avg_c_f, avg_m_f, 
                       avg_c_k2_f, avg_m_k2_f) %>% as.data.frame()
  names(results) <- c('avg_c_k2', 'avg_m_k2', 'avg_c_f', 'avg_m_f', 
                      'avg_c_k2_f', 'avg_m_k2_f')
  
  return(results)
  
}


tmp_results <- PET_avg(dat)



dat_results <- dat %>% 
  dplyr::select(-any_of(c("avg_c_f", "avg_c_k2", "avg_m_f", "avg_m_k2",
                          "avg_c_k2_f", "avg_m_k2_f"))) %>%
  bind_cols(tmp_results, tmp_results_vw)


dat_results <- dat_results %>% filter(!is.na(avg_c_k2))




# ============================================================================
# PET ANALYSIS: Hypermetabolism Across Disease Groups
# Boxplots, Sex-stratified boxplots, and Heatmaps
# Matching FSOC analysis style
# ============================================================================

# Additional libraries needed
library(car)
library(rstatix)
library(ggpubr)
library(ggsignif)

# Set output directory
OUTPUT_DIR <- "C:/Users/netio/Documents/UofW/Projects/Imaging_Shivani/PET_Analysis"
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ============================================================================
# DATA PREPARATION
# ============================================================================

# Define PET endpoints (non-VW versions)
pet_endpoints <- c("avg_c_k2", "avg_m_k2", "avg_c_f", "avg_m_f", 
                   "avg_c_k2_f", "avg_m_k2_f")

# Define nice labels for PET variables
pet_labels <- c(
  "avg_c_k2" = "Cortical K2",
  "avg_m_k2" = "Medullary K2",
  "avg_c_f" = "Cortical F (Perfusion)",
  "avg_m_f" = "Medullary F (Perfusion)",
  "avg_c_k2_f" = "Cortical K2/F Ratio",
  "avg_m_k2_f" = "Medullary K2/F Ratio"
)

# Group levels and colors - MATCHING FSOC ANALYSIS
group_levels <- c("Lean Control", "Obese Control", "Type 1 Diabetes", 
                  "Type 2 Diabetes", "PKD")

# Colors matching your FSOC analysis exactly
group_colors <- c("Lean Control" = "#3182BD",
                  "Obese Control" = "#9ECAE1",
                  "Type 1 Diabetes" = "#FDAE6B",
                  "Type 2 Diabetes" = "#E6550D",
                  "PKD" = "#9E9AC8")

# Prepare the data - set up group and sex factors
dat_results <- dat_results %>%
  mutate(group = factor(group, levels = group_levels),
         sex = factor(sex, levels = c("Female", "Male")))

# ============================================================================
# 1. BOXPLOTS BY GROUP WITH PAIRWISE P-VALUES (simplified, fast version)
# ============================================================================

plot_pet_by_group <- function(dat, endpoint, label) {
  
  dat_plot <- dat %>% 
    filter(!is.na(group), !is.na(.data[[endpoint]])) %>%
    mutate(group = factor(group, levels = group_levels))
  
  # Run ANOVA with covariates for adjusted analysis
  formula_str <- paste(endpoint, "~ group + age + sex + weight")
  model <- tryCatch({
    lm(as.formula(formula_str), data = dat_plot)
  }, error = function(e) {
    message("Model failed for ", endpoint, ": ", e$message)
    return(NULL)
  })
  
  # Get pairwise comparisons using emmeans (adjusted)
  pairs_df <- NULL
  if (!is.null(model)) {
    emm <- emmeans(model, ~ group)
    pairs_result <- pairs(emm, adjust = "tukey")
    pairs_df <- as.data.frame(pairs_result)
    
    # Print pairwise results
    cat("\n--- Pairwise comparisons for", label, "(adjusted) ---\n")
    print(pairs_df[, c("contrast", "estimate", "SE", "p.value")])
  }
  
  # Calculate unadjusted pairwise comparisons (Wilcoxon)
  stat_test <- dat_plot %>%
    wilcox_test(as.formula(paste(endpoint, "~ group"))) %>%
    mutate(p_label = case_when(
      p < 0.001 ~ "p<0.001",
      p < 0.01 ~ paste0("p=", formatC(p, format = "f", digits = 3)),
      TRUE ~ paste0("p=", formatC(p, format = "f", digits = 2))
    ))
  
  # Get y-axis range for positioning
  y_max <- max(dat_plot[[endpoint]], na.rm = TRUE)
  y_min <- min(dat_plot[[endpoint]], na.rm = TRUE)
  y_range <- y_max - y_min
  
  # Create boxplot (no p-value brackets for speed - they're printed to console)
  p <- ggplot(dat_plot, aes(x = group, y = .data[[endpoint]], fill = group)) +
    geom_boxplot(width = 0.7, outlier.shape = NA, alpha = 0.8) +
    geom_jitter(width = 0.2, alpha = 0.5, size = 1.5) +
    scale_fill_manual(values = group_colors) +
    labs(title = label,
         subtitle = "Raw pairwise Wilcoxon p-values printed to console",
         y = label, x = "") +
    theme_bw() +
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      plot.title = element_text(face = "bold", size = 12),
      plot.subtitle = element_text(size = 8, color = "gray50"),
      panel.grid.minor = element_blank()
    )
  
  return(list(plot = p, model = model, pairwise = pairs_df, stat_test = stat_test))
}

# Generate all group boxplots
generate_group_boxplots <- function(dat) {
  results <- list()
  plots <- list()
  
  for (endpoint in pet_endpoints) {
    label <- pet_labels[endpoint]
    result <- tryCatch({
      plot_pet_by_group(dat, endpoint, label)
    }, error = function(e) {
      message("Failed for ", endpoint, ": ", e$message)
      return(NULL)
    })
    
    if (!is.null(result)) {
      results[[endpoint]] <- result
      plots[[endpoint]] <- result$plot
    }
  }
  
  # Combine plots
  combined_plot <- cowplot::plot_grid(plotlist = plots, ncol = 3, nrow = 2)
  
  return(list(individual = results, combined = combined_plot))
}

# ============================================================================
# 2. BOXPLOTS BY GROUP SPLIT BY SEX (raw pairwise p-values)
# ============================================================================

plot_pet_by_group_sex <- function(dat, endpoint, label) {
  
  dat_plot <- dat %>% 
    filter(!is.na(group), !is.na(.data[[endpoint]]), !is.na(sex)) %>%
    mutate(
      group = factor(group, levels = group_levels),
      sex = factor(sex, levels = c("Female", "Male"))
    ) %>%
    droplevels()
  
  # Check we have data
  if (nrow(dat_plot) < 5) {
    message("Not enough data for ", endpoint)
    return(NULL)
  }
  
  # Calculate RAW pairwise comparisons within each group (between sexes)
  # No adjustment - raw Wilcoxon p-values
  stat_test <- NULL
  tryCatch({
    stat_test <- dat_plot %>%
      group_by(group) %>%
      filter(n_distinct(sex) == 2) %>%  # Only groups with both sexes
      wilcox_test(as.formula(paste(endpoint, "~ sex"))) %>%
      ungroup() %>%
      mutate(p_label = case_when(
        p < 0.001 ~ "p<0.001",
        p < 0.01 ~ paste0("p=", formatC(p, format = "f", digits = 3)),
        TRUE ~ paste0("p=", formatC(p, format = "f", digits = 2))
      ))
  }, error = function(e) {
    message("Wilcoxon test failed: ", e$message)
  })
  
  # Get y-axis range for positioning
  y_max <- max(dat_plot[[endpoint]], na.rm = TRUE)
  y_min <- min(dat_plot[[endpoint]], na.rm = TRUE)
  y_range <- y_max - y_min
  
  # Create boxplot - NO subtitle
  p <- ggplot(dat_plot, aes(x = group, y = .data[[endpoint]], fill = sex)) +
    geom_boxplot(width = 0.7, alpha = 0.8,
                 position = position_dodge(width = 0.75)) +
    scale_fill_manual(values = c("Female" = "#E78AC3", "Male" = "#66C2A5"),
                      name = "Sex") +
    labs(title = label,
         y = label, x = "Disease Group") +
    theme_bw() +
    theme(
      legend.position = "bottom",
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      plot.title = element_text(face = "bold", size = 12),
      panel.grid.minor = element_blank()
    )
  
  # Add raw p-values above each group if we have them
  if (!is.null(stat_test) && nrow(stat_test) > 0) {
    groups_in_plot <- levels(dat_plot$group)
    for (i in seq_along(groups_in_plot)) {
      grp <- groups_in_plot[i]
      p_row <- stat_test %>% filter(group == grp)
      if (nrow(p_row) > 0) {
        p <- p + annotate("text", 
                          x = i, 
                          y = y_max + 0.08 * y_range,
                          label = p_row$p_label[1],
                          size = 2.5)
      }
    }
    # Expand y-axis to fit p-values
    p <- p + coord_cartesian(ylim = c(y_min - 0.02 * y_range, y_max + 0.18 * y_range))
  }
  
  return(list(plot = p, stat_test = stat_test))
}

# Generate all sex-stratified boxplots
generate_sex_boxplots <- function(dat) {
  results <- list()
  plots <- list()
  
  cat("Starting sex-stratified boxplots...\n")
  
  for (endpoint in pet_endpoints) {
    cat("  Processing:", endpoint, "\n")
    label <- pet_labels[endpoint]
    result <- tryCatch({
      plot_pet_by_group_sex(dat, endpoint, label)
    }, error = function(e) {
      cat("    ERROR for ", endpoint, ": ", e$message, "\n")
      return(NULL)
    })
    
    if (!is.null(result)) {
      results[[endpoint]] <- result
      plots[[endpoint]] <- result$plot
      cat("    Success!\n")
    } else {
      cat("    Returned NULL\n")
    }
  }
  
  if (length(plots) == 0) {
    cat("No sex-stratified plots were generated!\n")
    return(NULL)
  }
  
  cat("Combining", length(plots), "plots...\n")
  combined_plot <- cowplot::plot_grid(plotlist = plots, ncol = 3, nrow = 2)
  
  return(list(individual = results, combined = combined_plot))
}

# ============================================================================
# 3. HEATMAPS FOR DXA ASSOCIATIONS
# ============================================================================

analyze_dxa_associations_pet <- function(dat) {
  
  # DXA variables - UPDATE THESE to match your actual column names
  dxa_vars <- c("dexa_fat_kg", "dexa_trunk_kg", "dexa_lean_kg", "dexa_body_fat")
  
  # Nice labels for DXA
  dxa_labels <- c(
    "dexa_fat_kg" = "Fat Mass (kg)",
    "dexa_trunk_kg" = "Trunk Fat (kg)",
    "dexa_lean_kg" = "Lean Mass (kg)",
    "dexa_body_fat" = "Body Fat (%)"
  )
  
  # Filter to only existing variables
  existing_dxa <- dxa_vars[dxa_vars %in% names(dat)]
  if (length(existing_dxa) == 0) {
    message("No DXA variables found! Check variable names.")
    return(NULL)
  }
  
  cat("DXA variables found:", paste(existing_dxa, collapse = ", "), "\n")
  
  # Create matrices for coefficients and p-values
  coef_matrix <- matrix(NA, nrow = length(pet_endpoints), ncol = length(existing_dxa))
  p_matrix <- matrix(NA, nrow = length(pet_endpoints), ncol = length(existing_dxa))
  
  rownames(coef_matrix) <- pet_labels[pet_endpoints]
  colnames(coef_matrix) <- dxa_labels[existing_dxa]
  rownames(p_matrix) <- pet_labels[pet_endpoints]
  colnames(p_matrix) <- dxa_labels[existing_dxa]
  
  # Run GLMs adjusted for age, sex, disease group
  for (i in seq_along(pet_endpoints)) {
    for (j in seq_along(existing_dxa)) {
      formula_str <- paste(pet_endpoints[i], "~", existing_dxa[j], "+ age + sex + group")
      
      tryCatch({
        model <- lm(as.formula(formula_str), data = dat)
        coef_summary <- summary(model)$coefficients
        
        dxa_row <- which(rownames(coef_summary) == existing_dxa[j])
        if (length(dxa_row) > 0) {
          coef_matrix[i, j] <- coef_summary[dxa_row, "Estimate"]
          p_matrix[i, j] <- coef_summary[dxa_row, "Pr(>|t|)"]
        }
      }, error = function(e) {
        message(paste("Error fitting model for", pet_endpoints[i], "~", existing_dxa[j]))
      })
    }
  }
  
  return(list(coefficient_matrix = coef_matrix, p_value_matrix = p_matrix))
}

plot_dxa_heatmap_pet <- function(dxa_results) {
  
  if (is.null(dxa_results)) return(NULL)
  
  # Create significance labels
  p_mat <- dxa_results$p_value_matrix
  sig_labels <- matrix("", nrow = nrow(p_mat), ncol = ncol(p_mat))
  sig_labels[p_mat < 0.001] <- "***"
  sig_labels[p_mat >= 0.001 & p_mat < 0.01] <- "**"
  sig_labels[p_mat >= 0.01 & p_mat < 0.05] <- "*"
  
  # Create heatmap (removed angle_col which was causing error)
  pheatmap(dxa_results$coefficient_matrix,
           display_numbers = sig_labels,
           cluster_rows = FALSE,
           cluster_cols = FALSE,
           color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
           main = "PET Associations with DXA Parameters\n(Adjusted for age, sex, disease group)",
           fontsize = 10,
           fontsize_number = 14)
}

# ============================================================================
# 4. HEATMAPS FOR CLINICAL ASSOCIATIONS
# ============================================================================

analyze_clinical_associations_pet <- function(dat) {
  
  # Clinical variables - UPDATE to match your column names
  clinical_vars <- c("diabetes_duration", "hba1c", "eGFR_CKD_epi", "gfr_raw_plasma", "gfr_bsa_plasma",
                     "erpf_raw_plasma", "erpf_bsa_plasma", "ff",
                     "acr_u", "sbp", "dbp")
  
  clinical_labels <- c(
    "diabetes_duration" = "DM Duration",
    "hba1c" = "HbA1c",
    "eGFR_CKD_epi" = "eGFR (CKD-EPI)",
    "gfr_raw_plasma" = "mGFR (raw)",
    "gfr_bsa_plasma" = "mGFR (BSA)",
    "erpf_raw_plasma" = "ERPF (raw)",
    "erpf_bsa_plasma" = "ERPF (BSA)",
    "ff" = "FF",
    "acr_u" = "log(UACR)",
    "sbp" = "SBP",
    "dbp" = "DBP"
  )
  
  # Filter to only existing variables
  existing_vars <- clinical_vars[clinical_vars %in% names(dat)]
  if (length(existing_vars) == 0) {
    message("No clinical variables found! Check variable names.")
    return(NULL)
  }
  
  cat("Clinical variables found:", paste(existing_vars, collapse = ", "), "\n")
  existing_labels <- clinical_labels[existing_vars]
  
  coef_matrix <- matrix(NA, nrow = length(pet_endpoints), ncol = length(existing_vars))
  p_matrix <- matrix(NA, nrow = length(pet_endpoints), ncol = length(existing_vars))
  
  rownames(coef_matrix) <- pet_labels[pet_endpoints]
  colnames(coef_matrix) <- existing_labels
  rownames(p_matrix) <- pet_labels[pet_endpoints]
  colnames(p_matrix) <- existing_labels
  
  for (i in seq_along(pet_endpoints)) {
    for (j in seq_along(existing_vars)) {
      
      # For UACR, use log transformation
      if (existing_vars[j] == "acr_u") {
        dat$log_uacr <- log(dat$acr_u + 1)
        formula_str <- paste(pet_endpoints[i], "~ log_uacr + age + sex + group")
        var_name <- "log_uacr"
      } else {
        formula_str <- paste(pet_endpoints[i], "~", existing_vars[j], "+ age + sex + group")
        var_name <- existing_vars[j]
      }
      
      tryCatch({
        model <- lm(as.formula(formula_str), data = dat)
        coef_summary <- summary(model)$coefficients
        
        var_row <- which(rownames(coef_summary) == var_name)
        if (length(var_row) > 0) {
          coef_matrix[i, j] <- coef_summary[var_row, "Estimate"]
          p_matrix[i, j] <- coef_summary[var_row, "Pr(>|t|)"]
        }
      }, error = function(e) {
        message(paste("Error fitting model for", pet_endpoints[i], "~", existing_vars[j]))
      })
    }
  }
  
  return(list(coefficient_matrix = coef_matrix, p_value_matrix = p_matrix))
}

plot_clinical_heatmap_pet <- function(clinical_results) {
  
  if (is.null(clinical_results)) return(NULL)
  
  p_mat <- clinical_results$p_value_matrix
  sig_labels <- matrix("", nrow = nrow(p_mat), ncol = ncol(p_mat))
  sig_labels[p_mat < 0.001] <- "***"
  sig_labels[p_mat >= 0.001 & p_mat < 0.01] <- "**"
  sig_labels[p_mat >= 0.01 & p_mat < 0.05] <- "*"
  
  # Create heatmap (removed angle_col which was causing error)
  pheatmap(clinical_results$coefficient_matrix,
           display_numbers = sig_labels,
           cluster_rows = FALSE,
           cluster_cols = FALSE,
           color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
           main = "PET Associations with Clinical Parameters\n(Adjusted for age, sex, disease group)",
           fontsize = 10,
           fontsize_number = 14)
}

# ============================================================================
# 5. MAIN WORKFLOW
# ============================================================================

run_pet_analysis <- function(dat) {
  
  cat("Starting PET Analysis...\n")
  cat("Sample size:", nrow(dat), "\n\n")
  
  # Step 1: Group boxplots with pairwise p-values
  cat("=== Generating group boxplots with pairwise comparisons ===\n")
  group_results <- generate_group_boxplots(dat)
  
  # Step 2: Sex-stratified boxplots
  cat("\n=== Generating sex-stratified boxplots ===\n")
  sex_results <- tryCatch({
    generate_sex_boxplots(dat)
  }, error = function(e) {
    cat("Sex analysis failed:", e$message, "\n")
    return(NULL)
  })
  
  # Step 3: DXA heatmap
  cat("\n=== Analyzing DXA associations ===\n")
  dxa_results <- analyze_dxa_associations_pet(dat)
  
  # Step 4: Clinical heatmap
  cat("\n=== Analyzing clinical associations ===\n")
  clinical_results <- analyze_clinical_associations_pet(dat)
  
  # Save outputs
  cat("\n=== Saving outputs to:", OUTPUT_DIR, "===\n")
  
  # Save group boxplots
  pdf(file.path(OUTPUT_DIR, "pet_by_group.pdf"), width = 14, height = 10)
  print(group_results$combined)
  dev.off()
  cat("Saved: pet_by_group.pdf\n")
  
  # Save sex-stratified boxplots
  if (!is.null(sex_results)) {
    pdf(file.path(OUTPUT_DIR, "pet_by_group_sex.pdf"), width = 14, height = 10)
    print(sex_results$combined)
    dev.off()
    cat("Saved: pet_by_group_sex.pdf\n")
  }
  
  # Save DXA heatmap
  if (!is.null(dxa_results) && !all(is.na(dxa_results$coefficient_matrix))) {
    tryCatch({
      pdf(file.path(OUTPUT_DIR, "pet_dxa_heatmap.pdf"), width = 10, height = 8)
      plot_dxa_heatmap_pet(dxa_results)
      dev.off()
      cat("Saved: pet_dxa_heatmap.pdf\n")
    }, error = function(e) {
      cat("Error saving DXA heatmap:", e$message, "\n")
      dev.off()
    })
    
    write.csv(dxa_results$coefficient_matrix, 
              file.path(OUTPUT_DIR, "pet_dxa_coefficients.csv"), row.names = TRUE)
    write.csv(dxa_results$p_value_matrix, 
              file.path(OUTPUT_DIR, "pet_dxa_pvalues.csv"), row.names = TRUE)
    cat("Saved: pet_dxa_coefficients.csv and pet_dxa_pvalues.csv\n")
  } else {
    cat("DXA results are NULL or all NA - skipping heatmap\n")
  }
  
  # Save clinical heatmap
  if (!is.null(clinical_results) && !all(is.na(clinical_results$coefficient_matrix))) {
    tryCatch({
      pdf(file.path(OUTPUT_DIR, "pet_clinical_heatmap.pdf"), width = 12, height = 8)
      plot_clinical_heatmap_pet(clinical_results)
      dev.off()
      cat("Saved: pet_clinical_heatmap.pdf\n")
    }, error = function(e) {
      cat("Error saving clinical heatmap:", e$message, "\n")
      dev.off()
    })
    
    write.csv(clinical_results$coefficient_matrix, 
              file.path(OUTPUT_DIR, "pet_clinical_coefficients.csv"), row.names = TRUE)
    write.csv(clinical_results$p_value_matrix, 
              file.path(OUTPUT_DIR, "pet_clinical_pvalues.csv"), row.names = TRUE)
    cat("Saved: pet_clinical_coefficients.csv and pet_clinical_pvalues.csv\n")
  } else {
    cat("Clinical results are NULL or all NA - skipping heatmap\n")
  }
  
  cat("\n=== Analysis complete! ===\n")
  
  return(list(
    group_boxplots = group_results,
    sex_boxplots = sex_results,
    dxa_associations = dxa_results,
    clinical_associations = clinical_results
  ))
}

# ============================================================================
# RUN THE ANALYSIS
# ============================================================================

# Check your data first
cat("Sample size:", nrow(dat_results), "\n")
cat("Groups:\n")
print(table(dat_results$group, useNA = "ifany"))
cat("\nSex:\n")
print(table(dat_results$sex, useNA = "ifany"))

# Run the full analysis
results <- run_pet_analysis(dat_results)

# View the boxplots interactively
print(results$group_boxplots$combined)
if (!is.null(results$sex_boxplots)) print(results$sex_boxplots$combined)

# View heatmaps interactively
if (!is.null(results$dxa_associations)) plot_dxa_heatmap_pet(results$dxa_associations)
if (!is.null(results$clinical_associations)) plot_clinical_heatmap_pet(results$clinical_associations)







library(gtsummary)
library(gt)

demographics <- dat_results %>%
  select(age, sex, race_ethnicity, bmi, hba1c, diabetes_duration, eGFR_CKD_epi, acr_u, group) %>%
  tbl_summary(
    by = group,
    type = list(
      age ~ "continuous",
      bmi ~ "continuous", 
      hba1c ~ "continuous",
      diabetes_duration ~ "continuous",
      eGFR_CKD_epi ~ "continuous",
      acr_u ~ "continuous",
      sex ~ "categorical",
      race_ethnicity ~ "categorical"
    ),
    statistic = list(
      all_continuous() ~ "{mean} ({sd})",
      all_categorical() ~ "{n} ({p}%)"
    ),
    digits = list(
      age ~ 1,
      bmi ~ 1,
      hba1c ~ 2,
      diabetes_duration ~ 1,
      eGFR_CKD_epi ~ 1,
      acr_u ~ 1,
      all_categorical() ~ c(0, 1)
    ),
    label = list(
      age ~ "Age, years",
      sex ~ "Sex", 
      race_ethnicity ~ "Race/Ethnicity",
      bmi ~ "BMI, kg/m²",
      hba1c ~ "HbA1c, %",
      diabetes_duration ~ "Diabetes Duration, years",
      eGFR_CKD_epi ~ "eGFR, mL/min/1.73m²",
      acr_u ~ "UACR, mg/g"
    ),
    missing_text = "Missing"
  ) %>%
  add_p(test = list(
    all_continuous() ~ "kruskal.test",
    all_categorical() ~ "chisq.test"
  )) %>%
  add_overall(col_label = "**Overall**\nN = {N}") %>%
  modify_header(label ~ "**Characteristic**") %>%
  modify_spanning_header(all_stat_cols() ~ "**Group**")

# View it
demographics

# Save it
demographics %>%
  as_gt() %>%
  gtsave(file.path(OUTPUT_DIR, "pet_demographics_table.png"), vwidth = 1400, vheight = 600)
















