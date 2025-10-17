############### Brain Data Analysis

library(purrr)
library(dplyr)
library(stringr)
library(tidyverse)
library(ggplot2)



harmonized_data <- read.csv("C:/Users/netio/OneDrive - UW/Laura Pyle's files - Biostatistics Core Shared Drive/Data Harmonization/Data Clean/harmonized_dataset.csv", na = '')

#date_of_screen
#screen_date

#dat <- harmonized_data %>% dplyr::select(-dob) %>% 
#  arrange(date_of_screen) %>% 
#  dplyr::summarise(across(where(negate(is.numeric)), ~ ifelse(all(is.na(.x)), NA_character_, first(na.omit(.x)))),
#                   across(where(is.numeric), ~ ifelse(all(is.na(.x)), NA_real_, first(na.omit(.x)))),
#                   .by = c(record_id, visit))



dat <- harmonized_data %>% dplyr::select(-dob) %>% 
  arrange(date_of_screen) %>% 
  dplyr::summarise(across(where(negate(is.numeric)), ~ ifelse(all(is.na(.x)), NA_character_, last(na.omit(.x)))),
                   across(where(is.numeric), ~ ifelse(all(is.na(.x)), NA_real_, mean(na.omit(.x), na.rm=T))),
                   .by = c(record_id, visit))



mri_ids <- c('CRC-10', 'CRC-11', 'CRC-12', 'CRC-13', 'CRC-26', 'CRC-39', 'CRC-46', 'CRC-51', 'CRC-53', 
             'CRC-55', 'CRC-58', 'CRC-60', 
             'RH2-01-O', 'RH2-03-O', 'RH2-08-T', 'RH2-10-L', 'RH2-11-O', 'RH2-13-O', 'RH2-16-O', 'RH2-17-L', 
             'RH2-18-O', 'RH2-19-T', 'RH2-22-T', 'RH2-24-L', 'RH2-27-L', 'RH2-28-L', 'RH2-29-L', 'RH2-33-L', 
             'RH2-34-O', 'RH2-35-T', 'RH2-38-T', 'RH2-39-O', 'RH2-41-T', 'RH2-42-T', 'RH2-43-T', 
             'RH2-44-T', 'RH2-45-T', 'RH2-48-T', 'RH2-49-T', 'RH2-50-L', 'RH2-52-T', 'RH2-53-T', 
             'RH2-55-T')

mri_ids_df <- data.frame(ID = mri_ids, file_id = c(
  'crc_10', 'crc_11', 'crc_12', 'crc_13', 'crc_26', 'crc_39', 'crc_46', 'crc_51', 'crc_53', 
  'crc_55', 'crc_58', 'crc_60', 
  'RH2_01_O', 'RH2_03_O', 'RH2_08_T', 'RH2_10_L', 'RH2_11_O', 'RH2_13_O', 'RH2_16_O', 'RH2_17_L', 
  'RH2_18_O', 'RH2_19_T', 'RH2_22_T', 'RH2_24_L', 'RH2_27_L', 'RH2_28_L', 'RH2_29_L', 'RH2_33_L', 
  'RH2_34_O', 'RH2_35_T', 'RH2_38_T', 'RH2_39_O', 'RH2_41_T', 'RH2_42_T', 'RH2_43_T', 
  'RH2_44_T', 'RH2_45_T', 'RH2_48_T', 'RH2_49_T', 'RH2_50_L', 'RH2_52_T', 'RH2_53_T', 
  'RH2_55_T'
))


small_dat <- dat %>% 
  filter(record_id %in% mri_ids)

#find missing
#mri_ids[which(!mri_ids %in% small_dat$record_id)]

#missing_dat <- dat %>% filter(rh2_id == 'RH2-38-O')


small_dat$group[which(small_dat$record_id == 'RH2-38-T')] <- 'Obese Control'
small_dat$record_id[which(small_dat$record_id == 'RH2-38-T')] <- 'RH2-38-O'





qx_var <- c("ab40_avg_conc","ab42_avg_conc","tau_avg_conc",
            "nfl_avg_conc","gfap_avg_conc","ptau_181_avg_conc","ptau_217_avg_conc")





small_dat <- small_dat %>% 
  dplyr::select(record_id, group, ab40_avg_conc, ab42_avg_conc, tau_avg_conc, 
                nfl_avg_conc, gfap_avg_conc, ptau_181_avg_conc, ptau_217_avg_conc)



t2d_ids <- small_dat$record_id[which(small_dat$group == 'Type 2 Diabetes')]
lc_ids <- small_dat$record_id[which(small_dat$group == 'Lean Control')]



small_dat <- small_dat %>% left_join(mri_ids_df, by = c('record_id'='ID'))





#Data 
subcortical <- data.table::fread("/Users/netio/Documents/UofW/Projects/Brain_fMRI_Analysis/t1_MRI/subcortical_volumes.txt")
names(subcortical)[1] <- 'file_id'
subcortical$file_id <- str_remove(subcortical$file_id, "_t1\\.?/$")

l_thick <- data.table::fread("/Users/netio/Documents/UofW/Projects/Brain_fMRI_Analysis/t1_MRI/lh_thickness.txt") 
names(l_thick)[1] <- 'file_id'
l_thick$file_id <- str_remove(l_thick$file_id, "_t1\\.?/$")

r_thick <- data.table::fread("/Users/netio/Documents/UofW/Projects/Brain_fMRI_Analysis/t1_MRI/rh_thickness.txt")
names(r_thick)[1] <- 'file_id'
r_thick$file_id <- str_remove(r_thick$file_id, "_t1\\.?/$")



t1_analysis <- small_dat %>% left_join(subcortical) %>% 
  left_join(l_thick, by = 'file_id') %>% 
  left_join(r_thick, by = 'file_id') %>% 
  mutate(group2 = ifelse(group == 'Type 2 Diabetes', 'Type 2 Diabetes', 'Control'))


ggplot(t1_analysis, aes(x= group, y = `Left-Putamen`))+geom_boxplot()






#### Plotting and Analyzing

setwd('/Users/netio/Documents/UofW/Projects/Brain_fMRI_Analysis/t1_MRI/')

library(dplyr)
library(tidyr)
library(ggplot2)
library(broom)
library(patchwork)

# Set your data
brain_data <- t1_analysis

# ===== Step 1: Identify left/right paired variables =====

# Get all column names
all_cols <- names(brain_data)

# Identify left hemisphere variables
lh_vars <- grep("^lh_|^Left-", all_cols, value = TRUE)

# Identify right hemisphere variables  
rh_vars <- grep("^rh_|^Right-", all_cols, value = TRUE)

# Create pairs (match lh_ with rh_, Left- with Right-)
lh_base <- gsub("^lh_|^Left-", "", lh_vars)
rh_base <- gsub("^rh_|^Right-", "", rh_vars)

# Find matching pairs
paired_regions <- intersect(lh_base, rh_base)

# ===== Step 2: Create combined (average) variables =====

for (region in paired_regions) {
  # Find the lh and rh column names
  lh_col <- lh_vars[grep(paste0("^lh_|^Left-", region), lh_vars)]
  rh_col <- rh_vars[grep(paste0("^rh_|^Right-", region), rh_vars)]
  
  if (length(lh_col) == 1 & length(rh_col) == 1) {
    # Create bilateral average
    brain_data[[paste0("bilateral_", region)]] <- 
      (brain_data[[lh_col]] + brain_data[[rh_col]]) / 2
  }
}

# ===== Step 3: Statistical testing for ALL variables (ROBUST) =====

# Get all numeric columns to test (exclude IDs and group)
test_vars <- names(brain_data)[sapply(brain_data, is.numeric)]
test_vars <- setdiff(test_vars, c("record_id", "file_id"))

# Initialize results dataframe
results <- data.frame()

# Test each variable
for (var in test_vars) {
  
  # Extract data for this variable
  test_data <- brain_data %>%
    select(group2, all_of(var)) %>%
    filter(!is.na(group2) & !is.na(.data[[var]]))
  
  # Skip if not enough data
  if (nrow(test_data) < 3 | length(unique(test_data$group2)) < 2) {
    next
  }
  
  # Get values for each group
  control_vals <- test_data[[var]][test_data$group2 == "Control"]
  t2d_vals <- test_data[[var]][test_data$group2 == "Type 2 Diabetes"]
  
  # Skip if not enough observations per group
  if (length(control_vals) < 2 | length(t2d_vals) < 2) {
    next
  }
  
  # Check if values are constant (no variance)
  control_has_variance <- length(unique(control_vals)) > 1
  t2d_has_variance <- length(unique(t2d_vals)) > 1
  
  # Skip if either group has no variance (all identical values)
  if (!control_has_variance | !t2d_has_variance) {
    cat("Skipping", var, "- no variance in one or both groups\n")
    next
  }
  
  # Check normality with error handling
  normality_control <- tryCatch({
    if (length(control_vals) >= 3 & length(control_vals) <= 5000) {
      shapiro.test(control_vals)$p.value > 0.05
    } else {
      TRUE
    }
  }, error = function(e) {
    # If Shapiro test fails, assume non-normal
    FALSE
  })
  
  normality_t2d <- tryCatch({
    if (length(t2d_vals) >= 3 & length(t2d_vals) <= 5000) {
      shapiro.test(t2d_vals)$p.value > 0.05
    } else {
      TRUE
    }
  }, error = function(e) {
    # If Shapiro test fails, assume non-normal
    FALSE
  })
  
  # Test for equal variances with error handling
  var_test_result <- tryCatch({
    var.test(control_vals, t2d_vals)
  }, error = function(e) {
    list(p.value = 0.5)  # Default to assuming equal variance if test fails
  })
  equal_var <- var_test_result$p.value > 0.05
  
  # Choose appropriate test
  if (normality_control & normality_t2d) {
    # T-test (parametric)
    test <- tryCatch({
      t.test(control_vals, t2d_vals, var.equal = equal_var)
    }, error = function(e) {
      # If t-test fails, use Mann-Whitney
      wilcox.test(control_vals, t2d_vals, exact = FALSE)
    })
    test_name <- ifelse(equal_var, "t-test", "Welch t-test")
    
  } else {
    # Mann-Whitney U (non-parametric)
    test <- tryCatch({
      wilcox.test(control_vals, t2d_vals, exact = FALSE)
    }, error = function(e) {
      # If Wilcox fails, try t-test as fallback
      t.test(control_vals, t2d_vals, var.equal = FALSE)
    })
    test_name <- "Mann-Whitney U"
  }
  
  # Skip if test completely failed
  if (is.null(test) | !("p.value" %in% names(test))) {
    cat("Skipping", var, "- statistical test failed\n")
    next
  }
  
  # Calculate means and SDs
  control_mean <- mean(control_vals, na.rm = TRUE)
  t2d_mean <- mean(t2d_vals, na.rm = TRUE)
  control_sd <- sd(control_vals, na.rm = TRUE)
  t2d_sd <- sd(t2d_vals, na.rm = TRUE)
  n_control <- length(control_vals)
  n_t2d <- length(t2d_vals)
  
  # Calculate effect size (Cohen's d) with error handling
  pooled_sd <- tryCatch({
    sqrt(((n_control - 1) * control_sd^2 + (n_t2d - 1) * t2d_sd^2) / 
           (n_control + n_t2d - 2))
  }, error = function(e) {
    # If pooled SD calculation fails, use average SD
    (control_sd + t2d_sd) / 2
  })
  
  # Avoid division by zero
  if (is.na(pooled_sd) | pooled_sd == 0) {
    cohens_d <- NA
  } else {
    cohens_d <- (t2d_mean - control_mean) / pooled_sd
  }
  
  # Store results
  results <- rbind(results, data.frame(
    variable = var,
    test_used = test_name,
    control_mean = control_mean,
    control_sd = control_sd,
    n_control = n_control,
    t2d_mean = t2d_mean,
    t2d_sd = t2d_sd,
    n_t2d = n_t2d,
    mean_difference = t2d_mean - control_mean,
    percent_change = ((t2d_mean - control_mean) / control_mean) * 100,
    p_value = test$p.value,
    cohens_d = cohens_d,
    statistic = ifelse("statistic" %in% names(test), test$statistic, NA),
    stringsAsFactors = FALSE
  ))
}

cat("\nTotal variables successfully tested:", nrow(results), "\n")

# ===== Step 4: Multiple testing correction =====

if (nrow(results) > 0) {
  results$p_adj_BH <- p.adjust(results$p_value, method = "BH")
  results$p_adj_bonferroni <- p.adjust(results$p_value, method = "bonferroni")
  
  # Add significance labels
  results$sig_raw <- ifelse(results$p_value < 0.001, "***",
                            ifelse(results$p_value < 0.01, "**",
                                   ifelse(results$p_value < 0.05, "*", "ns")))
  
  results$sig_BH <- ifelse(results$p_adj_BH < 0.001, "***",
                           ifelse(results$p_adj_BH < 0.01, "**",
                                  ifelse(results$p_adj_BH < 0.05, "*", "ns")))
  
  # Sort by p-value
  results <- results[order(results$p_value), ]
  
  # ===== Step 5: Save results table =====
  
  write.csv(results, "brain_statistics_all_variables.csv", row.names = FALSE)
  
  # Summary of significant findings
  sig_results <- results[results$p_adj_BH < 0.05, ]
  
  if (nrow(sig_results) > 0) {
    write.csv(sig_results, "brain_statistics_significant_only.csv", row.names = FALSE)
  }
  
  cat("\nTotal variables tested:", nrow(results))
  cat("\nSignificant (p < 0.05):", sum(results$p_value < 0.05))
  cat("\nSignificant after BH correction:", sum(results$p_adj_BH < 0.05))
  
  if (nrow(results) >= 10) {
    cat("\n\nTop 10 most significant:\n")
    print(results[1:10, c("variable", "p_value", "p_adj_BH", "mean_difference", "percent_change")])
  } else {
    cat("\n\nAll results:\n")
    print(results[, c("variable", "p_value", "p_adj_BH", "mean_difference", "percent_change")])
  }
  
  # ===== Step 6: Create plots for significant variables =====
  
  # Filter significant results
  sig_vars <- results$variable[results$p_adj_BH < 0.05]
  
  if (length(sig_vars) > 0) {
    # Create output directory
    dir.create("brain_plots", showWarnings = FALSE)
    
    # Plot each significant variable
    plot_list <- list()
    
    for (var in sig_vars) {
      
      # Get significance info
      var_result <- results[results$variable == var, ]
      
      # Create plot
      p <- ggplot(brain_data, aes(x = group2, y = .data[[var]], fill = group2)) +
        geom_boxplot(alpha = 0.7, outlier.shape = NA) +
        geom_jitter(width = 0.2, alpha = 0.5, size = 2) +
        stat_summary(fun = mean, geom = "point", shape = 23, size = 3, 
                     fill = "red", color = "black") +
        scale_fill_manual(values = c("Control" = "#00BFC4", 
                                     "Type 2 Diabetes" = "#F8766D")) +
        labs(
          title = gsub("_", " ", var),
          subtitle = sprintf("p = %.4f (adj = %.4f), d = %.2f",
                             var_result$p_value, 
                             var_result$p_adj_BH,
                             var_result$cohens_d),
          x = "",
          y = "Value",
          fill = "Group"
        ) +
        theme_classic(base_size = 12) +
        theme(
          plot.title = element_text(face = "bold", hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5, size = 10),
          legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 1)
        )
      
      # Add significance bracket
      y_max <- max(brain_data[[var]], na.rm = TRUE)
      y_min <- min(brain_data[[var]], na.rm = TRUE)
      y_range <- y_max - y_min
      
      p <- p + 
        annotate("segment", x = 1, xend = 2, 
                 y = y_max + y_range * 0.05, 
                 yend = y_max + y_range * 0.05) +
        annotate("text", x = 1.5, y = y_max + y_range * 0.08,
                 label = var_result$sig_BH, size = 6)
      
      # Save individual plot
      ggsave(filename = paste0("brain_plots/", gsub("[^A-Za-z0-9_]", "_", var), ".png"),
             plot = p, width = 6, height = 5, dpi = 300)
      
      plot_list[[var]] <- p
    }
    
    # ===== Step 7: Create combined plot of top findings =====
    
    n_plots <- min(length(sig_vars), 12)  # Show top 12
    
    combined_plot <- wrap_plots(plot_list[1:n_plots], ncol = 3)
    
    ggsave("brain_plots/combined_significant_variables.png",
           plot = combined_plot, width = 15, height = 12, dpi = 300)
  }
  
  # ===== Step 8: Summary statistics by category =====
  
  # Categorize variables
  results$category <- ifelse(grepl("^bilateral_", results$variable), "Bilateral (Combined)",
                             ifelse(grepl("^lh_|^Left-", results$variable), "Left Hemisphere",
                                    ifelse(grepl("^rh_|^Right-", results$variable), "Right Hemisphere",
                                           ifelse(grepl("thickness", results$variable, ignore.case = TRUE), "Thickness",
                                                  ifelse(grepl("volume|Vol", results$variable, ignore.case = TRUE), "Volume",
                                                         ifelse(grepl("area", results$variable, ignore.case = TRUE), "Surface Area", "Other"))))))
  
  # Summary by category
  category_summary <- results %>%
    group_by(category) %>%
    summarise(
      n_tested = n(),
      n_sig_raw = sum(p_value < 0.05),
      n_sig_corrected = sum(p_adj_BH < 0.05),
      mean_effect_size = mean(abs(cohens_d), na.rm = TRUE),
      .groups = 'drop'
    )
  
  print(category_summary)
  write.csv(category_summary, "brain_statistics_by_category.csv", row.names = FALSE)
  
  # ===== Step 9: Create volcano plot =====
  
  # Most standardized approach - uses Cohen's d
  volcano_plot <- ggplot(results, aes(x = cohens_d, y = -log10(p_value))) +
    geom_point(aes(color = sig_BH), alpha = 0.6, size = 2) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
    geom_vline(xintercept = c(-0.5, 0.5), linetype = "dotted", color = "gray50") +  # Medium effect size
    scale_color_manual(values = c("ns" = "gray", "*" = "orange", 
                                  "**" = "darkorange", "***" = "red"),
                       name = "Significance\n(BH corrected)") +
    labs(
      title = "Volcano Plot: Type 2 Diabetes vs Control",
      x = "Effect Size (Cohen's d)",
      y = "-log10(p-value)"
    ) +
    theme_classic(base_size = 12) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5))
  
  ggsave("brain_plots/volcano_plot.png", plot = volcano_plot, 
         width = 10, height = 7, dpi = 300)
  
  print(volcano_plot)
  
  # ===== Step 10: Effect size plot =====
  
  sig_for_effect <- results[results$p_adj_BH < 0.05, ]
  
  if (nrow(sig_for_effect) > 0) {
    n_show <- min(20, nrow(sig_for_effect))
    
    effect_plot_data <- sig_for_effect %>%
      arrange(desc(abs(cohens_d))) %>%
      head(n_show)
    
    effect_plot <- ggplot(effect_plot_data, 
                          aes(x = reorder(variable, cohens_d), y = cohens_d, fill = cohens_d > 0)) +
      geom_bar(stat = "identity") +
      coord_flip() +
      scale_fill_manual(values = c("TRUE" = "#F8766D", "FALSE" = "#00BFC4"),
                        labels = c("Decreased in T2D", "Increased in T2D"),
                        name = "") +
      labs(
        title = "Effect Sizes of Significant Variables",
        x = "",
        y = "Cohen's d"
      ) +
      theme_classic(base_size = 11) +
      theme(
        plot.title = element_text(face = "bold", hjust = 0.5),
        axis.text.y = element_text(size = 9)
      )
    
    ggsave("brain_plots/effect_sizes.png", plot = effect_plot,
           width = 10, height = 8, dpi = 300)
    
    print(effect_plot)
  }
  
  cat("\n\n===== Analysis Complete! =====\n")
  cat("Files saved:\n")
  cat("  - brain_statistics_all_variables.csv\n")
  if (nrow(sig_results) > 0) {
    cat("  - brain_statistics_significant_only.csv\n")
  }
  cat("  - brain_statistics_by_category.csv\n")
  if (length(sig_vars) > 0) {
    cat("  - brain_plots/ folder with individual plots\n")
    cat("  - brain_plots/combined_significant_variables.png\n")
  }
  cat("  - brain_plots/volcano_plot.png\n")
  if (nrow(sig_for_effect) > 0) {
    cat("  - brain_plots/effect_sizes.png\n")
  }
  
} else {
  cat("\nNo variables could be tested. Check your data.\n")
}





################# Compare brain biomarkers with imaging traits 
library(dplyr)
library(tidyr)
library(ggplot2)
library(pheatmap)
library(corrplot)
library(broom)

# ===== Define biomarkers and imaging traits =====

biomarkers <- c("ab40_avg_conc", "ab42_avg_conc", "tau_avg_conc", 
                "nfl_avg_conc", "gfap_avg_conc", "ptau_181_avg_conc", 
                "ptau_217_avg_conc")

# Get all imaging traits (numeric columns that aren't biomarkers or IDs)
imaging_traits <- names(brain_data)[sapply(brain_data, is.numeric)]
imaging_traits <- setdiff(imaging_traits, c("record_id", "file_id", biomarkers))

cat("Number of biomarkers:", length(biomarkers), "\n")
cat("Number of imaging traits:", length(imaging_traits), "\n")

# ===== PART 1: Overall Correlations (Heatmap) =====

# Calculate p-values for each correlation using a data frame approach
cor_results_detailed <- data.frame()

for (biomarker in biomarkers) {
  for (imaging_trait in imaging_traits) {
    
    test_data <- brain_data %>%
      select(all_of(c(biomarker, imaging_trait))) %>%
      na.omit()
    
    n_obs <- nrow(test_data)
    
    if (n_obs >= 3) {
      tryCatch({
        cor_val <- cor(test_data[[biomarker]], test_data[[imaging_trait]])
        test_result <- cor.test(test_data[[biomarker]], test_data[[imaging_trait]])
        p_val <- test_result$p.value
      }, error = function(e) {
        cor_val <- NA
        p_val <- NA
      })
    } else {
      cor_val <- NA
      p_val <- NA
    }
    
    cor_results_detailed <- rbind(cor_results_detailed, data.frame(
      biomarker = biomarker,
      imaging_trait = imaging_trait,
      correlation = cor_val,
      p_value = p_val,
      n_obs = n_obs,
      stringsAsFactors = FALSE
    ))
  }
}

# Adjust p-values for multiple testing
cor_results_detailed$p_adj <- p.adjust(cor_results_detailed$p_value, method = "BH")

# Add clean biomarker names
biomarker_labels <- c("ab40_avg_conc" = "Aβ40", 
                      "ab42_avg_conc" = "Aβ42", 
                      "tau_avg_conc" = "Tau", 
                      "nfl_avg_conc" = "NfL",
                      "gfap_avg_conc" = "GFAP", 
                      "ptau_181_avg_conc" = "pTau-181", 
                      "ptau_217_avg_conc" = "pTau-217")

cor_results_detailed$biomarker_name <- biomarker_labels[cor_results_detailed$biomarker]

# Create matrices for heatmap
cor_matrix_plot <- cor_results_detailed %>%
  select(biomarker, imaging_trait, correlation) %>%
  pivot_wider(names_from = imaging_trait, values_from = correlation) %>%
  column_to_rownames("biomarker") %>%
  as.matrix()

p_matrix_adj <- cor_results_detailed %>%
  select(biomarker, imaging_trait, p_adj) %>%
  pivot_wider(names_from = imaging_trait, values_from = p_adj) %>%
  column_to_rownames("biomarker") %>%
  as.matrix()

# Clean up row names
rownames(cor_matrix_plot) <- biomarker_labels[rownames(cor_matrix_plot)]
rownames(p_matrix_adj) <- biomarker_labels[rownames(p_matrix_adj)]

# Replace NAs with 0 for visualization (but keep track of them)
cor_matrix_plot_vis <- cor_matrix_plot
cor_matrix_plot_vis[is.na(cor_matrix_plot_vis)] <- 0

p_matrix_adj_vis <- p_matrix_adj
p_matrix_adj_vis[is.na(p_matrix_adj_vis)] <- 1

cat("Heatmap dimensions:", nrow(cor_matrix_plot), "x", ncol(cor_matrix_plot), "\n")

# Plot heatmap with pheatmap
png("brain_plots/biomarker_imaging_correlations.png", 
    width = 2400, height = 1000, res = 150)

# Only cluster if we have enough data points
cluster_rows_flag <- nrow(cor_matrix_plot_vis) > 2
cluster_cols_flag <- ncol(cor_matrix_plot_vis) > 2

# Create significance annotation matrix
sig_annotation <- matrix(
  ifelse(is.na(p_matrix_adj), "",
         ifelse(p_matrix_adj < 0.001, "***",
                ifelse(p_matrix_adj < 0.01, "**",
                       ifelse(p_matrix_adj < 0.05, "*", "")))),
  nrow = nrow(p_matrix_adj),
  dimnames = dimnames(p_matrix_adj)
)

pheatmap(cor_matrix_plot_vis,
         cluster_rows = cluster_rows_flag,
         cluster_cols = cluster_cols_flag,
         color = colorRampPalette(c("blue", "white", "red"))(100),
         breaks = seq(-1, 1, length.out = 101),
         display_numbers = sig_annotation,
         number_color = "black",
         fontsize_number = 6,
         main = "Biomarker-Brain Imaging Correlations (pairwise complete obs)",
         fontsize = 8,
         fontsize_row = 9,
         fontsize_col = 6,
         angle_col = 45,
         na_col = "grey90")

dev.off()

# Save results
write.csv(cor_results_detailed, "brain_biomarker_imaging_correlations_all.csv", 
          row.names = FALSE)

# Filter and sort significant correlations
top_correlations <- cor_results_detailed %>%
  filter(!is.na(correlation) & !is.na(p_adj) & p_adj < 0.05) %>%
  arrange(desc(abs(correlation)))

if (nrow(top_correlations) > 0) {
  write.csv(top_correlations, "brain_biomarker_imaging_correlations_significant.csv", 
            row.names = FALSE)
  
  cat("\nTop 10 significant biomarker-imaging correlations:\n")
  print(head(top_correlations[, c("biomarker_name", "imaging_trait", 
                                  "correlation", "p_adj", "n_obs")], 10))
} else {
  cat("\nNo significant correlations found after correction.\n")
}

# Print summary statistics
cat("\nCorrelation summary:\n")
cat("Total tests:", nrow(cor_results_detailed), "\n")
cat("Valid correlations:", sum(!is.na(cor_results_detailed$correlation)), "\n")
cat("Significant (p_adj < 0.05):", nrow(top_correlations), "\n")
cat("Mean sample size:", round(mean(cor_results_detailed$n_obs, na.rm = TRUE), 1), "\n")
cat("Median sample size:", median(cor_results_detailed$n_obs, na.rm = TRUE), "\n")

# ===== PART 2: Interaction Analysis =====
# Test: imaging ~ biomarker * group2

interaction_results <- data.frame()

for (biomarker in biomarkers) {
  for (imaging_trait in imaging_traits) {
    
    # Prepare data
    model_data <- brain_data %>%
      select(group2, all_of(biomarker), all_of(imaging_trait)) %>%
      filter(!is.na(group2) & 
               !is.na(.data[[biomarker]]) & 
               !is.na(.data[[imaging_trait]]))
    
    # Skip if insufficient data
    if (nrow(model_data) < 10) next
    
    # Check if there's variance in predictor
    if (length(unique(model_data[[biomarker]])) < 2) next
    
    # Fit interaction model
    tryCatch({
      # Fit model - use column indexing to avoid name issues
      model <- lm(model_data[[imaging_trait]] ~ 
                    model_data[[biomarker]] * model_data$group2)
      
      # Extract results
      model_summary <- summary(model)
      coefs <- coef(model_summary)
      
      # Get interaction term (last row)
      if (nrow(coefs) >= 4) {
        interaction_coef <- coefs[4, "Estimate"]
        interaction_se <- coefs[4, "Std. Error"]
        interaction_p <- coefs[4, "Pr(>|t|)"]
      } else {
        interaction_coef <- NA
        interaction_se <- NA
        interaction_p <- NA
      }
      
      # Get biomarker main effect (second row)
      biomarker_coef <- coefs[2, "Estimate"]
      biomarker_p <- coefs[2, "Pr(>|t|)"]
      
      # Store results
      interaction_results <- rbind(interaction_results, data.frame(
        biomarker = biomarker,
        imaging_trait = imaging_trait,
        n_obs = nrow(model_data),
        r_squared = model_summary$r.squared,
        adj_r_squared = model_summary$adj.r.squared,
        biomarker_coef = biomarker_coef,
        biomarker_p = biomarker_p,
        interaction_coef = interaction_coef,
        interaction_se = interaction_se,
        interaction_p = interaction_p,
        stringsAsFactors = FALSE
      ))
      
    }, error = function(e) {
      # Skip if model fails
    })
  }
}

cat("\nTotal interaction models tested:", nrow(interaction_results), "\n")

# Multiple testing correction
if (nrow(interaction_results) > 0) {
  interaction_results$biomarker_p_adj <- p.adjust(interaction_results$biomarker_p, 
                                                  method = "BH")
  interaction_results$interaction_p_adj <- p.adjust(interaction_results$interaction_p, 
                                                    method = "BH")
  
  # Filter significant interactions
  sig_interactions <- interaction_results %>%
    filter(interaction_p_adj < 0.05) %>%
    arrange(interaction_p_adj)
  
  cat("Significant interactions found:", nrow(sig_interactions), "\n")
  
  # Save results
  write.csv(interaction_results, "biomarker_imaging_interactions_all.csv", 
            row.names = FALSE)
  
  if (nrow(sig_interactions) > 0) {
    write.csv(sig_interactions, "biomarker_imaging_interactions_significant.csv", 
              row.names = FALSE)
  }
  
  # ===== PART 3: Visualize Top Interactions =====
  
  if (nrow(sig_interactions) > 0) {
    
    # Plot top 6 interactions
    n_plot <- min(6, nrow(sig_interactions))
    
    plot_list <- list()
    
    for (i in 1:n_plot) {
      biomarker <- sig_interactions$biomarker[i]
      imaging_trait <- sig_interactions$imaging_trait[i]
      
      # Clean names for display
      biomarker_clean <- biomarker_labels[biomarker]
      imaging_clean <- gsub("_", " ", imaging_trait)
      
      p <- ggplot(brain_data, 
                  aes(x = .data[[biomarker]], 
                      y = .data[[imaging_trait]], 
                      color = group2)) +
        geom_point(alpha = 0.6, size = 2.5) +
        geom_smooth(method = "lm", se = TRUE, linewidth = 1.2) +
        scale_color_manual(values = c("Control" = "#00BFC4", 
                                      "Type 2 Diabetes" = "#F8766D")) +
        labs(
          title = paste(biomarker_clean, "×", imaging_clean),
          subtitle = sprintf("Interaction p(adj) = %.4f, R² = %.3f", 
                             sig_interactions$interaction_p_adj[i],
                             sig_interactions$r_squared[i]),
          x = paste(biomarker_clean, "(pg/mL)"),
          y = imaging_clean,
          color = "Group"
        ) +
        theme_classic(base_size = 11) +
        theme(
          plot.title = element_text(face = "bold", hjust = 0.5, size = 12),
          plot.subtitle = element_text(hjust = 0.5, size = 9),
          legend.position = "bottom"
        )
      
      plot_list[[i]] <- p
      
      # Save individual plot
      ggsave(paste0("brain_plots/interaction_", i, "_", 
                    gsub("[^A-Za-z0-9]", "_", biomarker), "_", 
                    gsub("[^A-Za-z0-9]", "_", imaging_trait), ".png"),
             plot = p, width = 6, height = 5, dpi = 300)
    }
    
    # Combined plot
    combined_interactions <- wrap_plots(plot_list, ncol = 2)
    
    ggsave("brain_plots/top_biomarker_imaging_interactions.png",
           plot = combined_interactions, width = 12, height = 9, dpi = 300)
  }
  
  # ===== PART 4: Heatmap of Interaction P-values =====
  
  interaction_p_long <- interaction_results %>%
    select(biomarker, imaging_trait, interaction_p_adj) %>%
    filter(!is.na(interaction_p_adj)) %>%
    mutate(neg_log_p = -log10(interaction_p_adj + 1e-10))
  
  interaction_p_matrix <- interaction_p_long %>%
    select(biomarker, imaging_trait, neg_log_p) %>%
    pivot_wider(names_from = imaging_trait, values_from = neg_log_p) %>%
    column_to_rownames("biomarker") %>%
    as.matrix()
  
  # Replace NAs with 0
  interaction_p_matrix[is.na(interaction_p_matrix)] <- 0
  
  # Clean row names
  biomarker_names_present <- rownames(interaction_p_matrix)
  rownames(interaction_p_matrix) <- biomarker_labels[biomarker_names_present]
  
  png("brain_plots/interaction_pvalues_heatmap.png", 
      width = 2400, height = 1000, res = 150)
  
  pheatmap(interaction_p_matrix,
           cluster_rows = nrow(interaction_p_matrix) > 2,
           cluster_cols = ncol(interaction_p_matrix) > 2,
           color = colorRampPalette(c("white", "yellow", "orange", "red", "darkred"))(100),
           main = "Biomarker × Group Interactions (-log10 adjusted p-value)",
           fontsize = 8,
           fontsize_row = 9,
           fontsize_col = 6,
           angle_col = 45,
           breaks = seq(0, max(interaction_p_matrix, na.rm = TRUE), 
                        length.out = 101))
  
  dev.off()
  
  # ===== PART 5: Summary Statistics =====
  
  # Create summary table
  biomarker_summary <- interaction_results %>%
    group_by(biomarker) %>%
    summarise(
      n_imaging_traits = n(),
      n_sig_main_effect = sum(biomarker_p_adj < 0.05, na.rm = TRUE),
      n_sig_interaction = sum(interaction_p_adj < 0.05, na.rm = TRUE),
      mean_r_squared = mean(r_squared, na.rm = TRUE),
      .groups = 'drop'
    ) %>%
    mutate(
      biomarker_name = biomarker_labels[biomarker]
    ) %>%
    select(biomarker_name, everything(), -biomarker)
  
  print("\n=== Biomarker Summary ===")
  print(biomarker_summary)
  
  write.csv(biomarker_summary, "biomarker_imaging_summary.csv", row.names = FALSE)
  
  cat("\n\n===== Analysis Complete! =====\n")
  cat("Files created:\n")
  cat("  - brain_plots/biomarker_imaging_correlations.png\n")
  cat("  - brain_biomarker_imaging_correlations_all.csv\n")
  if (nrow(top_correlations) > 0) {
    cat("  - brain_biomarker_imaging_correlations_significant.csv\n")
  }
  cat("  - biomarker_imaging_interactions_all.csv\n")
  if (nrow(sig_interactions) > 0) {
    cat("  - biomarker_imaging_interactions_significant.csv\n")
    cat("  - brain_plots/top_biomarker_imaging_interactions.png\n")
    cat("  - Individual interaction plots\n")
  }
  cat("  - brain_plots/interaction_pvalues_heatmap.png\n")
  cat("  - biomarker_imaging_summary.csv\n")
} else {
  cat("\nNo interaction models could be fitted.\n")
}



biomarker_imaging_interactions <- data.table::fread("biomarker_imaging_interactions_all.csv")





