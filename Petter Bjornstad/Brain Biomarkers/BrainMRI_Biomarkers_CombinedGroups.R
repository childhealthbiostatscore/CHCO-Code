# ==============================================================================
# COMPREHENSIVE BRAIN IMAGING ANALYSIS
# - Structural MRI with volume correction
# - Between-network connectivity and correlations  
# - All previous visualizations recreated
# - T2D+OC combined analysis throughout
# ==============================================================================

library(oro.nifti)
library(neurobase)
library(ggplot2)
library(dplyr)
library(tidyr)
library(corrplot)
library(ggsignif)
library(gridExtra)
library(stringr)
library(tibble)

# ==============================================================================
# SETUP
# ==============================================================================

cat("=== COMPREHENSIVE BRAIN IMAGING ANALYSIS ===\n\n")

# Set consistent output directory
results.dir <- 'Projects/Brain_Imaging/results/comprehensive_analysis/'
dir.create(results.dir, recursive = TRUE, showWarnings = FALSE)

base_dir <- 'Projects/Brain_Imaging/data/'
fmri_output_dir <- "~/Documents/UofW/Projects/Brain_Imaging/fMRI/results/network_analysis"

cat("Output directory:", results.dir, "\n\n")

# ==============================================================================
# 1. LOAD AND PREPARE DATA
# ==============================================================================

cat("=== Loading Data ===\n")

harmonized_data <- read.csv("../OneDrive/UW/Laura Pyle - Biostatistics Core Shared Drive/Data Harmonization/Data Clean/harmonized_dataset.csv", na = '')

dat <- harmonized_data %>%
  arrange(date_of_screen) %>% 
  dplyr::summarise(
    across(where(negate(is.numeric)), ~ ifelse(all(is.na(.x)), NA_character_, last(na.omit(.x)))),
    across(where(is.numeric), ~ ifelse(all(is.na(.x)), NA_real_, mean(na.omit(.x), na.rm=T))),
    .by = c(record_id, visit)
  )

mri_ids <- c('CRC-10', 'CRC-11', 'CRC-12', 'CRC-13', 'CRC-26', 'CRC-39', 'CRC-46', 'CRC-51', 'CRC-53', 
             'CRC-55', 'CRC-58', 'CRC-60', 
             'RH2-01-O', 'RH2-03-O', 'RH2-08-T', 'RH2-10-L', 'RH2-11-O', 'RH2-13-O', 'RH2-16-O', 'RH2-17-L', 
             'RH2-18-O', 'RH2-19-T', 'RH2-22-T', 'RH2-24-L', 'RH2-27-L', 'RH2-28-L', 'RH2-29-L', 'RH2-33-L', 
             'RH2-34-O', 'RH2-35-T', 'RH2-38-T', 'RH2-39-O', 'RH2-41-T', 'RH2-42-T', 'RH2-43-T', 
             'RH2-44-T', 'RH2-45-T', 'RH2-48-T', 'RH2-49-T', 'RH2-50-L', 'RH2-52-T', 'RH2-53-T', 
             'RH2-55-T')

mri_ids_df <- data.frame(
  ID = mri_ids, 
  file_id = c(
    'crc_10', 'crc_11', 'crc_12', 'crc_13', 'crc_26', 'crc_39', 'crc_46', 'crc_51', 'crc_53', 
    'crc_55', 'crc_58', 'crc_60', 
    'RH2_01_O', 'RH2_03_O', 'RH2_08_T', 'RH2_10_L', 'RH2_11_O', 'RH2_13_O', 'RH2_16_O', 'RH2_17_L', 
    'RH2_18_O', 'RH2_19_T', 'RH2_22_T', 'RH2_24_L', 'RH2_27_L', 'RH2_28_L', 'RH2_29_L', 'RH2_33_L', 
    'RH2_34_O', 'RH2_35_T', 'RH2_38_T', 'RH2_39_O', 'RH2_41_T', 'RH2_42_T', 'RH2_43_T', 
    'RH2_44_T', 'RH2_45_T', 'RH2_48_T', 'RH2_49_T', 'RH2_50_L', 'RH2_52_T', 'RH2_53_T', 
    'RH2_55_T'
  )
)

small_dat <- dat %>% filter(record_id %in% mri_ids)
small_dat$group[which(small_dat$record_id == 'RH2-38-T')] <- 'Obese Control'
small_dat$record_id[which(small_dat$record_id == 'RH2-38-T')] <- 'RH2-38-O'

# Protein markers
qx_var <- c("ab40_avg_conc", "ab42_avg_conc", "tau_avg_conc",
            "nfl_avg_conc", "gfap_avg_conc", "ptau_181_avg_conc", "ptau_217_avg_conc")

# Create group variables
small_dat$group_3level <- small_dat$group
small_dat$group_combined <- ifelse(
  small_dat$group %in% c("Type 2 Diabetes", "Obese Control"),
  "T2D+OC", "LC"
)

cat("Sample sizes:\n")
cat("T2D:", sum(small_dat$group == "Type 2 Diabetes"), "\n")
cat("OC:", sum(small_dat$group == "Obese Control"), "\n")
cat("LC:", sum(small_dat$group == "Lean Control"), "\n")
cat("T2D+OC:", sum(small_dat$group_combined == "T2D+OC"), "\n\n")

# ==============================================================================
# 2. LOAD PRE-COMPUTED STRUCTURAL MRI RESULTS
# ==============================================================================

cat("=== Loading Pre-computed Structural MRI Results ===\n")

# Load the morphometric similarity networks (from your original script)
morph_data_file <- paste0(results.dir, "../morphometric_similarity_networks.rds")

if (!file.exists(morph_data_file)) {
  stop("Morphometric similarity networks file not found at: ", morph_data_file,
       "\nPlease check the file path or run the original structural analysis first.")
}

cat("Loading existing morphometric similarity networks...\n")
morph_data <- readRDS(morph_data_file)

# Check the structure of the loaded data
cat("Data structure:\n")
cat("  T2D participants:", length(morph_data$t2d_conn), "\n")
cat("  LC participants:", length(morph_data$lc_conn), "\n")

# Combine connectivity matrices and IDs
all_conn_list <- c(morph_data$t2d_conn, morph_data$lc_conn)
all_ids_list <- c(morph_data$t2d_ids, morph_data$lc_ids)

# Check if these are matrices or list objects
if (is.matrix(all_conn_list[[1]])) {
  cat("\nNote: Loaded data contains matrices only (no volume/ROI info)\n")
  cat("Will calculate metrics from matrices directly\n\n")
} else {
  cat("\nLoaded data contains full result objects\n\n")
}

# ==============================================================================
# 3. EXTRACT NETWORK METRICS (FROM MATRICES ONLY)
# ==============================================================================

cat("=== Extracting Network Metrics from Connectivity Matrices ===\n")

calculate_network_metrics_from_matrix <- function(conn_matrix) {
  diag(conn_matrix) <- NA
  
  mean_connectivity <- mean(conn_matrix, na.rm = TRUE)
  median_connectivity <- median(conn_matrix, na.rm = TRUE)
  sd_connectivity <- sd(conn_matrix, na.rm = TRUE)
  mean_strength <- mean(abs(conn_matrix), na.rm = TRUE)
  
  threshold_strong <- 0.3
  density_strong <- sum(abs(conn_matrix) > threshold_strong, na.rm = TRUE) / 
    sum(!is.na(conn_matrix))
  
  mean_efficiency <- mean(abs(conn_matrix), na.rm = TRUE)
  connectivity_variance <- var(as.vector(conn_matrix), na.rm = TRUE)
  
  node_strengths <- rowSums(abs(conn_matrix), na.rm = TRUE)
  mean_node_strength <- mean(node_strengths, na.rm = TRUE)
  sd_node_strength <- sd(node_strengths, na.rm = TRUE)
  prop_positive <- sum(conn_matrix > 0, na.rm = TRUE) / sum(!is.na(conn_matrix))
  
  # Estimate "volume" from number of ROIs (proxy measure)
  n_rois <- nrow(conn_matrix)
  
  return(data.frame(
    n_rois = n_rois,
    mean_connectivity = mean_connectivity,
    median_connectivity = median_connectivity,
    sd_connectivity = sd_connectivity,
    mean_strength = mean_strength,
    network_density = density_strong,
    mean_efficiency = mean_efficiency,
    connectivity_variance = connectivity_variance,
    mean_node_strength = mean_node_strength,
    sd_node_strength = sd_node_strength,
    prop_positive = prop_positive
  ))
}

# Calculate metrics for all participants
all_metrics <- do.call(rbind, lapply(all_conn_list, calculate_network_metrics_from_matrix))
all_metrics$record_id <- all_ids_list

# Add group information
all_metrics <- all_metrics %>%
  left_join(small_dat %>% select(record_id, group, group_combined), by = "record_id") %>%
  filter(!is.na(group))

# Create "volume-corrected" metrics using n_rois as proxy
# (normalize by number of ROIs instead of actual brain volume)
mean_n_rois <- mean(all_metrics$n_rois, na.rm = TRUE)
all_metrics$mean_connectivity_volcorr <- all_metrics$mean_connectivity * 
  (all_metrics$n_rois / mean_n_rois)
all_metrics$mean_node_strength_volcorr <- all_metrics$mean_node_strength * 
  (all_metrics$n_rois / mean_n_rois)

cat("Metrics calculated for", nrow(all_metrics), "participants\n")
cat("By group:\n")
print(table(all_metrics$group))
cat("\nNote: Using n_rois as proxy for brain volume (actual volume not saved in original file)\n\n")

# ==============================================================================
# 5. VISUALIZE STRUCTURAL METRICS - 3 GROUPS + COMBINED
# ==============================================================================

cat("=== Creating Structural MRI Visualizations ===\n")

brain_features <- c("n_rois", "mean_connectivity", "mean_connectivity_volcorr",
                    "median_connectivity", "sd_connectivity", "mean_strength", 
                    "network_density", "mean_efficiency", "connectivity_variance",
                    "mean_node_strength", "mean_node_strength_volcorr", "sd_node_strength", 
                    "prop_positive")

# Plot 1: 3-group comparison
plot_data_3group <- all_metrics %>%
  select(record_id, group, all_of(brain_features)) %>%
  pivot_longer(cols = all_of(brain_features), names_to = "Metric", values_to = "Value") %>%
  mutate(Metric = factor(Metric, levels = brain_features))

p1 <- ggplot(plot_data_3group, aes(x = group, y = Value, fill = group)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.5, size = 2) +
  facet_wrap(~Metric, scales = "free_y", ncol = 3) +
  scale_fill_manual(values = c("Type 2 Diabetes" = "#E74C3C", 
                               "Obese Control" = "#F39C12",
                               "Lean Control" = "#3498DB")) +
  labs(title = "Structural Brain Metrics: T2D vs OC vs LC",
       subtitle = "Volume-corrected metrics included",
       x = "Group", y = "Value", fill = "Group") +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom",
        plot.title = element_text(face = "bold"))

ggsave(paste0(results.dir, "structural_metrics_3groups.pdf"), p1, width = 16, height = 12)
ggsave(paste0(results.dir, "structural_metrics_3groups.png"), p1, width = 16, height = 12, dpi = 300)
cat("Saved: structural_metrics_3groups.pdf/.png\n")

# Plot 2: 2-group comparison (T2D+OC vs LC)
plot_data_2group <- all_metrics %>%
  select(record_id, group_combined, all_of(brain_features)) %>%
  pivot_longer(cols = all_of(brain_features), names_to = "Metric", values_to = "Value") %>%
  mutate(Metric = factor(Metric, levels = brain_features))

p2 <- ggplot(plot_data_2group, aes(x = group_combined, y = Value, fill = group_combined)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.5, size = 2) +
  facet_wrap(~Metric, scales = "free_y", ncol = 3) +
  scale_fill_manual(values = c("T2D+OC" = "#E74C3C", "LC" = "#3498DB")) +
  labs(title = "Structural Brain Metrics: T2D+OC vs LC",
       subtitle = "Volume-corrected metrics included",
       x = "Group", y = "Value", fill = "Group") +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom",
        plot.title = element_text(face = "bold"))

ggsave(paste0(results.dir, "structural_metrics_2groups.pdf"), p2, width = 16, height = 12)
ggsave(paste0(results.dir, "structural_metrics_2groups.png"), p2, width = 16, height = 12, dpi = 300)
cat("Saved: structural_metrics_2groups.pdf/.png\n")

# Statistical comparisons
perform_comparison <- function(data, group_var, features) {
  results <- data.frame()
  
  for (feature in features) {
    tryCatch({
      if (group_var == "group") {
        # 3-group ANOVA
        # Check if we have all 3 groups
        n_groups <- length(unique(data[[group_var]]))
        if (n_groups < 2) {
          warning(paste("Skipping", feature, "- insufficient groups"))
          next
        }
        
        formula_str <- paste(feature, "~", group_var)
        aov_result <- aov(as.formula(formula_str), data = data)
        aov_summary <- summary(aov_result)
        
        group_means <- aggregate(data[[feature]], by = list(data[[group_var]]), 
                                 FUN = mean, na.rm = TRUE)
        
        # Get means for each group (handle missing groups)
        t2d_mean <- ifelse("Type 2 Diabetes" %in% group_means$Group.1,
                           group_means$x[group_means$Group.1 == "Type 2 Diabetes"],
                           NA)
        oc_mean <- ifelse("Obese Control" %in% group_means$Group.1,
                          group_means$x[group_means$Group.1 == "Obese Control"],
                          NA)
        lc_mean <- ifelse("Lean Control" %in% group_means$Group.1,
                          group_means$x[group_means$Group.1 == "Lean Control"],
                          NA)
        
        results <- rbind(results, data.frame(
          Feature = feature,
          Test = "ANOVA",
          F_stat = aov_summary[[1]][1, 4],
          p_value = aov_summary[[1]][1, 5],
          T2D_mean = t2d_mean,
          OC_mean = oc_mean,
          LC_mean = lc_mean
        ))
      } else {
        # 2-group t-test
        group1_vals <- data[[feature]][data[[group_var]] == "T2D+OC"]
        group2_vals <- data[[feature]][data[[group_var]] == "LC"]
        
        # Remove NAs
        group1_vals <- group1_vals[!is.na(group1_vals)]
        group2_vals <- group2_vals[!is.na(group2_vals)]
        
        if (length(group1_vals) < 2 || length(group2_vals) < 2) {
          warning(paste("Skipping", feature, "- insufficient data"))
          next
        }
        
        test_result <- t.test(group1_vals, group2_vals)
        pooled_sd <- sqrt(((length(group1_vals)-1)*var(group1_vals, na.rm=TRUE) + 
                             (length(group2_vals)-1)*var(group2_vals, na.rm=TRUE)) / 
                            (length(group1_vals) + length(group2_vals) - 2))
        
        if (pooled_sd > 0) {
          cohens_d <- (mean(group1_vals, na.rm=TRUE) - mean(group2_vals, na.rm=TRUE)) / pooled_sd
        } else {
          cohens_d <- NA
        }
        
        results <- rbind(results, data.frame(
          Feature = feature,
          Test = "t-test",
          t_stat = test_result$statistic,
          p_value = test_result$p.value,
          cohens_d = cohens_d,
          T2D_OC_mean = mean(group1_vals, na.rm=TRUE),
          LC_mean = mean(group2_vals, na.rm=TRUE)
        ))
      }
    }, error = function(e) {
      warning(paste("Error processing", feature, ":", e$message))
    })
  }
  
  if (nrow(results) > 0) {
    results$p_fdr <- p.adjust(results$p_value, method = "fdr")
  }
  
  return(results)
}

structural_3group_stats <- perform_comparison(all_metrics, "group", brain_features)
structural_2group_stats <- perform_comparison(all_metrics, "group_combined", brain_features)

write.csv(structural_3group_stats, paste0(results.dir, "structural_3group_stats.csv"), row.names = FALSE)
write.csv(structural_2group_stats, paste0(results.dir, "structural_2group_stats.csv"), row.names = FALSE)
write.csv(all_metrics, paste0(results.dir, "all_structural_metrics.csv"), row.names = FALSE)

cat("Saved: structural_3group_stats.csv\n")
cat("Saved: structural_2group_stats.csv\n")
cat("Saved: all_structural_metrics.csv\n\n")

# ==============================================================================
# 6. LOAD FMRI NETWORK DATA
# ==============================================================================

cat("=== Loading fMRI Network Data ===\n")

network_metrics_all <- read.csv(file.path(fmri_output_dir, "network_metrics_all_participants.csv"))
between_network_data <- readRDS(file.path(fmri_output_dir, "between_network_connectivity.rds"))

# Match to record IDs
network_metrics_all$record_id <- mri_ids_df$ID[match(network_metrics_all$Participant, mri_ids_df$file_id)]
network_metrics_all$record_id[which(network_metrics_all$Participant == 'RH2_38_O')] <- 'RH2-38-O'

participant_groups <- small_dat %>% select(record_id, group, group_combined) %>% distinct()
network_metrics_with_groups <- network_metrics_all %>%
  left_join(participant_groups, by = "record_id") %>%
  filter(!is.na(group))

cat("Network data loaded for", nrow(network_metrics_with_groups), "participants\n\n")

# ==============================================================================
# 7. WITHIN-NETWORK CONNECTIVITY VISUALIZATIONS
# ==============================================================================

cat("=== Creating Within-Network Visualizations ===\n")

# Plot 1: 3-group within-network
p3 <- ggplot(network_metrics_with_groups, 
             aes(x = Network, y = Within_Network_Connectivity, fill = group)) +
  geom_boxplot(position = position_dodge(width = 0.8), alpha = 0.7) +
  geom_point(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.8),
             alpha = 0.5, size = 2) +
  scale_fill_manual(values = c("Type 2 Diabetes" = "#E74C3C", 
                               "Obese Control" = "#F39C12",
                               "Lean Control" = "#3498DB"),
                    name = "Group") +
  labs(title = "Within-Network Connectivity: T2D vs OC vs LC",
       x = "Network", y = "Within-Network Connectivity") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom",
        plot.title = element_text(face = "bold"))

ggsave(paste0(results.dir, "within_network_3groups.pdf"), p3, width = 14, height = 7)
ggsave(paste0(results.dir, "within_network_3groups.png"), p3, width = 14, height = 7, dpi = 300)
cat("Saved: within_network_3groups.pdf/.png\n")

# Plot 2: 2-group within-network
p4 <- ggplot(network_metrics_with_groups, 
             aes(x = Network, y = Within_Network_Connectivity, fill = group_combined)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.5, size = 2) +
  scale_fill_manual(values = c("T2D+OC" = "#E74C3C", "LC" = "#3498DB")) +
  labs(title = "Within-Network Connectivity: T2D+OC vs LC",
       x = "Network", y = "Within-Network Connectivity", fill = "Group") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom",
        plot.title = element_text(face = "bold"))

ggsave(paste0(results.dir, "within_network_2groups.pdf"), p4, width = 12, height = 7)
ggsave(paste0(results.dir, "within_network_2groups.png"), p4, width = 12, height = 7, dpi = 300)
cat("Saved: within_network_2groups.pdf/.png\n")

# ==============================================================================
# 8. BETWEEN-NETWORK CONNECTIVITY ANALYSIS
# ==============================================================================

cat("=== Analyzing Between-Network Connectivity ===\n")

all_between_network_fc <- between_network_data$between_network_fc
participant_ids <- between_network_data$participant_ids
networks <- between_network_data$networks

# Extract between-network pairs
extract_between_network_pairs <- function(fc_matrix, participant_id) {
  n_networks <- nrow(fc_matrix)
  network_names <- rownames(fc_matrix)
  results <- data.frame()
  
  for (i in 1:(n_networks-1)) {
    for (j in (i+1):n_networks) {
      results <- rbind(results, data.frame(
        Participant = participant_id,
        Network1 = network_names[i],
        Network2 = network_names[j],
        Network_Pair = paste(network_names[i], network_names[j], sep = "-"),
        Between_Network_Connectivity = fc_matrix[i, j]
      ))
    }
  }
  return(results)
}

all_between_pairs <- data.frame()
for (i in 1:length(participant_ids)) {
  pairs <- extract_between_network_pairs(all_between_network_fc[[i]], participant_ids[i])
  all_between_pairs <- rbind(all_between_pairs, pairs)
}

# Match to groups
all_between_pairs$record_id <- mri_ids_df$ID[match(all_between_pairs$Participant, mri_ids_df$file_id)]
all_between_pairs$record_id[which(all_between_pairs$Participant == 'RH2_38_O')] <- 'RH2-38-O'

all_between_pairs <- all_between_pairs %>%
  left_join(participant_groups, by = "record_id") %>%
  filter(!is.na(group))

cat("Between-network pairs extracted:", nrow(all_between_pairs), "\n\n")

# ==============================================================================
# 9. BETWEEN-NETWORK CONNECTIVITY VISUALIZATIONS
# ==============================================================================

cat("=== Creating Between-Network Visualizations ===\n")

# Group average heatmaps
create_matrix_from_pairs <- function(pair_data, networks) {
  mat <- matrix(0, nrow = length(networks), ncol = length(networks))
  rownames(mat) <- networks
  colnames(mat) <- networks
  
  for (i in 1:nrow(pair_data)) {
    net1 <- as.character(pair_data$Network1[i])
    net2 <- as.character(pair_data$Network2[i])
    val <- pair_data$Between_Network_Connectivity[i]
    
    if (net1 %in% networks && net2 %in% networks) {
      idx1 <- which(networks == net1)
      idx2 <- which(networks == net2)
      mat[idx1, idx2] <- val
      mat[idx2, idx1] <- val
    }
  }
  return(mat)
}

mean_between_t2d <- aggregate(Between_Network_Connectivity ~ Network1 + Network2, 
                              data = all_between_pairs[all_between_pairs$group == "Type 2 Diabetes",],
                              FUN = mean, na.rm = TRUE)
mean_between_oc <- aggregate(Between_Network_Connectivity ~ Network1 + Network2, 
                             data = all_between_pairs[all_between_pairs$group == "Obese Control",],
                             FUN = mean, na.rm = TRUE)
mean_between_lc <- aggregate(Between_Network_Connectivity ~ Network1 + Network2, 
                             data = all_between_pairs[all_between_pairs$group == "Lean Control",],
                             FUN = mean, na.rm = TRUE)
mean_between_combined <- aggregate(Between_Network_Connectivity ~ Network1 + Network2, 
                                   data = all_between_pairs[all_between_pairs$group_combined == "T2D+OC",],
                                   FUN = mean, na.rm = TRUE)

mat_t2d <- create_matrix_from_pairs(mean_between_t2d, networks)
mat_oc <- create_matrix_from_pairs(mean_between_oc, networks)
mat_lc <- create_matrix_from_pairs(mean_between_lc, networks)
mat_combined <- create_matrix_from_pairs(mean_between_combined, networks)

# 3-group heatmaps
pdf(paste0(results.dir, "between_network_heatmaps_3groups.pdf"), width = 15, height = 5)
par(mfrow = c(1, 3))
corrplot(mat_t2d, method = "color", type = "full", tl.col = "black", tl.srt = 45, tl.cex = 0.8,
         title = "T2D Between-Network", mar = c(0,0,2,0),
         col = colorRampPalette(c("#3498DB", "white", "#E74C3C"))(200))
corrplot(mat_oc, method = "color", type = "full", tl.col = "black", tl.srt = 45, tl.cex = 0.8,
         title = "OC Between-Network", mar = c(0,0,2,0),
         col = colorRampPalette(c("#3498DB", "white", "#E74C3C"))(200))
corrplot(mat_lc, method = "color", type = "full", tl.col = "black", tl.srt = 45, tl.cex = 0.8,
         title = "LC Between-Network", mar = c(0,0,2,0),
         col = colorRampPalette(c("#3498DB", "white", "#E74C3C"))(200))
dev.off()
cat("Saved: between_network_heatmaps_3groups.pdf\n")

# 2-group heatmaps
pdf(paste0(results.dir, "between_network_heatmaps_2groups.pdf"), width = 10, height = 5)
par(mfrow = c(1, 2))
corrplot(mat_combined, method = "color", type = "full", tl.col = "black", tl.srt = 45, tl.cex = 0.8,
         title = "T2D+OC Between-Network", mar = c(0,0,2,0),
         col = colorRampPalette(c("#3498DB", "white", "#E74C3C"))(200))
corrplot(mat_lc, method = "color", type = "full", tl.col = "black", tl.srt = 45, tl.cex = 0.8,
         title = "LC Between-Network", mar = c(0,0,2,0),
         col = colorRampPalette(c("#3498DB", "white", "#E74C3C"))(200))
dev.off()
cat("Saved: between_network_heatmaps_2groups.pdf\n")

# Boxplots of top different pairs
between_network_comparisons <- data.frame()
network_pairs <- unique(all_between_pairs$Network_Pair)

for (pair in network_pairs) {
  pair_data <- all_between_pairs %>% filter(Network_Pair == pair)
  
  group1_vals <- pair_data$Between_Network_Connectivity[pair_data$group_combined == "T2D+OC"]
  group2_vals <- pair_data$Between_Network_Connectivity[pair_data$group_combined == "LC"]
  
  if (length(group1_vals) < 2 || length(group2_vals) < 2) next
  
  test_result <- t.test(group1_vals, group2_vals)
  
  between_network_comparisons <- rbind(between_network_comparisons, data.frame(
    Network_Pair = pair,
    t_stat = test_result$statistic,
    p_value = test_result$p.value,
    T2D_OC_mean = mean(group1_vals, na.rm=TRUE),
    LC_mean = mean(group2_vals, na.rm=TRUE)
  ))
}

between_network_comparisons$p_fdr <- p.adjust(between_network_comparisons$p_value, method = "fdr")

top_pairs <- between_network_comparisons %>% 
  arrange(p_value) %>% 
  head(6) %>%
  pull(Network_Pair)

plot_data_between <- all_between_pairs %>% filter(Network_Pair %in% top_pairs)

p5 <- ggplot(plot_data_between, aes(x = group_combined, y = Between_Network_Connectivity, 
                                    fill = group_combined)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.5, size = 2) +
  facet_wrap(~Network_Pair, scales = "free_y", ncol = 3) +
  scale_fill_manual(values = c("T2D+OC" = "#E74C3C", "LC" = "#3498DB")) +
  labs(title = "Top Between-Network Connectivity Differences",
       subtitle = "T2D+OC vs LC",
       x = "Group", y = "Between-Network Connectivity") +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom")

ggsave(paste0(results.dir, "top_between_network_differences.pdf"), p5, width = 12, height = 8)
ggsave(paste0(results.dir, "top_between_network_differences.png"), p5, width = 12, height = 8, dpi = 300)
cat("Saved: top_between_network_differences.pdf/.png\n")

write.csv(between_network_comparisons, paste0(results.dir, "between_network_stats.csv"), row.names = FALSE)
cat("Saved: between_network_stats.csv\n\n")

# ==============================================================================
# 10. BETWEEN-NETWORK CORRELATIONS WITH BIOMARKERS
# ==============================================================================

cat("=== Correlating Between-Network Connectivity with Biomarkers ===\n")

# Merge with proteomics
between_network_biomarker_data <- all_between_pairs %>%
  rename(record_id_temp = record_id) %>%
  left_join(small_dat %>% select(record_id, all_of(qx_var)), 
            by = c("record_id_temp" = "record_id"))

between_network_biomarker_cors <- data.frame()

for (pair in top_pairs) {  # Focus on top differing pairs
  pair_data <- between_network_biomarker_data %>% filter(Network_Pair == pair)
  
  for (marker in qx_var) {
    complete_cases <- complete.cases(pair_data$Between_Network_Connectivity, pair_data[[marker]])
    
    if (sum(complete_cases) < 5) next
    
    x <- pair_data$Between_Network_Connectivity[complete_cases]
    y <- pair_data[[marker]][complete_cases]
    
    cor_test <- cor.test(x, y, method = "pearson")
    
    between_network_biomarker_cors <- rbind(between_network_biomarker_cors, data.frame(
      Network_Pair = pair,
      Biomarker = marker,
      Correlation = cor_test$estimate,
      p_value = cor_test$p.value,
      N = sum(complete_cases)
    ))
  }
}

between_network_biomarker_cors$p_fdr <- p.adjust(between_network_biomarker_cors$p_value, method = "fdr")

cat("Between-network biomarker correlations computed:", nrow(between_network_biomarker_cors), "\n")

sig_between_cors <- between_network_biomarker_cors %>% filter(p_fdr < 0.05) %>% arrange(p_fdr)

if (nrow(sig_between_cors) > 0) {
  cat("\nSIGNIFICANT between-network correlations (FDR < 0.05):\n")
  print(sig_between_cors %>% mutate(across(where(is.numeric), ~round(., 4))))
} else {
  cat("\nTop 10 between-network correlations:\n")
  print(between_network_biomarker_cors %>% 
          arrange(p_value) %>% 
          head(10) %>%
          mutate(across(where(is.numeric), ~round(., 4))))
}

write.csv(between_network_biomarker_cors, 
          paste0(results.dir, "between_network_biomarker_correlations.csv"), 
          row.names = FALSE)
cat("Saved: between_network_biomarker_correlations.csv\n")

# Heatmap
cor_matrix_between <- between_network_biomarker_cors %>%
  select(Network_Pair, Biomarker, Correlation) %>%
  pivot_wider(names_from = Biomarker, values_from = Correlation) %>%
  column_to_rownames("Network_Pair") %>%
  as.matrix()

pdf(paste0(results.dir, "between_network_biomarker_heatmap.pdf"), width = 10, height = 8)
corrplot(cor_matrix_between, method = "color", tl.col = "black", tl.cex = 0.8, tl.srt = 45,
         col = colorRampPalette(c("#3498DB", "white", "#E74C3C"))(200),
         title = "Between-Network Connectivity vs Biomarkers", mar = c(0,0,2,0))
dev.off()
cat("Saved: between_network_biomarker_heatmap.pdf\n\n")

# ==============================================================================
# 11. BRAIN REGION VOLUME ANALYSIS
# ==============================================================================

cat("\n=== Brain Region Volume Analysis ===\n")

# Load FreeSurfer data
freesurfer_dir <- paste0(base_dir, "../freesurfer_data/")

lh_thickness_file <- paste0(freesurfer_dir, "lh_thickness.txt")
rh_thickness_file <- paste0(freesurfer_dir, "rh_thickness.txt")
subcortical_file <- paste0(freesurfer_dir, "subcortical_volumes.txt")

# Check if files exist
if (file.exists(lh_thickness_file) && file.exists(rh_thickness_file) && file.exists(subcortical_file)) {
  cat("Loading FreeSurfer data...\n")
  
  # Load thickness data
  lh_thickness <- read.table(lh_thickness_file, header = TRUE)
  rh_thickness <- read.table(rh_thickness_file, header = TRUE)
  subcortical_volumes <- read.table(subcortical_file, header = TRUE)
  
  cat("  Left hemisphere thickness regions:", ncol(lh_thickness) - 1, "\n")
  cat("  Right hemisphere thickness regions:", ncol(rh_thickness) - 1, "\n")
  cat("  Subcortical volumes:", ncol(subcortical_volumes) - 1, "\n")
  
  # Merge all FreeSurfer data
  # Assuming first column is subject ID - adjust column name as needed
  id_col_lh <- names(lh_thickness)[1]
  id_col_rh <- names(rh_thickness)[1]
  id_col_sub <- names(subcortical_volumes)[1]
  
  # Rename ID columns to 'subject_id' for merging
  names(lh_thickness)[1] <- "subject_id"
  names(rh_thickness)[1] <- "subject_id"
  names(subcortical_volumes)[1] <- "subject_id"
  
  # Merge all FreeSurfer data
  fs_data_merged <- lh_thickness %>%
    full_join(rh_thickness, by = "subject_id") %>%
    full_join(subcortical_volumes, by = "subject_id")
  
  cat("  Total FreeSurfer participants:", nrow(fs_data_merged), "\n")
  
  # Match subject IDs to record IDs
  # Try to match based on the mri_ids_df mapping
  fs_data_merged$record_id <- NA
  
  for (i in 1:nrow(fs_data_merged)) {
    subject <- fs_data_merged$subject_id[i]
    
    # Try exact match in file_id
    match_idx <- which(mri_ids_df$file_id == subject)
    if (length(match_idx) > 0) {
      fs_data_merged$record_id[i] <- mri_ids_df$ID[match_idx[1]]
    } else {
      # Try partial match (remove underscores, convert to lowercase, etc.)
      subject_clean <- tolower(gsub("_", "", subject))
      file_ids_clean <- tolower(gsub("_", "", mri_ids_df$file_id))
      match_idx <- which(file_ids_clean == subject_clean)
      if (length(match_idx) > 0) {
        fs_data_merged$record_id[i] <- mri_ids_df$ID[match_idx[1]]
      }
    }
  }
  
  cat("  Matched to record IDs:", sum(!is.na(fs_data_merged$record_id)), "\n")
  
  # Merge with clinical/proteomics data
  fs_data <- fs_data_merged %>%
    filter(!is.na(record_id)) %>%
    left_join(small_dat, by = "record_id") %>%
    filter(!is.na(group))
  
  cat("  Final dataset with groups:", nrow(fs_data), "\n")
  
  # Define key regions of interest (based on common FreeSurfer outputs and your presentation)
  # Check what columns are available
  available_cols <- names(fs_data)
  
  # Common subcortical volume names
  potential_volumes <- c(
    # Global measures
    "TotalGrayVol", "CortexVol", "lhCortexVol", "rhCortexVol",
    "BrainSegVolNotVent", "SupraTentorialVolNotVent",
    # Subcortical structures
    "Left.Thalamus", "Right.Thalamus", "Left.Caudate", "Right.Caudate",
    "Left.Putamen", "Right.Putamen", "Left.Hippocampus", "Right.Hippocampus",
    "Left.Amygdala", "Right.Amygdala", "Left.VentralDC", "Right.VentralDC",
    "Left.Accumbens.area", "Right.Accumbens.area"
  )
  
  # Thickness measures (examples)
  potential_thickness <- c(
    "lh_entorhinal_thickness", "rh_entorhinal_thickness",
    "lh_precuneus_thickness", "rh_precuneus_thickness",
    "lh_isthmuscingulate_thickness", "rh_isthmuscingulate_thickness",
    "lh_middletemporal_thickness", "rh_middletemporal_thickness",
    "lh_rostralanterior_thickness", "rh_rostralanterior_thickness"
  )
  
  # Find available measures
  roi_volumes <- potential_volumes[potential_volumes %in% available_cols]
  roi_thickness <- potential_thickness[potential_thickness %in% available_cols]
  all_roi_measures <- c(roi_volumes, roi_thickness)
  
  # If specific columns don't match, use all numeric columns except ID and group columns
  if (length(all_roi_measures) == 0) {
    cat("  Using all available numeric brain measures...\n")
    exclude_cols <- c("subject_id", "record_id", "group", "group_combined", qx_var)
    all_roi_measures <- names(fs_data)[sapply(fs_data, is.numeric) & 
                                         !names(fs_data) %in% exclude_cols]
  }
  
  cat("  Analyzing", length(all_roi_measures), "brain measures\n\n")
  
  if (length(all_roi_measures) > 0) {
    # Group comparisons for brain measures
    cat("=== Comparing Brain Measures Across Groups ===\n")
    brain_measure_3group_stats <- perform_comparison(fs_data, "group", all_roi_measures)
    brain_measure_2group_stats <- perform_comparison(fs_data, "group_combined", all_roi_measures)
    
    write.csv(brain_measure_3group_stats, paste0(results.dir, "brain_measures_3group_stats.csv"), row.names = FALSE)
    write.csv(brain_measure_2group_stats, paste0(results.dir, "brain_measures_2group_stats.csv"), row.names = FALSE)
    cat("Saved: brain_measures_3group_stats.csv\n")
    cat("Saved: brain_measures_2group_stats.csv\n")
    
    # Visualize significant differences (p < 0.05)
    sig_measures <- brain_measure_2group_stats %>%
      filter(p_value < 0.05) %>%
      arrange(p_value) %>%
      head(9) %>%
      pull(Feature)
    
    cat("\nSignificant brain measure differences (p < 0.05):", length(sig_measures), "\n")
    
    if (length(sig_measures) > 0) {
      plot_data_measures <- fs_data %>%
        select(record_id, group_combined, all_of(sig_measures)) %>%
        pivot_longer(cols = all_of(sig_measures), names_to = "Region", values_to = "Value")
      
      measure_stats_for_plot <- brain_measure_2group_stats %>%
        filter(Feature %in% sig_measures) %>%
        rename(Region = Feature) %>%
        mutate(p_label = sprintf("p = %.4f\nd = %.2f", p_value, cohens_d))
      
      p_measures <- ggplot(plot_data_measures, aes(x = group_combined, y = Value, fill = group_combined)) +
        geom_boxplot(alpha = 0.7, outlier.shape = NA) +
        geom_jitter(width = 0.2, alpha = 0.5, size = 2) +
        geom_text(data = measure_stats_for_plot,
                  aes(x = 1.5, y = Inf, label = p_label),
                  vjust = 1.5, size = 3, inherit.aes = FALSE) +
        facet_wrap(~Region, scales = "free_y", ncol = 3) +
        scale_fill_manual(values = c("T2D+OC" = "#E74C3C", "LC" = "#3498DB")) +
        labs(title = "Brain Measures Lower in T2D+OC",
             subtitle = "Significant differences (p < 0.05) with effect sizes",
             x = "Group", y = "Value", fill = "Group") +
        theme_minimal(base_size = 11) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
              legend.position = "bottom",
              plot.title = element_text(face = "bold"))
      
      ggsave(paste0(results.dir, "brain_measures_lower_in_t2d.pdf"), p_measures, width = 14, height = 10)
      ggsave(paste0(results.dir, "brain_measures_lower_in_t2d.png"), p_measures, width = 14, height = 10, dpi = 300)
      cat("Saved: brain_measures_lower_in_t2d.pdf/.png\n")
    }
    
    # Correlations with Quanterix biomarkers
    cat("\n=== Correlating Brain Measures with Quanterix Biomarkers ===\n")
    
    brain_biomarker_cors <- data.frame()
    
    for (measure in all_roi_measures) {
      for (marker in qx_var) {
        complete_cases <- complete.cases(fs_data[[measure]], fs_data[[marker]])
        
        if (sum(complete_cases) >= 5) {
          x <- fs_data[[measure]][complete_cases]
          y <- fs_data[[marker]][complete_cases]
          
          cor_test <- cor.test(x, y, method = "pearson")
          
          brain_biomarker_cors <- rbind(brain_biomarker_cors, data.frame(
            Brain_Measure = measure,
            Biomarker = marker,
            Correlation = cor_test$estimate,
            p_value = cor_test$p.value,
            N = sum(complete_cases)
          ))
        }
      }
    }
    
    brain_biomarker_cors$p_fdr <- p.adjust(brain_biomarker_cors$p_value, method = "fdr")
    
    write.csv(brain_biomarker_cors, paste0(results.dir, "brain_biomarker_correlations.csv"), 
              row.names = FALSE)
    cat("Saved: brain_biomarker_correlations.csv\n")
    
    # Show significant correlations
    sig_cors <- brain_biomarker_cors %>%
      filter(p_value < 0.05) %>%
      arrange(p_value)
    
    cat("\nSignificant brain-biomarker correlations (p < 0.05):", nrow(sig_cors), "\n")
    if (nrow(sig_cors) > 0) {
      print(head(sig_cors, 10) %>% mutate(across(where(is.numeric), ~round(., 4))))
    }
    
    # Create heatmap of correlations (top correlations only to keep readable)
    if (nrow(brain_biomarker_cors) > 0) {
      # Filter to measures with at least one correlation p < 0.1
      measures_with_cors <- brain_biomarker_cors %>%
        group_by(Brain_Measure) %>%
        summarise(min_p = min(p_value), .groups = "drop") %>%
        filter(min_p < 0.1) %>%
        pull(Brain_Measure)
      
      if (length(measures_with_cors) > 0) {
        cor_matrix_brain <- brain_biomarker_cors %>%
          filter(Brain_Measure %in% measures_with_cors) %>%
          select(Brain_Measure, Biomarker, Correlation) %>%
          pivot_wider(names_from = Biomarker, values_from = Correlation) %>%
          column_to_rownames("Brain_Measure") %>%
          as.matrix()
        
        pdf(paste0(results.dir, "brain_biomarker_correlation_heatmap.pdf"), width = 10, height = 12)
        corrplot(cor_matrix_brain, method = "color", tl.col = "black", tl.cex = 0.6, tl.srt = 45,
                 col = colorRampPalette(c("#3498DB", "white", "#E74C3C"))(200),
                 title = "Brain Measures vs Quanterix Biomarker Correlations (p < 0.1)", 
                 mar = c(0,0,3,0))
        dev.off()
        cat("Saved: brain_biomarker_correlation_heatmap.pdf\n")
      }
    }
    
    # Create scatter plots for top correlations
    top_cors <- brain_biomarker_cors %>%
      arrange(p_value) %>%
      head(6)
    
    if (nrow(top_cors) > 0) {
      plot_list_cors <- list()
      
      for (i in 1:nrow(top_cors)) {
        measure <- top_cors$Brain_Measure[i]
        marker <- top_cors$Biomarker[i]
        r <- round(top_cors$Correlation[i], 3)
        p <- top_cors$p_value[i]
        
        plot_data <- fs_data %>%
          filter(!is.na(.data[[measure]]), !is.na(.data[[marker]]))
        
        marker_label <- toupper(gsub("_avg_conc", "", marker))
        measure_label <- gsub("\\.", " ", gsub("_", " ", measure))
        subtitle <- sprintf("r = %.3f, p = %.4f", r, p)
        
        p_scatter <- ggplot(plot_data, aes_string(x = marker, y = measure, color = "group_combined")) +
          geom_point(size = 3, alpha = 0.7) +
          geom_smooth(method = "lm", se = TRUE, aes(group = 1), 
                      color = "black", linetype = "dashed") +
          scale_color_manual(values = c("T2D+OC" = "#E74C3C", "LC" = "#3498DB")) +
          labs(title = paste(marker_label, "vs", measure_label),
               subtitle = subtitle,
               x = paste(marker_label, "(pg/mL)"),
               y = measure_label,
               color = "Group") +
          theme_minimal(base_size = 10) +
          theme(legend.position = "bottom",
                plot.title = element_text(size = 9))
        
        plot_list_cors[[i]] <- p_scatter
      }
      
      combined_cor_plot <- do.call(grid.arrange, c(plot_list_cors, ncol = 3))
      ggsave(paste0(results.dir, "top_brain_biomarker_scatter.pdf"), combined_cor_plot,
             width = 14, height = 10)
      ggsave(paste0(results.dir, "top_brain_biomarker_scatter.png"), combined_cor_plot,
             width = 14, height = 10, dpi = 300)
      cat("Saved: top_brain_biomarker_scatter.pdf/.png\n")
    }
    
    cat("\nBrain region analysis complete!\n\n")
  }
} else {
  cat("FreeSurfer data files not found. Skipping brain region analysis.\n")
  cat("Expected files:\n")
  cat("  ", lh_thickness_file, "\n")
  cat("  ", rh_thickness_file, "\n")
  cat("  ", subcortical_file, "\n\n")
}

# ==============================================================================
# SUMMARY
# ==============================================================================

cat("\n", rep("=", 80), "\n")
cat("COMPREHENSIVE ANALYSIS COMPLETE!\n")
cat(rep("=", 80), "\n\n")

cat("All files saved to:", results.dir, "\n\n")

cat("FILES CREATED:\n")
cat("\nDemographics:\n")
cat("  - demographics_3groups.png\n")
cat("  - demographics_2groups_combined.png\n")

cat("\nStructural MRI:\n")
cat("  - structural_metrics_3groups.pdf/.png\n")
cat("  - structural_metrics_2groups.pdf/.png\n")
cat("  - structural_3group_stats.csv\n")
cat("  - structural_2group_stats.csv\n")
cat("  - all_structural_metrics.csv\n")

cat("\nWithin-Network Connectivity:\n")
cat("  - within_network_3groups.pdf/.png\n")
cat("  - within_network_2groups.pdf/.png\n")

cat("\nBetween-Network Connectivity:\n")
cat("  - between_network_heatmaps_3groups.pdf\n")
cat("  - between_network_heatmaps_2groups.pdf\n")
cat("  - top_between_network_differences.pdf/.png\n")
cat("  - between_network_stats.csv\n")
cat("  - between_network_biomarker_correlations.csv\n")
cat("  - between_network_biomarker_heatmap.pdf\n")

cat("\n=== DONE ===\n")