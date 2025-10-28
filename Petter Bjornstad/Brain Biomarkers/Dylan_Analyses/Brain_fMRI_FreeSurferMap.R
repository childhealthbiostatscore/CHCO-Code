# ==============================================================================
# FUNCTIONAL CONNECTIVITY USING FREESURFER PARCELLATION
# This uses your FreeSurfer output to define ROIs for fMRI analysis
# ==============================================================================

library(oro.nifti)
library(neurobase)
library(ggplot2)
library(dplyr)
library(corrplot)
library(stringr)

# ==============================================================================
# APPROACH 1: Use FreeSurfer aseg (subcortical) parcellation
# ==============================================================================

# This approach registers fMRI to FreeSurfer space and extracts time series
# from FreeSurfer-defined regions

# However, this requires coregistration between T1 and fMRI
# For a simpler approach that doesn't require perfect registration,
# we can use a standard MNI space atlas

# ==============================================================================
# APPROACH 2: Use Desikan-Killiany Atlas (FreeSurfer's default)
# Implemented in MNI space for easy application to fMRI
# ==============================================================================

# Desikan-Killiany regions (68 cortical + 14 subcortical = 82 total)
desikan_killiany_regions <- c(
  # Left hemisphere cortical (34)
  "L_BanksSTS", "L_CaudalAnteriorCingulate", "L_CaudalMiddleFrontal",
  "L_Cuneus", "L_Entorhinal", "L_Fusiform", "L_InferiorParietal",
  "L_InferiorTemporal", "L_IsthmusCingulate", "L_LateralOccipital",
  "L_LateralOrbitofrontal", "L_Lingual", "L_MedialOrbitofrontal",
  "L_MiddleTemporal", "L_Parahippocampal", "L_Paracentral",
  "L_ParsOpercularis", "L_ParsOrbitalis", "L_ParsTriangularis",
  "L_Pericalcarine", "L_Postcentral", "L_PosteriorCingulate",
  "L_Precentral", "L_Precuneus", "L_RostralAnteriorCingulate",
  "L_RostralMiddleFrontal", "L_SuperiorFrontal", "L_SuperiorParietal",
  "L_SuperiorTemporal", "L_Supramarginal", "L_FrontalPole",
  "L_TemporalPole", "L_TransverseTemporal", "L_Insula",
  
  # Right hemisphere cortical (34)
  "R_BanksSTS", "R_CaudalAnteriorCingulate", "R_CaudalMiddleFrontal",
  "R_Cuneus", "R_Entorhinal", "R_Fusiform", "R_InferiorParietal",
  "R_InferiorTemporal", "R_IsthmusCingulate", "R_LateralOccipital",
  "R_LateralOrbitofrontal", "R_Lingual", "R_MedialOrbitofrontal",
  "R_MiddleTemporal", "R_Parahippocampal", "R_Paracentral",
  "R_ParsOpercularis", "R_ParsOrbitalis", "R_ParsTriangularis",
  "R_Pericalcarine", "R_Postcentral", "R_PosteriorCingulate",
  "R_Precentral", "R_Precuneus", "R_RostralAnteriorCingulate",
  "R_RostralMiddleFrontal", "R_SuperiorFrontal", "R_SuperiorParietal",
  "R_SuperiorTemporal", "R_Supramarginal", "R_FrontalPole",
  "R_TemporalPole", "R_TransverseTemporal", "R_Insula",
  
  # Subcortical (14)
  "L_Thalamus", "L_Caudate", "L_Putamen", "L_Pallidum",
  "L_Hippocampus", "L_Amygdala", "L_Accumbens",
  "R_Thalamus", "R_Caudate", "R_Putamen", "R_Pallidum",
  "R_Hippocampus", "R_Amygdala", "R_Accumbens"
)

# ==============================================================================
# SIMPLIFIED APPROACH: Use Major Brain Regions
# Based on FreeSurfer nomenclature but simplified for fMRI resolution
# ==============================================================================

create_freesurfer_style_parcellation <- function(brain_coords, dims) {
  
  n_voxels <- nrow(brain_coords)
  
  # Normalize coordinates
  x_norm <- (brain_coords[,1] - min(brain_coords[,1])) / 
    (max(brain_coords[,1]) - min(brain_coords[,1]))
  y_norm <- (brain_coords[,2] - min(brain_coords[,2])) / 
    (max(brain_coords[,2]) - min(brain_coords[,2]))
  z_norm <- (brain_coords[,3] - min(brain_coords[,3])) / 
    (max(brain_coords[,3]) - min(brain_coords[,3]))
  
  regions <- character(n_voxels)
  
  for (i in 1:n_voxels) {
    # Hemisphere (left/right)
    hemi <- ifelse(x_norm[i] < 0.5, "L", "R")
    
    # Anterior-Posterior division (4 zones)
    if (y_norm[i] < 0.25) {
      ap_region <- "Frontal"
    } else if (y_norm[i] < 0.5) {
      ap_region <- "CentralParietal"
    } else if (y_norm[i] < 0.75) {
      ap_region <- "Temporal"
    } else {
      ap_region <- "Occipital"
    }
    
    # Superior-Inferior division
    if (z_norm[i] < 0.35) {
      si_region <- "Inferior"
    } else if (z_norm[i] < 0.65) {
      si_region <- "Middle"
    } else {
      si_region <- "Superior"
    }
    
    # Special regions based on location
    if (y_norm[i] < 0.3 && z_norm[i] < 0.4) {
      # Likely orbitofrontal
      region_name <- paste0(hemi, "_Orbitofrontal")
    } else if (y_norm[i] > 0.4 && y_norm[i] < 0.6 && z_norm[i] < 0.4) {
      # Likely hippocampus/amygdala region
      region_name <- paste0(hemi, "_MedialTemporal")
    } else if (y_norm[i] < 0.35 && z_norm[i] > 0.6) {
      # Likely superior frontal
      region_name <- paste0(hemi, "_SuperiorFrontal")
    } else if (y_norm[i] > 0.35 && y_norm[i] < 0.5 && z_norm[i] > 0.5) {
      # Likely parietal
      region_name <- paste0(hemi, "_Parietal")
    } else if (y_norm[i] > 0.75) {
      # Occipital
      region_name <- paste0(hemi, "_Occipital")
    } else {
      # General region
      region_name <- paste0(hemi, "_", ap_region, "_", si_region)
    }
    
    regions[i] <- region_name
  }
  
  return(regions)
}

# ==============================================================================
# COMPUTE FC WITH FREESURFER-STYLE LABELS
# ==============================================================================

compute_fc_freesurfer_style <- function(fmri_file, participant_id) {
  
  tryCatch({
    cat("\n=== Processing:", participant_id, "===\n")
    
    # Load fMRI data
    fmri_data <- readNIfTI(fmri_file, reorient = FALSE)
    dims <- dim(fmri_data)
    
    if (length(dims) != 4) {
      warning("Not 4D data")
      return(NULL)
    }
    
    n_timepoints <- dims[4]
    cat("  Dimensions:", paste(dims, collapse = " x "), "\n")
    cat("  Time points:", n_timepoints, "\n")
    
    # Create brain mask
    mean_img <- apply(fmri_data, c(1,2,3), mean, na.rm = TRUE)
    sd_img <- apply(fmri_data, c(1,2,3), sd, na.rm = TRUE)
    
    threshold_mean <- quantile(mean_img[mean_img > 0], 0.25, na.rm = TRUE)
    threshold_sd <- quantile(sd_img[sd_img > 0], 0.1, na.rm = TRUE)
    
    brain_mask <- (mean_img > threshold_mean) & (sd_img > threshold_sd)
    brain_coords <- which(brain_mask, arr.ind = TRUE)
    
    cat("  Brain voxels:", nrow(brain_coords), "\n")
    
    # Create FreeSurfer-style parcellation
    cat("  Creating FreeSurfer-style parcellation...\n")
    region_labels <- create_freesurfer_style_parcellation(brain_coords, dims)
    
    # Get unique regions
    unique_regions <- unique(region_labels)
    n_regions <- length(unique_regions)
    cat("  Brain regions:", n_regions, "\n")
    
    # Extract time series for each region
    roi_timeseries <- matrix(NA, nrow = n_timepoints, ncol = n_regions)
    colnames(roi_timeseries) <- unique_regions
    
    cat("  Extracting ROI time series...\n")
    for (r in 1:n_regions) {
      region_name <- unique_regions[r]
      region_idx <- region_labels == region_name
      region_voxels <- brain_coords[region_idx, , drop = FALSE]
      
      # Extract and average time series
      region_signals <- matrix(NA, nrow = n_timepoints, ncol = nrow(region_voxels))
      
      for (v in 1:nrow(region_voxels)) {
        x <- region_voxels[v, 1]
        y <- region_voxels[v, 2]
        z <- region_voxels[v, 3]
        region_signals[, v] <- fmri_data[x, y, z, ]
      }
      
      roi_timeseries[, r] <- rowMeans(region_signals, na.rm = TRUE)
    }
    
    # Detrend
    cat("  Preprocessing time series...\n")
    for (r in 1:n_regions) {
      if (!all(is.na(roi_timeseries[, r]))) {
        time_points <- 1:n_timepoints
        lm_fit <- lm(roi_timeseries[, r] ~ time_points)
        roi_timeseries[, r] <- residuals(lm_fit)
      }
    }
    
    # Remove invalid ROIs
    valid_rois <- apply(roi_timeseries, 2, function(x) {
      !all(is.na(x)) && 
        sum(!is.na(x)) > (n_timepoints * 0.8) &&
        var(x, na.rm = TRUE) > 1e-10
    })
    
    roi_timeseries <- roi_timeseries[, valid_rois]
    roi_names <- colnames(roi_timeseries)
    
    cat("  Valid regions:", length(roi_names), "\n")
    
    # Compute functional connectivity
    cat("  Computing functional connectivity...\n")
    fc_matrix <- cor(roi_timeseries, use = "pairwise.complete.obs")
    
    # Clean up
    fc_matrix[is.na(fc_matrix)] <- 0
    fc_matrix[is.infinite(fc_matrix)] <- 0
    diag(fc_matrix) <- 0
    
    # Add row/column names
    rownames(fc_matrix) <- roi_names
    colnames(fc_matrix) <- roi_names
    
    cat("  Final FC matrix:", nrow(fc_matrix), "x", ncol(fc_matrix), "\n")
    cat("  Mean connectivity:", round(mean(fc_matrix, na.rm = TRUE), 3), "\n\n")
    
    return(list(
      fc_matrix = fc_matrix,
      roi_names = roi_names
    ))
    
  }, error = function(e) {
    cat("\n  ERROR:", e$message, "\n")
    return(NULL)
  })
}

# ==============================================================================
# PROCESS ALL PARTICIPANTS
# ==============================================================================

cat("=== FUNCTIONAL CONNECTIVITY WITH FREESURFER-STYLE REGIONS ===\n\n")

base_dir <- '~/Documents/UofW/Projects/Brain_Imaging/fMRI/data/'

# Get all fMRI files
fmri_files <- list.files(base_dir, 
                         pattern = "rest_brain\\.nii\\.gz$",
                         full.names = TRUE)

# Extract participant IDs
participant_ids <- basename(fmri_files) %>%
  str_remove("_rest_brain\\.nii\\.gz$")

cat("Found", length(participant_ids), "participants\n\n")

# Process all participants
all_results <- list()

for (i in 1:length(fmri_files)) {
  result <- compute_fc_freesurfer_style(fmri_files[i], participant_ids[i])
  
  if (!is.null(result)) {
    all_results[[participant_ids[i]]] <- result
  }
}

cat("\n=== Processing Complete ===\n")
cat("Successfully processed:", length(all_results), "participants\n")

# ==============================================================================
# SAVE INDIVIDUAL AND COMBINED FILES
# ==============================================================================

output_dir <- "~/Documents/UofW/Projects/Brain_Imaging/fMRI/results/individual_fc_freesurfer"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

cat("\n=== Saving Individual Files ===\n")

# Get common ROI names
common_roi_names <- all_results[[1]]$roi_names

# Region information
regions_info <- data.frame(
  ROI_Number = 1:length(common_roi_names),
  ROI_Name = common_roi_names,
  Hemisphere = ifelse(grepl("^L_", common_roi_names), "Left", "Right"),
  Region_Type = sapply(strsplit(common_roi_names, "_"), function(x) {
    if (length(x) > 1) paste(x[-1], collapse = "_") else x[1]
  }),
  stringsAsFactors = FALSE
)

# Save each participant
for (participant_id in names(all_results)) {
  
  participant_data <- list(
    participant_id = participant_id,
    fc_matrix = all_results[[participant_id]]$fc_matrix,
    roi_names = all_results[[participant_id]]$roi_names,
    regions = regions_info,
    atlas = "FreeSurfer-style parcellation",
    date_created = Sys.time()
  )
  
  # Save as RDS
  output_file_rds <- file.path(output_dir, paste0(participant_id, "_fc_freesurfer.rds"))
  saveRDS(participant_data, output_file_rds)
  
  # Save as CSV
  output_file_csv <- file.path(output_dir, paste0(participant_id, "_fc_matrix.csv"))
  write.csv(all_results[[participant_id]]$fc_matrix, output_file_csv, row.names = TRUE)
  
  cat("  Saved:", participant_id, "\n")
}

# Save combined file
output_data <- list(
  fc_matrices = lapply(all_results, function(x) x$fc_matrix),
  roi_names = common_roi_names,
  participant_ids = names(all_results),
  regions = regions_info,
  atlas = "FreeSurfer-style parcellation",
  note = "Regions named following FreeSurfer conventions"
)

saveRDS(output_data, "fc_freesurfer_style_all.rds")

# Save region labels
write.csv(regions_info, "brain_regions_freesurfer_style.csv", row.names = FALSE)

cat("\nFiles saved to:", output_dir, "\n")
cat("Combined file: fc_freesurfer_style_all.rds\n")
cat("Region labels: brain_regions_freesurfer_style.csv\n")

# Print regions
cat("\n=== Brain Regions (FreeSurfer-style) ===\n")
print(regions_info)

cat("\n=== Analysis Complete ===\n")






















########## Resting State Networks 
# ==============================================================================
# RESTING-STATE NETWORK ANALYSIS
# Analyze major brain networks: DMN, FPN, SAN, DAN, VAN, SMN, VIS
# ==============================================================================

library(oro.nifti)
library(neurobase)
library(ggplot2)
library(dplyr)
library(tidyr)
library(corrplot)
library(stringr)

# ==============================================================================
# 1. DEFINE RESTING-STATE NETWORKS
# ==============================================================================

# Based on Yeo et al. (2011) 7-network parcellation
# and Power et al. (2011) functional network definitions

# Network definitions based on anatomical regions
define_resting_state_networks <- function(roi_names) {
  
  network_assignments <- data.frame(
    ROI_Name = roi_names,
    Network = NA_character_,
    stringsAsFactors = FALSE
  )
  
  for (i in 1:length(roi_names)) {
    roi <- roi_names[i]
    
    # LIMBIC NETWORK (LIN)
    if (grepl("MedialTemporal|Orbitofrontal", roi, ignore.case = TRUE)) {
      network_assignments$Network[i] <- "LIN"
    }
    
    # VISUAL NETWORK (VIS)
    else if (grepl("Occipital", roi, ignore.case = TRUE)) {
      network_assignments$Network[i] <- "VIS"
    }
    
    # SOMATOMOTOR NETWORK (SMN)
    else if (grepl("CentralParietal", roi, ignore.case = TRUE)) {
      network_assignments$Network[i] <- "SMN"
    }
    
    # DEFAULT MODE NETWORK (DMN)
    else if (grepl("Temporal_Middle|Temporal_Superior", roi, ignore.case = TRUE)) {
      network_assignments$Network[i] <- "DMN"
    }
    
    # FRONTOPARIETAL NETWORK (FPN)
    else if (grepl("Frontal_Middle|SuperiorFrontal|^[LR]_Parietal$", roi, ignore.case = TRUE)) {
      network_assignments$Network[i] <- "FPN"
    }
    
    # SALIENCE NETWORK (SAN)
    else if (grepl("Temporal_Inferior", roi, ignore.case = TRUE)) {
      network_assignments$Network[i] <- "SAN"
    }
    
    else {
      network_assignments$Network[i] <- "Other"
    }
  }
  
  return(network_assignments)
}

# NOW RE-RUN THE ASSIGNMENT
network_assignments <- define_resting_state_networks(roi_names)

# CHECK THE RESULTS
table(network_assignments$Network)

# ==============================================================================
# 2. COMPUTE NETWORK-LEVEL METRICS
# ==============================================================================

compute_network_metrics <- function(fc_matrix, network_assignments) {
  
  # Get unique networks (excluding "Other")
  networks <- unique(network_assignments$Network)
  networks <- networks[networks != "Other" & !is.na(networks)]
  
  # Initialize results
  network_metrics <- data.frame(
    Network = networks,
    Within_Network_Connectivity = NA,
    N_Regions = NA,
    Mean_Degree = NA,
    stringsAsFactors = FALSE
  )
  
  # Compute metrics for each network
  for (net_idx in 1:length(networks)) {
    network_name <- networks[net_idx]
    
    # Get ROIs in this network
    network_rois <- network_assignments$ROI_Name[network_assignments$Network == network_name]
    
    if (length(network_rois) > 1) {
      # Get indices
      roi_indices <- which(rownames(fc_matrix) %in% network_rois)
      
      # Within-network connectivity (average correlation within network)
      if (length(roi_indices) > 1) {
        network_fc <- fc_matrix[roi_indices, roi_indices]
        # Get upper triangle (exclude diagonal)
        within_conn <- network_fc[upper.tri(network_fc)]
        
        network_metrics$Within_Network_Connectivity[net_idx] <- mean(within_conn, na.rm = TRUE)
        network_metrics$N_Regions[net_idx] <- length(roi_indices)
        
        # Mean degree (average number of strong connections per node)
        threshold <- 0.3  # Define "strong" connection
        strong_conn <- abs(network_fc) > threshold
        diag(strong_conn) <- FALSE
        network_metrics$Mean_Degree[net_idx] <- mean(rowSums(strong_conn), na.rm = TRUE)
      }
    }
  }
  
  return(network_metrics)
}

# ==============================================================================
# 3. COMPUTE BETWEEN-NETWORK CONNECTIVITY
# ==============================================================================

compute_between_network_connectivity <- function(fc_matrix, network_assignments) {
  
  networks <- unique(network_assignments$Network)
  networks <- networks[networks != "Other" & !is.na(networks)]
  
  n_networks <- length(networks)
  
  # Initialize between-network connectivity matrix
  between_network_fc <- matrix(NA, nrow = n_networks, ncol = n_networks)
  rownames(between_network_fc) <- networks
  colnames(between_network_fc) <- networks
  
  for (i in 1:n_networks) {
    for (j in 1:n_networks) {
      net_i <- networks[i]
      net_j <- networks[j]
      
      # Get ROIs in each network
      rois_i <- network_assignments$ROI_Name[network_assignments$Network == net_i]
      rois_j <- network_assignments$ROI_Name[network_assignments$Network == net_j]
      
      if (length(rois_i) > 0 && length(rois_j) > 0) {
        # Get indices
        idx_i <- which(rownames(fc_matrix) %in% rois_i)
        idx_j <- which(rownames(fc_matrix) %in% rois_j)
        
        if (i == j) {
          # Within-network (upper triangle only)
          if (length(idx_i) > 1) {
            within_fc <- fc_matrix[idx_i, idx_i]
            between_network_fc[i, j] <- mean(within_fc[upper.tri(within_fc)], na.rm = TRUE)
          }
        } else {
          # Between-network
          between_fc <- fc_matrix[idx_i, idx_j]
          between_network_fc[i, j] <- mean(between_fc, na.rm = TRUE)
        }
      }
    }
  }
  
  return(between_network_fc)
}

# ==============================================================================
# 4. PROCESS ALL PARTICIPANTS WITH NETWORK ANALYSIS
# ==============================================================================

cat("=== RESTING-STATE NETWORK ANALYSIS ===\n\n")

base_dir <- '~/Documents/UofW/Projects/Brain_Imaging/fMRI/data/'

# Get all fMRI files
fmri_files <- list.files(base_dir, 
                         pattern = "rest_brain\\.nii\\.gz$",
                         full.names = TRUE)

participant_ids <- basename(fmri_files) %>%
  str_remove("_rest_brain\\.nii\\.gz$")

cat("Found", length(participant_ids), "participants\n\n")

# Load or create FC matrices (using previous script output)
# Option 1: If you already have FC matrices
if (file.exists("fc_freesurfer_style_all.rds")) {
  cat("Loading existing FC matrices...\n")
  fc_data <- readRDS("fc_freesurfer_style_all.rds")
  all_fc_matrices <- fc_data$fc_matrices
  roi_names <- fc_data$roi_names
  participant_ids <- fc_data$participant_ids
} else {
  stop("Please run the FreeSurfer-style FC analysis first to generate fc_freesurfer_style_all.rds")
}

# Define network assignments
cat("\nAssigning regions to resting-state networks...\n")
network_assignments <- define_resting_state_networks(roi_names)

# Print network composition
cat("\n=== Network Composition ===\n")
network_summary <- network_assignments %>%
  filter(!is.na(Network)) %>%
  group_by(Network) %>%
  summarise(N_Regions = n(), .groups = "drop") %>%
  arrange(desc(N_Regions))
print(network_summary)

cat("\n=== Regions by Network ===\n")
for (net in unique(network_assignments$Network)) {
  if (!is.na(net) && net != "Other") {
    regions <- network_assignments$ROI_Name[network_assignments$Network == net]
    cat("\n", net, "(", length(regions), "regions):\n")
    cat("  ", paste(regions, collapse = ", "), "\n")
  }
}

# ==============================================================================
# 5. COMPUTE NETWORK METRICS FOR ALL PARTICIPANTS
# ==============================================================================

cat("\n=== Computing network metrics for all participants ===\n")

all_network_metrics <- list()
all_between_network_fc <- list()

for (i in 1:length(participant_ids)) {
  participant_id <- participant_ids[i]
  fc_matrix <- all_fc_matrices[[i]]
  
  cat("Processing:", participant_id, "\n")
  
  # Compute network metrics
  network_metrics <- compute_network_metrics(fc_matrix, network_assignments)
  network_metrics$Participant <- participant_id
  all_network_metrics[[participant_id]] <- network_metrics
  
  # Compute between-network connectivity
  between_network_fc <- compute_between_network_connectivity(fc_matrix, network_assignments)
  all_between_network_fc[[participant_id]] <- between_network_fc
}

# Combine network metrics across participants
network_metrics_all <- do.call(rbind, all_network_metrics)

cat("\n=== Network Metrics Summary ===\n")
network_summary_stats <- network_metrics_all %>%
  group_by(Network) %>%
  summarise(
    Mean_Within_Connectivity = mean(Within_Network_Connectivity, na.rm = TRUE),
    SD_Within_Connectivity = sd(Within_Network_Connectivity, na.rm = TRUE),
    Mean_Degree = mean(Mean_Degree, na.rm = TRUE),
    .groups = "drop"
  )
print(network_summary_stats)

# ==============================================================================
# 6. SAVE NETWORK ANALYSIS RESULTS
# ==============================================================================

cat("\n=== Saving Results ===\n")

output_dir <- "~/Documents/UofW/Projects/Brain_Imaging/fMRI/results/network_analysis"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Save network assignments
write.csv(network_assignments, 
          file.path(output_dir, "network_assignments.csv"), 
          row.names = FALSE)

# Save network metrics for all participants
write.csv(network_metrics_all, 
          file.path(output_dir, "network_metrics_all_participants.csv"), 
          row.names = FALSE)

# Save between-network connectivity matrices
saveRDS(list(
  between_network_fc = all_between_network_fc,
  participant_ids = participant_ids,
  networks = unique(network_assignments$Network[network_assignments$Network != "Other"])
), file.path(output_dir, "between_network_connectivity.rds"))

# ==============================================================================
# 7. VISUALIZATIONS
# ==============================================================================

cat("\n=== Creating Visualizations ===\n")

# Plot 1: Average within-network connectivity by network
p1 <- ggplot(network_metrics_all, 
             aes(x = reorder(Network, Within_Network_Connectivity, FUN = median), 
                 y = Within_Network_Connectivity)) +
  geom_boxplot(aes(fill = Network), alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  coord_flip() +
  labs(title = "Within-Network Connectivity by Resting-State Network",
       x = "Network",
       y = "Within-Network Connectivity (Mean FC)") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none")

ggsave(file.path(output_dir, "within_network_connectivity.pdf"), p1, 
       width = 10, height = 6)

# Plot 2: Between-network connectivity heatmap (average across participants)
mean_between_network_fc <- Reduce("+", all_between_network_fc) / length(all_between_network_fc)

pdf(file.path(output_dir, "between_network_connectivity_heatmap.pdf"), 
    width = 10, height = 9)

# Make sure corrplot is loaded
if (!require("corrplot")) {
  install.packages("corrplot")
  library(corrplot)
}

corrplot::corrplot(mean_between_network_fc,
                   method = "color",
                   type = "full",
                   tl.col = "black",
                   tl.srt = 45,
                   tl.cex = 1,
                   addCoef.col = "black",
                   number.cex = 0.8,
                   col = colorRampPalette(c("#3498DB", "white", "#E74C3C"))(200),
                   title = "Between-Network Connectivity (Group Average)",
                   mar = c(0, 0, 2, 0))
dev.off()

# Plot 3: Network size distribution
p3 <- ggplot(network_summary, aes(x = reorder(Network, N_Regions), y = N_Regions)) +
  geom_bar(stat = "identity", aes(fill = Network), alpha = 0.8) +
  coord_flip() +
  labs(title = "Number of Regions per Network",
       x = "Network",
       y = "Number of Brain Regions") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none")

ggsave(file.path(output_dir, "network_size_distribution.pdf"), p3,
       width = 8, height = 6)

cat("\nVisualizations saved to:", output_dir, "\n")

# ==============================================================================
# 8. PREPARE DATA FOR GROUP COMPARISONS AND PROTEOMICS
# ==============================================================================

cat("\n=== Preparing data for integration ===\n")

# Reshape network metrics to wide format (one row per participant)
network_metrics_wide <- network_metrics_all %>%
  select(Participant, Network, Within_Network_Connectivity) %>%
  pivot_wider(names_from = Network, 
              values_from = Within_Network_Connectivity,
              names_prefix = "Network_")

# Save for integration with clinical data
write.csv(network_metrics_wide,
          file.path(output_dir, "network_metrics_wide_format.csv"),
          row.names = FALSE)

cat("Wide format data saved for merging with proteomics\n")

# ==============================================================================
# 9. NETWORK STATISTICS SUMMARY
# ==============================================================================


cat("Networks identified:\n")
for (net in network_summary$Network) {
  n_regions <- network_summary$N_Regions[network_summary$Network == net]
  cat(sprintf("  %-6s: %2d regions\n", net, n_regions))
}

cat("\nFiles created:\n")
cat("1. network_assignments.csv - ROI to network mapping\n")
cat("2. network_metrics_all_participants.csv - Metrics for each participant\n")
cat("3. network_metrics_wide_format.csv - Wide format for analysis\n")
cat("4. between_network_connectivity.rds - Network-to-network FC\n")
cat("5. Visualization PDFs\n")

cat("\n=== Next Steps ===\n")
cat("1. Merge network_metrics_wide_format.csv with your proteomics data\n")
cat("2. Compare network metrics between T2D and LC groups\n")
cat("3. Correlate network connectivity with biomarkers\n")
cat("4. Focus on DMN (associated with Alzheimer's) and LIN (hippocampus/memory)\n")

cat("\n=== Key Networks for Your Study ===\n")
cat("DMN - Default Mode Network: Most studied in AD/cognition\n")
cat("LIN - Limbic Network: Includes hippocampus (memory)\n")
cat("FPN - Frontoparietal Network: Executive function\n")
cat("SAN - Salience/Ventral Attention: Related to diabetes complications\n")







########## Network Comparisons

# ==============================================================================
# NETWORK GROUP COMPARISONS: T2D vs OC vs LC
# ==============================================================================

library(ggplot2)
library(dplyr)
library(tidyr)
library(corrplot)
library(stringr)
library(gridExtra)
library(ggpubr)

# ==============================================================================
# 1. LOAD DATA AND ASSIGN GROUPS
# ==============================================================================

cat("=== Loading Network Analysis Data ===\n\n")
harmonized_data <- read.csv("../OneDrive/UW/Laura Pyle - Biostatistics Core Shared Drive/Data Harmonization/Data Clean/harmonized_dataset.csv", na = '')


dat <- harmonized_data %>%
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



# ==============================================================================
# NETWORK GROUP COMPARISONS: T2D vs OC vs LC
# ==============================================================================

library(ggplot2)
library(dplyr)
library(tidyr)
library(corrplot)
library(stringr)
library(gridExtra)
library(ggpubr)

# Set output directory
output_dir <- "~/Documents/UofW/Projects/Brain_Imaging/fMRI/results/network_analysis"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
cat("=== Output directory:", output_dir, "===\n\n")

# ==============================================================================
# 1. LOAD DATA AND ASSIGN GROUPS
# ==============================================================================

cat("=== Loading Network Analysis Data ===\n\n")
network_metrics_all <- read.csv(file.path(output_dir, "network_metrics_all_participants.csv"))
between_network_data <- readRDS(file.path(output_dir, "between_network_connectivity.rds"))

# Match participants to groups
participant_groups <- small_dat %>% select(record_id, group) %>% distinct()
network_metrics_all$record_id <- mri_ids_df$ID[match(network_metrics_all$Participant, mri_ids_df$file_id)]
network_metrics_all$record_id[which(network_metrics_all$Participant == 'RH2_38_O')] <- 'RH2-38-0'

network_metrics_with_groups <- network_metrics_all %>%
  left_join(participant_groups, by = "record_id") %>%
  filter(!is.na(group))

cat("Sample sizes:\n")
print(table(network_metrics_with_groups$group))

# ==============================================================================
# 2. WITHIN-NETWORK CONNECTIVITY COMPARISONS
# ==============================================================================

cat("\n=== Comparing Within-Network Connectivity ===\n\n")

networks <- unique(network_metrics_with_groups$Network)
within_network_results <- data.frame()
pairwise_results_all <- data.frame()

for (net in networks) {
  net_data <- network_metrics_with_groups %>% filter(Network == net)
  cat("\n--- Network:", net, "---\n")
  
  # Check normality
  normality_results <- data.frame()
  for (grp in unique(net_data$group)) {
    grp_data <- net_data$Within_Network_Connectivity[net_data$group == grp]
    if (length(grp_data) >= 3) {
      shapiro_test <- shapiro.test(grp_data)
      normality_results <- rbind(normality_results, data.frame(
        Network = net, Group = grp,
        Shapiro_W = shapiro_test$statistic,
        Shapiro_p = shapiro_test$p.value,
        Normal = ifelse(shapiro_test$p.value > 0.05, "Yes", "No")
      ))
    }
  }
  cat("Normality tests:\n")
  print(normality_results)
  all_normal <- all(normality_results$Shapiro_p > 0.05, na.rm = TRUE)
  
  # Overall test
  if (length(unique(net_data$group)) == 3) {
    if (all_normal) {
      cat("Using parametric tests\n")
      aov_result <- aov(Within_Network_Connectivity ~ group, data = net_data)
      aov_summary <- summary(aov_result)
      overall_test <- "ANOVA"
      overall_statistic <- aov_summary[[1]][1, 4]
      overall_p <- aov_summary[[1]][1, 5]
      
      # Pairwise t-tests
      pairs <- combn(unique(net_data$group), 2, simplify = FALSE)
      for (pair in pairs) {
        grp1 <- net_data$Within_Network_Connectivity[net_data$group == pair[1]]
        grp2 <- net_data$Within_Network_Connectivity[net_data$group == pair[2]]
        test_result <- t.test(grp1, grp2)
        pooled_sd <- sqrt(((length(grp1)-1)*var(grp1) + (length(grp2)-1)*var(grp2)) / 
                            (length(grp1) + length(grp2) - 2))
        cohens_d <- (mean(grp1) - mean(grp2)) / pooled_sd
        
        pairwise_results_all <- rbind(pairwise_results_all, data.frame(
          Network = net, Group1 = pair[1], Group2 = pair[2],
          Test_Used = "t-test", p_value = test_result$p.value,
          Effect_Size = cohens_d, N_Group1 = length(grp1), N_Group2 = length(grp2)
        ))
      }
    } else {
      cat("Using non-parametric tests\n")
      kw_result <- kruskal.test(Within_Network_Connectivity ~ group, data = net_data)
      overall_test <- "Kruskal-Wallis"
      overall_statistic <- kw_result$statistic
      overall_p <- kw_result$p.value
      
      # Pairwise Mann-Whitney U
      pairs <- combn(unique(net_data$group), 2, simplify = FALSE)
      for (pair in pairs) {
        grp1 <- net_data$Within_Network_Connectivity[net_data$group == pair[1]]
        grp2 <- net_data$Within_Network_Connectivity[net_data$group == pair[2]]
        test_result <- wilcox.test(grp1, grp2, exact = FALSE)
        n1 <- length(grp1)
        n2 <- length(grp2)
        effect_size <- 1 - (2*test_result$statistic) / (n1 * n2)
        
        pairwise_results_all <- rbind(pairwise_results_all, data.frame(
          Network = net, Group1 = pair[1], Group2 = pair[2],
          Test_Used = "Mann-Whitney U", p_value = test_result$p.value,
          Effect_Size = effect_size, N_Group1 = n1, N_Group2 = n2
        ))
      }
    }
    
    # Store overall results
    group_stats <- net_data %>% group_by(group) %>%
      summarise(Mean = mean(Within_Network_Connectivity, na.rm = TRUE),
                SD = sd(Within_Network_Connectivity, na.rm = TRUE), .groups = "drop")
    
    within_network_results <- rbind(within_network_results, data.frame(
      Network = net, Test = overall_test, Statistic = overall_statistic,
      p_value = overall_p, All_Normal = all_normal,
      T2D_mean = group_stats$Mean[group_stats$group == "Type 2 Diabetes"],
      OC_mean = group_stats$Mean[group_stats$group == "Obese Control"],
      LC_mean = group_stats$Mean[group_stats$group == "Lean Control"]
    ))
  }
}

# Print and save results
cat("\n=== Within-Network Connectivity: Overall Tests ===\n")
print(within_network_results %>% mutate(across(where(is.numeric), ~round(., 4))))

cat("\n=== Pairwise Comparisons (uncorrected) ===\n")
print(pairwise_results_all %>% mutate(across(where(is.numeric), ~round(., 4))))

sig_pairwise <- pairwise_results_all %>% filter(p_value < 0.05) %>% arrange(p_value)
cat("\n=== SIGNIFICANT (p < 0.05) ===\n")
if (nrow(sig_pairwise) > 0) print(sig_pairwise %>% mutate(across(where(is.numeric), ~round(., 4))))

write.csv(within_network_results, file.path(output_dir, "network_overall_tests.csv"), row.names = FALSE)
write.csv(pairwise_results_all, file.path(output_dir, "network_pairwise_uncorrected.csv"), row.names = FALSE)

cat("\n=== FILES SAVED TO:", output_dir, "===\n")
cat("- network_overall_tests.csv\n")
cat("- network_pairwise_uncorrected.csv\n")
cat("\n=== ANALYSIS COMPLETE ===\n")

# Add this after the line "cat("\n=== ANALYSIS COMPLETE ===\n")"

# ==============================================================================
# 3. VISUALIZATIONS
# ==============================================================================

cat("\n=== Creating Visualizations ===\n")

# Plot 1: Main comparison with overall p-values
p1 <- ggplot(network_metrics_with_groups, 
             aes(x = Network, y = Within_Network_Connectivity, fill = group)) +
  geom_boxplot(position = position_dodge(width = 0.8), alpha = 0.7) +
  geom_point(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.8),
             alpha = 0.5, size = 2) +
  scale_fill_manual(values = c("Type 2 Diabetes" = "#E74C3C", 
                               "Obese Control" = "#F39C12",
                               "Lean Control" = "#3498DB"),
                    name = "Group") +
  labs(title = "Within-Network Connectivity by Group",
       subtitle = "ANOVA p-values shown above each network",
       x = "Network",
       y = "Within-Network Connectivity") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
        legend.position = "bottom",
        plot.title = element_text(face = "bold"))

# Add p-values above each network
for (i in 1:nrow(within_network_results)) {
  net <- within_network_results$Network[i]
  p_val <- within_network_results$p_value[i]
  p_label <- sprintf("p = %.3f", p_val)
  
  y_pos <- max(network_metrics_with_groups$Within_Network_Connectivity[
    network_metrics_with_groups$Network == net], na.rm = TRUE) * 1.05
  
  p1 <- p1 + annotate("text", 
                      x = which(unique(network_metrics_with_groups$Network) == net),
                      y = y_pos, label = p_label, size = 4, fontface = "bold")
}

ggsave(file.path(output_dir, "within_network_by_group.pdf"), p1, width = 14, height = 7)
ggsave(file.path(output_dir, "within_network_by_group.png"), p1, width = 14, height = 7, dpi = 300)
cat("Saved: within_network_by_group.pdf/.png\n")

# Plot 2: Heatmap of pairwise p-values
p_matrix_data <- pairwise_results_all %>%
  mutate(Comparison = paste0(Group1, " vs\n", Group2))

p2 <- ggplot(p_matrix_data, aes(x = Network, y = Comparison, fill = -log10(p_value))) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.3f", p_value)), size = 3) +
  scale_fill_gradient2(low = "white", mid = "#FFA500", high = "#E74C3C",
                       midpoint = -log10(0.05),
                       name = "-log10(p)") +
  labs(title = "Pairwise Comparison P-values (Uncorrected)",
       subtitle = "Darker = more significant",
       x = "Network",
       y = "Group Comparison") +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(output_dir, "pairwise_pvalue_heatmap.pdf"), p2, width = 10, height = 6)
ggsave(file.path(output_dir, "pairwise_pvalue_heatmap.png"), p2, width = 10, height = 6, dpi = 300)
cat("Saved: pairwise_pvalue_heatmap.pdf/.png\n")

# Plot 3: Individual network panels with significance brackets
plot_list <- list()
for (net in unique(network_metrics_with_groups$Network)) {
  net_data <- network_metrics_with_groups %>% filter(Network == net)
  net_pairwise <- pairwise_results_all %>% filter(Network == net)
  
  p_net <- ggplot(net_data, aes(x = group, y = Within_Network_Connectivity, fill = group)) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA) +
    geom_jitter(width = 0.2, alpha = 0.6, size = 2) +
    scale_fill_manual(values = c("Type 2 Diabetes" = "#E74C3C", 
                                 "Obese Control" = "#F39C12",
                                 "Lean Control" = "#3498DB")) +
    labs(title = net, x = NULL, y = "Within-Network Connectivity") +
    theme_minimal(base_size = 11) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none",
          plot.title = element_text(face = "bold", hjust = 0.5))
  
  # Add p-value annotations for significant comparisons
  if (nrow(net_pairwise) > 0) {
    y_max <- max(net_data$Within_Network_Connectivity, na.rm = TRUE)
    y_range <- diff(range(net_data$Within_Network_Connectivity, na.rm = TRUE))
    
    for (i in 1:nrow(net_pairwise)) {
      if (net_pairwise$p_value[i] < 0.10) {
        y_position <- y_max + y_range * (0.05 + i * 0.08)
        sig_label <- ifelse(net_pairwise$p_value[i] < 0.001, "***",
                            ifelse(net_pairwise$p_value[i] < 0.01, "**",
                                   ifelse(net_pairwise$p_value[i] < 0.05, "*", "†")))
        
        grp1_x <- which(levels(factor(net_data$group)) == net_pairwise$Group1[i])
        grp2_x <- which(levels(factor(net_data$group)) == net_pairwise$Group2[i])
        
        p_net <- p_net +
          annotate("segment", x = grp1_x, xend = grp2_x,
                   y = y_position, yend = y_position) +
          annotate("text", x = mean(c(grp1_x, grp2_x)), y = y_position + y_range * 0.02,
                   label = sprintf("p=%.3f %s", net_pairwise$p_value[i], sig_label), size = 3)
      }
    }
  }
  
  plot_list[[net]] <- p_net
}

# Combine all network plots
combined_plot <- do.call(grid.arrange, c(plot_list, ncol = 3))
ggsave(file.path(output_dir, "within_network_detailed.pdf"), combined_plot, width = 14, height = 10)
ggsave(file.path(output_dir, "within_network_detailed.png"), combined_plot, width = 14, height = 10, dpi = 300)
cat("Saved: within_network_detailed.pdf/.png\n")

cat("\n=== ALL VISUALIZATIONS CREATED ===\n")













################## Brain biomarker associations


library(ggplot2)
library(dplyr)
library(tidyr)
library(corrplot)
library(stringr)
library(gridExtra)
library(ggpubr)

# Set output directory
output_dir <- "~/Documents/UofW/Projects/Brain_Imaging/fMRI/results/network_analysis"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
cat("=== Output directory:", output_dir, "===\n\n")

# ==============================================================================
# 1. LOAD DATA AND ASSIGN GROUPS
# ==============================================================================

harmonized_data <- read.csv("../OneDrive/UW/Laura Pyle - Biostatistics Core Shared Drive/Data Harmonization/Data Clean/harmonized_dataset.csv", na = '')

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


dat <- harmonized_data %>%
  arrange(date_of_screen) %>% 
  dplyr::summarise(across(where(negate(is.numeric)), ~ ifelse(all(is.na(.x)), NA_character_, last(na.omit(.x)))),
                   across(where(is.numeric), ~ ifelse(all(is.na(.x)), NA_real_, mean(na.omit(.x), na.rm=T))),
                   .by = c(record_id, visit))

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



cat("=== Loading Network Analysis Data ===\n\n")
network_metrics_all <- read.csv(file.path(output_dir, "network_metrics_all_participants.csv"))
between_network_data <- readRDS(file.path(output_dir, "between_network_connectivity.rds"))

# Match participants to groups
participant_groups <- small_dat %>% select(record_id, group) %>% distinct()
network_metrics_all$record_id <- mri_ids_df$ID[match(network_metrics_all$Participant, mri_ids_df$file_id)]
network_metrics_all$record_id[which(network_metrics_all$Participant == 'RH2_38_O')] <- 'RH2-38-0'

network_metrics_with_groups <- network_metrics_all %>%
  left_join(participant_groups, by = "record_id") %>%
  filter(!is.na(group))

cat("Sample sizes:\n")
print(table(network_metrics_with_groups$group))


network_w_bbiomarkers <- network_metrics_with_groups %>%
  left_join(small_dat %>% dplyr::select(record_id, 3:9), by='record_id')







########## AAL Mapping 


