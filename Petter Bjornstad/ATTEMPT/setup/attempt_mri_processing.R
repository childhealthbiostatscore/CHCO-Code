# MRI reports to readable files

library(dplyr)
library(readxl)

user <- Sys.info()[["user"]]

if (user == "choiyej") { # local version
  root_path <- "/Users/choiyej/Library/CloudStorage/OneDrive-SharedLibraries-UW/Laura Pyle - Bjornstad/Biostatistics Core Shared Drive"
  git_path <- "/Users/choiyej/GitHub/CHCO-Code/Petter Bjornstad"
} else if (user == "yejichoi") { # hyak version
  root_path <- ""
  git_path <- "/mmfs1/gscratch/togo/yejichoi/CHCO-Code/Petter Bjornstad"
} else if (user == "pylell") {
  root_path <- "/Users/pylell/Library/CloudStorage/OneDrive-SharedLibraries-UW/Bjornstad/Biostatistics Core Shared Drive"
  git_path <- "/Users/pylell/Documents/GitHub/CHCO-Code/Petter Bjornstad"
} else {
  stop("Unknown user: please specify root path for this user.")
}

extract_mri_data <- function(file_path) {
  # Read the Excel file without column names
  df <- read_excel(file_path, col_names = FALSE)
  
  # Initialize results list
  results <- list()
  
  # Add source file name
  results$source_file <- basename(file_path)
  
  # Extract basic information (all in column 3)
  results$participant_id <- as.character(df[2, 3])
  results$date_of_scan <- as.numeric(df[3, 3])
  results$date_mri_read <- as.numeric(df[4, 3])
  results$complete_bold <- as.numeric(df[5, 3])
  results$complete_perfusion <- as.numeric(df[6, 3])
  
  # Check if this file has T1 data by looking for "Complete T1" in column 2
  has_t1 <- any(grepl("Complete T1 MRI data", df[[2]], ignore.case = TRUE), na.rm = TRUE)
  
  if (has_t1) {
    # File with T1 data - adjust row numbers
    results$complete_t1 <- as.numeric(df[7, 3])
    results$complete_diffusion <- as.numeric(df[8, 3])
    results$kidney_scanned <- as.numeric(df[9, 3])
    
    # Staff member might be in row 10, but could be a date
    staff_value <- as.character(df[10, 3])
    results$staff_member <- ifelse(grepl("/", staff_value) || !is.na(suppressWarnings(as.numeric(staff_value))), 
                                   NA, staff_value)
    
    results$date_form_completed <- as.numeric(df[11, 3])
    results$date_form_updated <- as.numeric(df[12, 3])
  } else {
    # File without T1 data - original row numbers
    results$complete_t1 <- NA
    results$complete_diffusion <- as.numeric(df[7, 3])
    results$kidney_scanned <- as.numeric(df[8, 3])
    results$staff_member <- as.character(df[9, 3])
    results$date_form_completed <- as.numeric(df[10, 3])
    results$date_form_updated <- as.numeric(df[11, 3])
  }
  
  # Extract MRI measurements by searching for specific patterns
  # Find rows by looking in column 2 for variable names
  
  # Helper function to find a row by pattern in column 2
  find_measurement_row <- function(pattern) {
    row_idx <- which(grepl(pattern, df[[2]], ignore.case = TRUE))
    if (length(row_idx) > 0) row_idx[1] else NA
  }
  
  # ASL Perfusion 3D_pCASL
  row <- find_measurement_row("Perfusion 3D_pCASL.*cortex")
  if (!is.na(row)) {
    results$asl_perfusion_3d_pcasl_cortex_right <- as.numeric(df[row, 3])
    results$asl_perfusion_3d_pcasl_cortex_left <- as.numeric(df[row, 4])
  }
  
  # ASL Perfusion 2D_pASL
  row <- find_measurement_row("Perfusion 2D_pASL.*cortex")
  if (!is.na(row)) {
    results$asl_perfusion_2d_pasl_cortex_right <- as.numeric(df[row, 3])
    results$asl_perfusion_2d_pasl_cortex_left <- as.numeric(df[row, 4])
  }
  
  # ADC
  row <- find_measurement_row("ADC.*cortex")
  if (!is.na(row)) {
    results$adc_cortex_right <- as.numeric(df[row, 3])
    results$adc_cortex_left <- as.numeric(df[row, 4])
  }
  
  # Kidney volume
  row <- find_measurement_row("kidney volume")
  if (!is.na(row)) {
    results$vibe_kidney_volume_right <- as.numeric(df[row, 3])
    results$vibe_kidney_volume_left <- as.numeric(df[row, 4])
  }
  
  # BOLD-pre measurements
  # Find BOLD-pre row in column 1
  bold_pre_start <- which(df[[1]] == "BOLD-pre")[1]
  if (!is.na(bold_pre_start)) {
    # R2* cortex is on the same row as BOLD-pre
    results$bold_pre_r2_cortex_right <- as.numeric(df[bold_pre_start, 3])
    results$bold_pre_r2_cortex_left <- as.numeric(df[bold_pre_start, 4])
    
    # R2* kidney is the next row
    if (bold_pre_start + 1 <= nrow(df)) {
      results$bold_pre_r2_kidney_right <- as.numeric(df[bold_pre_start + 1, 3])
      results$bold_pre_r2_kidney_left <- as.numeric(df[bold_pre_start + 1, 4])
    }
    
    # R2* medulla is two rows after
    if (bold_pre_start + 2 <= nrow(df)) {
      results$bold_pre_r2_medulla_right <- as.numeric(df[bold_pre_start + 2, 3])
      results$bold_pre_r2_medulla_left <- as.numeric(df[bold_pre_start + 2, 4])
    }
  }
  
  # T1 measurements (if present)
  t1_start <- which(df[[1]] == "T1")[1]
  if (!is.na(t1_start)) {
    # Check for T1 cortex
    t1_cortex_row <- find_measurement_row("T1.*cortex")
    if (!is.na(t1_cortex_row)) {
      results$t1_cortex_right <- as.numeric(df[t1_cortex_row, 3])
      results$t1_cortex_left <- as.numeric(df[t1_cortex_row, 4])
    }
    
    # Check for T1 kidney (not cortex)
    for (i in t1_start:(t1_start + 3)) {
      if (i <= nrow(df) && grepl("T1.*kidney", df[i, 2]) && !grepl("cortex", df[i, 2])) {
        results$t1_kidney_right <- as.numeric(df[i, 3])
        results$t1_kidney_left <- as.numeric(df[i, 4])
        break
      }
    }
  }
  
  # BOLD-post measurements
  bold_post_start <- which(df[[1]] == "BOLD-post")[1]
  if (!is.na(bold_post_start)) {
    # R2* cortex
    val_right <- df[bold_post_start, 3]
    val_left <- df[bold_post_start, 4]
    
    results$bold_post_r2_cortex_right <- ifelse(
      grepl("not acq", as.character(val_right), ignore.case = TRUE) || is.na(val_right),
      NA,
      as.numeric(val_right)
    )
    results$bold_post_r2_cortex_left <- ifelse(
      grepl("not acq", as.character(val_left), ignore.case = TRUE) || is.na(val_left),
      NA,
      as.numeric(val_left)
    )
    
    # R2* kidney (next row)
    if (bold_post_start + 1 <= nrow(df)) {
      val_right <- df[bold_post_start + 1, 3]
      val_left <- df[bold_post_start + 1, 4]
      
      results$bold_post_r2_kidney_right <- ifelse(
        is.na(val_right) || as.character(val_right) == "NA",
        NA,
        as.numeric(val_right)
      )
      results$bold_post_r2_kidney_left <- ifelse(
        is.na(val_left) || as.character(val_left) == "NA",
        NA,
        as.numeric(val_left)
      )
    }
  }
  
  # Convert to data frame
  return(as.data.frame(results))
}

# specify folder path
mri_folder_path <- file.path(root_path, "ATTEMPT/Data Clean/MR2star results folder")

xlsx_files <- list.files(mri_folder_path, pattern = "*.xlsx", full.names = TRUE)

all_mri_data <- map_df(xlsx_files, function(file) {
  cat("Processing:", basename(file), "\n")
  tryCatch({
    extract_mri_data(file)
  }, error = function(e) {
    cat("  Error:", e$message, "\n")
    return(NULL)
  })
})

all_mri_data <- all_mri_data %>%
  dplyr::mutate(
    date_of_scan = as.Date(date_of_scan, origin = "1899-12-30"),
    date_mri_read = as.Date(date_mri_read, origin = "1899-12-30"),
    date_form_completed = as.Date(date_form_completed, origin = "1899-12-30"),
    date_form_updated = as.Date(date_form_updated, origin = "1899-12-30")
  )

colnames(all_mri_data)

rename_mri_columns <- function(data) {
  data %>%
    rename(
      # Date field
      mri_lab_date = date_of_scan,
      
      # Right kidney measurements
      volume_right = vibe_kidney_volume_right,
      bold_r_bl_cortex = bold_pre_r2_cortex_right,
      bold_r_bl_kidney = bold_pre_r2_kidney_right,
      t1_r_cortex = t1_cortex_right,
      t1_r_kidney = t1_kidney_right,
      pcasl3d_right = asl_perfusion_3d_pcasl_cortex_right,
      adc_right = adc_cortex_right,
      
      # Left kidney measurements
      volume_left = vibe_kidney_volume_left,
      bold_l_bl_cortex = bold_pre_r2_cortex_left,
      bold_l_bl_kidney = bold_pre_r2_kidney_left,
      t1_l_cortex = t1_cortex_left,
      t1_l_kidney = t1_kidney_left,
      pcasl3d_left = asl_perfusion_3d_pcasl_cortex_left,
      adc_left = adc_cortex_left
    ) %>%
    # Add completion status based on available data
    mutate(
      study_visit_boldasl_mri_complete = case_when(
        # If key measurements are present, mark as complete
        !is.na(bold_r_bl_cortex) & !is.na(pcasl3d_right) ~ 2,  # Complete
        # If some data exists but not all key measurements
        !is.na(participant_id) & (is.na(bold_r_bl_cortex) | is.na(pcasl3d_right)) ~ 1,  # Unverified
        # Otherwise incomplete
        TRUE ~ 0  # Incomplete
      ),
      # Add placeholder for StO2 if not already present
      sto2_c = NA_real_
    )
}

# Apply the renaming to your data
all_mri_data <- rename_mri_columns(all_mri_data)

# rename participant ID
all_mri_data <- all_mri_data %>%
  dplyr::mutate(
    # Remove ATTEMPT prefix and create subject_id
    subject_id = gsub("^ATTEMPT[[:space:]-]+", "", participant_id, ignore.case = TRUE),
    
    # Extract visit based on suffix
    visit = case_when(
      grepl("-12$", subject_id) ~ "follow_up",
      grepl("-2$", subject_id) ~ "follow_up",
      grepl("-1$", subject_id) ~ "baseline",
      TRUE ~ "baseline"
    ),
    
    # Remove visit suffix from subject_id
    subject_id = gsub("-1$|-2$|-12$", "", subject_id)
  )

# Save the combined data
output_file <-file.path(root_path, "ATTEMPT/Data Clean/all_mri_data_extracted.csv")
write.csv(all_mri_data, output_file, row.names = FALSE, na = "")
cat("\nSaved combined data to:", output_file, "\n")
