# QTL Progress Tracker Functions
# File: qtl_progress_tracker.R
# Description: Functions to track and visualize QTL analysis progress across multiple folders

# Load required libraries
suppressPackageStartupMessages({
  require(tidyverse)
  require(cli)
  require(knitr)
  require(kableExtra)
})

# Function to read analysis parameters and extract phenotype file
read_phenotypes_from_params <- function(params_file) {
  if (!file.exists(params_file)) {
    return(NULL)
  }
  
  params <- tryCatch({
    read.table(params_file, 
               header = FALSE, 
               col.names = c("parameter", "value"),
               sep = "\t",  # Adjust separator if needed
               stringsAsFactors = FALSE)
  }, error = function(e) {
    # Try with different separator if tab doesn't work
    tryCatch({
      read.table(params_file, 
                 header = FALSE, 
                 col.names = c("parameter", "value"),
                 sep = " ",
                 stringsAsFactors = FALSE)
    }, error = function(e2) {
      # Try with = separator
      read.table(params_file, 
                 header = FALSE, 
                 col.names = c("parameter", "value"),
                 sep = "=",
                 stringsAsFactors = FALSE)
    })
  })
  
  # Clean up potential whitespace
  params$parameter <- trimws(params$parameter)
  params$value <- trimws(params$value)
  
  # Find phenotype_file parameter
  phenotype_file_row <- params[params$parameter == "phenotype_file", ]
  
  if (nrow(phenotype_file_row) == 0) {
    return(NULL)
  }
  
  phenotype_file_path <- phenotype_file_row$value[1]
  
  # Handle relative paths
  if (!file.exists(phenotype_file_path)) {
    # Try relative to params file directory
    phenotype_file_path <- file.path(dirname(params_file), phenotype_file_path)
  }
  
  # Read phenotype file to get list of phenotypes
  if (file.exists(phenotype_file_path)) {
    # Assuming phenotype file has phenotypes as column names (adjust as needed)
    phenotypes_data <- tryCatch({
      read.table(phenotype_file_path, 
                 header = TRUE, check.names = FALSE, 
                 nrows = 1)  # Just read header
    }, error = function(e) {
      # Try CSV format
      read.csv(phenotype_file_path, nrows = 1)
    })
    
    # Return column names except the first (assuming first is ID column)
    phenotypes <- colnames(phenotypes_data)[-1]
    
    message(paste(phenotypes, collapse = "; "))
    
    # Clean phenotype names (remove special characters that might not appear in filenames)
    phenotypes <- gsub("[()/]", "_", phenotypes)
    
    return(phenotypes)
  } else {
    warning(paste("Phenotype file not found:", phenotype_file_path))
    return(NULL)
  }
}

# Function to check completion status for a single QTL folder
check_qtl_folder_status <- function(qtl_folder, verbose = FALSE) {
  folder_name <- basename(qtl_folder)
  
  # Check if required subfolders exist
  summary_stats_dir <- file.path(qtl_folder, "summary_statistics")
  manhattan_dir <- file.path(qtl_folder, "manhattan_plots")
  qq_plots_dir <- file.path(qtl_folder, "qq_plots")
  params_file <- file.path(qtl_folder, "analysis_parameters.txt")
  
  # Initialize results
  results <- list(
    folder = folder_name,
    folder_path = qtl_folder,
    phenotypes = character(),
    summary_stats = character(),
    manhattan_plots = character(),
    qq_plots = character(),
    total_phenotypes = 0,
    completed_summary = 0,
    completed_manhattan = 0,
    completed_qq = 0
  )
  
  # Read phenotypes from parameters file
  phenotypes <- read_phenotypes_from_params(params_file)
  
  if (is.null(phenotypes)) {
    if (verbose) cli_alert_warning("Could not read phenotypes from {folder_name}")
    return(results)
  }
  
  results$phenotypes <- phenotypes
  results$total_phenotypes <- length(phenotypes)
  
  if (verbose) {
    cli_alert_info("Found {length(phenotypes)} phenotypes in {folder_name}")
  }
  
  # Check which phenotypes have summary statistics
  if (dir.exists(summary_stats_dir)) {
    summary_files <- list.files(summary_stats_dir, pattern = "\\.(rds|txt|csv|tsv|summary)$")
    for (pheno in phenotypes) {
      # More flexible matching - phenotype name might be part of filename
      if (any(grepl(pheno, summary_files, ignore.case = TRUE))) {
        results$summary_stats <- c(results$summary_stats, pheno)
        results$completed_summary <- results$completed_summary + 1
      }
    }
  }
  
  # Check which phenotypes have Manhattan plots
  if (dir.exists(manhattan_dir)) {
    manhattan_files <- list.files(manhattan_dir, pattern = "\\.(png|pdf|jpg|jpeg|svg|tiff)$")
    for (pheno in phenotypes) {
      if (any(grepl(pheno, manhattan_files, ignore.case = TRUE))) {
        results$manhattan_plots <- c(results$manhattan_plots, pheno)
        results$completed_manhattan <- results$completed_manhattan + 1
      }
    }
  }
  
  # Check which phenotypes have QQ plots
  if (dir.exists(qq_plots_dir)) {
    qq_files <- list.files(qq_plots_dir, pattern = "\\.(png|pdf|jpg|jpeg|svg|tiff)$")
    for (pheno in phenotypes) {
      if (any(grepl(pheno, qq_files, ignore.case = TRUE))) {
        results$qq_plots <- c(results$qq_plots, pheno)
        results$completed_qq <- results$completed_qq + 1
      }
    }
  }
  
  return(results)
}

# Function to create visual progress bar
create_progress_bar <- function(completed, total, label, color = "blue", width = 300) {
  if (total == 0) {
    pct <- 0
  } else {
    pct <- round(100 * completed / total)
  }
  
  bar_colors <- list(
    blue = "#3498db",
    green = "#2ecc71",
    orange = "#f39c12",
    red = "#e74c3c",
    purple = "#9b59b6",
    yellow = "#f1c40f"
  )
  
  bar_color <- ifelse(color %in% names(bar_colors), bar_colors[[color]], bar_colors$blue)
  
  # Create the HTML progress bar
  bar <- paste0(
    "<div style='display: flex; align-items: center; margin: 5px 0;'>",
    "<div style='flex: 0 0 ", width, "px; margin-right: 15px;'>",
    "<div style='background-color: #ecf0f1; border-radius: 10px; overflow: hidden; height: 22px; ",
    "box-shadow: inset 0 1px 3px rgba(0,0,0,0.2);'>",
    "<div style='background: linear-gradient(to bottom, ", bar_color, ", ", 
    adjustcolor(bar_color, red.f = 0.8, green.f = 0.8, blue.f = 0.8), "); ",
    "width: ", pct, "%; height: 100%; ",
    "display: flex; align-items: center; justify-content: center; ",
    "color: white; font-weight: bold; font-size: 12px; ",
    "text-shadow: 1px 1px 1px rgba(0,0,0,0.3); transition: width 0.3s ease;'>",
    ifelse(pct >= 10, paste0(pct, "%"), ""),
    "</div></div></div>",
    "<span style='font-size: 13px; color: #555;'>",
    "<strong>", completed, "/", total, "</strong> ", label,
    "</span></div>"
  )
  
  return(bar)
}

# Function to create a compact progress indicator
create_compact_progress <- function(completed, total, label) {
  pct <- ifelse(total == 0, 0, round(100 * completed / total))
  color <- case_when(
    pct == 100 ~ "✅",
    pct >= 75 ~ "🔵",
    pct >= 50 ~ "🟡",
    pct >= 25 ~ "🟠",
    TRUE ~ "🔴"
  )
  
  return(paste0(color, " ", completed, "/", total, " (", pct, "%) ", label))
}

# Main function to process multiple QTL folders with enhanced output
track_qtl_progress <- function(qtl_folders, 
                               show_missing = TRUE, 
                               show_table = TRUE,
                               compact = FALSE,
                               verbose = FALSE) {
  
  cat("\n")
  cli_rule(left = "QTL Analysis Progress Report", right = format(Sys.Date()))
  cat("\n")
  
  all_results <- list()
  
  # Validate folders
  valid_folders <- qtl_folders[file.exists(qtl_folders)]
  invalid_folders <- qtl_folders[!file.exists(qtl_folders)]
  
  if (length(invalid_folders) > 0) {
    cli_alert_danger("The following folders were not found:")
    for (folder in invalid_folders) {
      cat("  - ", folder, "\n")
    }
    cat("\n")
  }
  
  if (length(valid_folders) == 0) {
    cli_alert_danger("No valid folders to analyze")
    return(invisible(NULL))
  }
  
  # Process each folder
  for (folder in valid_folders) {
    folder_name <- basename(folder)
    
    # Create section header
    cli_h2(folder_name)
    
    # Check status
    status <- check_qtl_folder_status(folder, verbose = verbose)
    all_results[[folder_name]] <- status
    
    if (status$total_phenotypes == 0) {
      cli_alert_warning("No phenotypes found in {folder_name}")
      cat("\n")
      next
    }
    
    # Display progress
    if (compact) {
      # Compact view
      cat(create_compact_progress(status$completed_summary, 
                                  status$total_phenotypes, 
                                  "Summary Statistics"), "\n")
      cat(create_compact_progress(status$completed_manhattan, 
                                  status$total_phenotypes, 
                                  "Manhattan Plots"), "\n")
      cat(create_compact_progress(status$completed_qq, 
                                  status$total_phenotypes, 
                                  "QQ Plots"), "\n")
    } else {
      # Full progress bars
      cat("**Summary Statistics:**\n")
      cat(create_progress_bar(status$completed_summary, 
                              status$total_phenotypes, 
                              "phenotypes completed",
                              color = ifelse(status$completed_summary == status$total_phenotypes, 
                                             "green", "blue")))
      cat("\n\n")
      
      cat("**Manhattan Plots:**\n")
      cat(create_progress_bar(status$completed_manhattan, 
                              status$total_phenotypes, 
                              "phenotypes completed",
                              color = ifelse(status$completed_manhattan == status$total_phenotypes, 
                                             "green", "orange")))
      cat("\n\n")
      
      cat("**QQ Plots:**\n")
      cat(create_progress_bar(status$completed_qq, 
                              status$total_phenotypes, 
                              "phenotypes completed",
                              color = ifelse(status$completed_qq == status$total_phenotypes, 
                                             "green", "purple")))
      cat("\n\n")
    }
    
    # Show missing phenotypes if requested
    if (show_missing) {
      missing_summary <- setdiff(status$phenotypes, status$summary_stats)
      missing_manhattan <- setdiff(status$phenotypes, status$manhattan_plots)
      missing_qq <- setdiff(status$phenotypes, status$qq_plots)
      
      if (length(missing_summary) > 0 || length(missing_manhattan) > 0 || length(missing_qq) > 0) {
        cat("\n<details><summary><strong>Missing Analyses</strong> (click to expand)</summary>\n\n")
        
        if (length(missing_summary) > 0) {
          cat("**Summary Statistics missing for:**\n")
          cat("- ", paste(missing_summary, collapse = "\n- "), "\n\n")
        }
        
        if (length(missing_manhattan) > 0) {
          cat("**Manhattan plots missing for:**\n")
          cat("- ", paste(missing_manhattan, collapse = "\n- "), "\n\n")
        }
        
        if (length(missing_qq) > 0) {
          cat("**QQ plots missing for:**\n")
          cat("- ", paste(missing_qq, collapse = "\n- "), "\n\n")
        }
        
        cat("</details>\n")
      }
    }
    
    cat("\n---\n\n")
  }
  
  # Create summary table if requested
  if (show_table && length(all_results) > 0) {
    cli_h2("Summary Table")
    
    summary_df <- map_df(all_results, ~ {
      data.frame(
        Folder = .x$folder,
        `Total` = .x$total_phenotypes,
        `Summary Stats` = paste0(.x$completed_summary, " (", 
                                 round(100 * .x$completed_summary / max(.x$total_phenotypes, 1)), "%)"),
        `Manhattan` = paste0(.x$completed_manhattan, " (", 
                             round(100 * .x$completed_manhattan / max(.x$total_phenotypes, 1)), "%)"),
        `QQ Plots` = paste0(.x$completed_qq, " (", 
                            round(100 * .x$completed_qq / max(.x$total_phenotypes, 1)), "%)"),
        `Overall` = paste0(round(100 * (.x$completed_summary + .x$completed_manhattan + .x$completed_qq) / 
                                   (3 * max(.x$total_phenotypes, 1))), "%"),
        check.names = FALSE
      )
    })
    
    # Add color coding to overall completion
    summary_df$Overall <- sapply(summary_df$Overall, function(x) {
      pct <- as.numeric(gsub("%", "", x))
      if (pct == 100) {
        paste0("<span style='color: green; font-weight: bold;'>", x, " ✓</span>")
      } else if (pct >= 75) {
        paste0("<span style='color: blue;'>", x, "</span>")
      } else if (pct >= 50) {
        paste0("<span style='color: orange;'>", x, "</span>")
      } else {
        paste0("<span style='color: red;'>", x, "</span>")
      }
    })
    
    summary_df %>%
      kable(format = "html", escape = FALSE, align = c('l', 'c', 'c', 'c', 'c', 'c')) %>%
      kable_styling(bootstrap_options = c("striped", "hover", "condensed"),
                    full_width = FALSE,
                    position = "left") %>%
      column_spec(1, bold = TRUE) %>%
      column_spec(6, bold = TRUE) %>%
      print()
    
    # Add summary statistics
    cat("\n")
    total_phenotypes <- sum(sapply(all_results, function(x) x$total_phenotypes))
    total_summary <- sum(sapply(all_results, function(x) x$completed_summary))
    total_manhattan <- sum(sapply(all_results, function(x) x$completed_manhattan))
    total_qq <- sum(sapply(all_results, function(x) x$completed_qq))
    
    cli_alert_success("Overall: {total_summary}/{total_phenotypes} summary statistics, {total_manhattan}/{total_phenotypes} Manhattan plots, {total_qq}/{total_phenotypes} QQ plots")
  }
  
  return(invisible(all_results))
}

# Helper function to get detailed status for a specific folder
get_folder_details <- function(qtl_folder) {
  status <- check_qtl_folder_status(qtl_folder, verbose = TRUE)
  
  details <- list(
    folder = status$folder,
    path = status$folder_path,
    phenotypes = list(
      total = status$total_phenotypes,
      names = status$phenotypes
    ),
    completion = data.frame(
      Analysis = c("Summary Statistics", "Manhattan Plots", "QQ Plots"),
      Completed = c(status$completed_summary, status$completed_manhattan, status$completed_qq),
      Total = rep(status$total_phenotypes, 3),
      Percentage = c(
        round(100 * status$completed_summary / max(status$total_phenotypes, 1)),
        round(100 * status$completed_manhattan / max(status$total_phenotypes, 1)),
        round(100 * status$completed_qq / max(status$total_phenotypes, 1))
      )
    ),
    missing = list(
      summary = setdiff(status$phenotypes, status$summary_stats),
      manhattan = setdiff(status$phenotypes, status$manhattan_plots),
      qq = setdiff(status$phenotypes, status$qq_plots)
    )
  )
  
  return(details)
}

# Helper function to get detailed status for a specific folder
get_folder_details <- function(qtl_folder) {
  status <- check_qtl_folder_status(qtl_folder, verbose = TRUE)
  
  details <- list(
    folder = status$folder,
    path = status$folder_path,
    phenotypes = list(
      total = status$total_phenotypes,
      names = status$phenotypes
    ),
    completion = data.frame(
      Analysis = c("Summary Statistics", "Manhattan Plots", "QQ Plots"),
      Completed = c(status$completed_summary, status$completed_manhattan, status$completed_qq),
      Total = rep(status$total_phenotypes, 3),
      Percentage = c(
        round(100 * status$completed_summary / max(status$total_phenotypes, 1)),
        round(100 * status$completed_manhattan / max(status$total_phenotypes, 1)),
        round(100 * status$completed_qq / max(status$total_phenotypes, 1))
      )
    ),
    missing = list(
      summary = setdiff(status$phenotypes, status$summary_stats),
      manhattan = setdiff(status$phenotypes, status$manhattan_plots),
      qq = setdiff(status$phenotypes, status$qq_plots)
    )
  )
  
  return(details)
}

# Export functions
message("QTL Progress Tracker loaded successfully!")
message("Main function: track_qtl_progress(qtl_folders, show_missing = TRUE, show_table = TRUE, compact = FALSE, verbose = FALSE, debug = FALSE)")
message("Details function: get_folder_details(qtl_folder)")
message("Debug function: debug_file_matching(qtl_folder)")
message("For debugging, use: track_qtl_progress(qtl_folders, debug = TRUE) or debug_file_matching(qtl_folder)")