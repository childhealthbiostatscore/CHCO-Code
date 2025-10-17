fit_models_emm <- function(outcomes, 
                           formula_rhs, 
                           emm_var = "bmi_cat",
                           data, 
                           adjust = "tukey",
                           ...) {
  
  # Load required packages
  require(emmeans)
  require(dplyr)
  require(purrr)
  
  # Handle formula_rhs input - convert to character if needed
  if (inherits(formula_rhs, "formula")) {
    formula_rhs <- deparse(formula_rhs)
    # Remove leading ~ if present
    formula_rhs <- gsub("^\\s*~\\s*", "", formula_rhs)
  } else {
    # Remove leading ~ if present in string
    formula_rhs <- gsub("^\\s*~\\s*", "", formula_rhs)
  }
  
  # If outcomes is a named list (for custom naming), extract names and values
  if (is.list(outcomes)) {
    outcome_names <- names(outcomes)
    outcome_vars <- unlist(outcomes)
  } else {
    # If it's a vector, use the values as both names and variables
    outcome_vars <- outcomes
    outcome_names <- outcomes
  }
  
  # Create formulas and fit models
  models <- setNames(
    lapply(outcome_vars, function(y) {
      formula_str <- paste(y, "~", formula_rhs)
      lm(as.formula(formula_str), data = data)
    }),
    outcome_names
  )
  
  # Compute EMMs and contrasts for each model
  emm_results <- imap(models, function(model, model_name) {
    # Create formula for emmeans (handles both quoted and unquoted variable names)
    emm_formula <- as.formula(paste("~", emm_var))
    emm <- emmeans(model, emm_formula, ...)
    
    list(
      emms = as.data.frame(emm) %>% 
        mutate(outcome = model_name),
      contrasts = as.data.frame(pairs(emm, adjust = adjust)) %>% 
        mutate(outcome = model_name)
    )
  })
  
  # Combine results into tidy tables
  emm_table <- bind_rows(lapply(emm_results, `[[`, "emms")) %>%
    relocate(outcome) %>%
    arrange(outcome, !!sym(emm_var))
  
  contrast_table <- bind_rows(lapply(emm_results, `[[`, "contrasts")) %>%
    relocate(outcome) %>%
    arrange(outcome, contrast)
  
  # Return both tables
  list(
    emm_table = emm_table,
    contrast_table = contrast_table,
    models = models  # Also return models for potential further analysis
  )
}


library(tidyverse)
library(ggplot2)
library(emmeans)
library(lme4)
library(broom)

plot_mean_ci_stars_renal <- function(data, 
                                     y_var, 
                                     y_axis_title,
                                     x_var,
                                     group_var = "uacr_group",  # Grouping variable
                                     group_var_title = "CKD Category",
                                     legend_position = c(0.7, 0.9),
                                     studys_to_plot = c("RENAL-HEIR", "RENAL-HEIRitage"),
                                     test_method = "t.test",  # Can be "t.test" or "wilcox"
                                     paired = TRUE,  # For paired tests if same subjects
                                     show_individual_points = FALSE,
                                     colors = c("#9dbebb", "#468189")
                                     ) {
  
  dodge_val <- 0.15
  y_sym <- rlang::ensym(y_var)
  x_sym <- rlang::ensym(x_var)
  x_name <- rlang::as_name(x_sym)
  y_name <- rlang::as_name(y_sym)
  group_sym <- rlang::ensym(group_var)
  
  # Ensure study is a factor with correct order
  data <- data %>%
    mutate(x = factor(!!x_sym, levels = studys_to_plot))
  
  # Get unique group levels and filter NAs
  data <- data %>% filter(!is.na(!!group_sym))
  group_levels <- unique(data[[group_var]])
  
  # Check if we have enough groups
  if (length(group_levels) < 2) {
    stop(paste("Need at least 2 groups for comparison. Found:", paste(group_levels, collapse=", ")))
  }
  if (length(group_levels) > 2) {
    warning(paste("Function expects exactly 2 groups for comparison. Found:", paste(group_levels, collapse=", "), 
                  "\nUsing first 2 groups:", group_levels[1], "and", group_levels[2]))
    group_levels <- group_levels[1:2]
    data <- data %>% filter(!!group_sym %in% group_levels)
  }
  
  # Print group information
  cat("\nGroups being compared:", group_levels[1], "vs", group_levels[2], "\n")
  cat("Sample sizes by group and study:\n")
  print(data %>% 
          group_by(!!group_sym, !!x_sym) %>% 
          summarise(n = sum(!is.na(!!y_sym)), .groups = "drop") %>%
          pivot_wider(names_from = x_name, values_from = n))
  cat("\n")
  
  # --- Calculate Mean ± 95% CI ---
  mean_dat <- data %>%
    group_by(!!group_sym, x) %>%
    dplyr::summarise(
      n = sum(!is.na(!!y_sym)),
      mean_y = mean(!!y_sym, na.rm = TRUE),
      sd_y = sd(!!y_sym, na.rm = TRUE),
      se_y = sd_y / sqrt(n),
      t_crit = qt(0.975, df = n - 1),
      lower = mean_y - (t_crit * se_y),
      upper = mean_y + (t_crit * se_y),
      .groups = "drop"
    )
  
  # --- Statistical Testing ---
  contrast_results <- purrr::map_dfr(studys_to_plot, function(v) {
    tryCatch({
      study_data <- data %>% filter(x == v)
      
      # Get data for each group
      group1_data <- study_data %>% 
        filter(!!group_sym == group_levels[1]) %>% 
        pull(!!y_sym) %>% 
        na.omit()
      
      group2_data <- study_data %>% 
        filter(!!group_sym == group_levels[2]) %>% 
        pull(!!y_sym) %>% 
        na.omit()
      
      # Check if we have enough data
      if (length(group1_data) < 2 || length(group2_data) < 2) {
        warning(paste("Not enough observations for statistical test at study", v,
                      "\n  Group", group_levels[1], "n =", length(group1_data),
                      "\n  Group", group_levels[2], "n =", length(group2_data)))
        return(tibble(
          x = v,
          p.value = NA_real_,
          estimate = NA_real_,
          conf.low = NA_real_,
          conf.high = NA_real_,
          method = "Insufficient data",
          contrast = paste(group_levels[1], "-", group_levels[2])
        ))
      }
      
      if (test_method == "t.test") {
        # Perform t-test
        if (paired && length(group1_data) == length(group2_data)) {
          test_result <- t.test(group1_data, group2_data, paired = TRUE)
        } else {
          if (paired && length(group1_data) != length(group2_data)) {
            warning("Unequal group sizes at study ", v, ". Using unpaired test.")
          }
          test_result <- t.test(group1_data, group2_data, paired = FALSE)
        }
        
        tibble(
          x = v,
          p.value = test_result$p.value,
          estimate = test_result$estimate[1],  # Mean difference
          conf.low = test_result$conf.int[1],
          conf.high = test_result$conf.int[2],
          method = test_result$method,
          contrast = paste(group_levels[1], "-", group_levels[2])
        )
        
      } else if (test_method == "wilcox") {
        # Perform Wilcoxon test
        if (paired && length(group1_data) == length(group2_data)) {
          test_result <- wilcox.test(group1_data, group2_data, paired = TRUE, conf.int = TRUE)
        } else {
          if (paired) warning("Unequal group sizes at study ", v, ". Using unpaired test.")
          test_result <- wilcox.test(group1_data, group2_data, paired = FALSE, conf.int = TRUE)
        }
        
        tibble(
          study = v,
          p.value = test_result$p.value,
          estimate = ifelse(!is.null(test_result$estimate), test_result$estimate, NA_real_),
          conf.low = ifelse(!is.null(test_result$conf.int), test_result$conf.int[1], NA_real_),
          conf.high = ifelse(!is.null(test_result$conf.int), test_result$conf.int[2], NA_real_),
          method = test_result$method,
          contrast = paste(group_levels[1], "-", group_levels[2])
        )
      }
    }, error = function(e) {
      message(glue::glue("Statistical test failed for study {v}: {e$message}"))
      tibble(x = v, p.value = NA_real_, estimate = NA_real_,
             conf.low = NA_real_, conf.high = NA_real_, 
             method = NA_character_, contrast = NA_character_)
    })
  }) %>%
    dplyr::mutate(
      stars = case_when(
        p.value < 0.001 ~ "***",
        p.value < 0.01  ~ "**",
        p.value < 0.05  ~ "*",
        p.value < 0.1   ~ ".",
        TRUE            ~ ""
      )
    )
  
  # --- Calculate Y position for significance stars ---
  star_positions <- contrast_results %>%
    left_join(
      mean_dat %>% 
        group_by(x) %>%
        dplyr::summarise(
          max_upper = max(upper, na.rm = TRUE),
          .groups = "drop"
        ) %>%
        mutate(max_upper = ifelse(is.infinite(max_upper) | is.na(max_upper), 0, max_upper)),
      by = "x"
    ) %>%
    dplyr::mutate(y_position = max_upper + 0.08 * abs(ifelse(max_upper == 0, 1, max_upper)))
  
  # Print test results for inspection
  cat("\n--- Statistical Test Results ---\n")
  print(contrast_results %>% dplyr::select(x, contrast, p.value, stars, method))
  cat("\n")
  
  # --- Create Plot ---
  p <- ggplot(data, aes(x = x, y = !!y_sym, color = !!group_sym)) +
    geom_hline(yintercept = 0, color = "grey", linetype = "dashed", size = 0.5) +
    # Error bars
    geom_errorbar(data = mean_dat, 
                  aes(x = x, ymin = lower, ymax = upper, group = !!group_sym, color = !!group_sym), 
                  inherit.aes = F,
                  width = 0.1, size = 0.4, position = position_dodge(width = dodge_val)) +
    # Mean lines connecting studys
    geom_line(data = mean_dat, 
              aes(y = mean_y, group = !!group_sym), 
              position = position_dodge(width = dodge_val),
              size = 1.2, alpha = 0.8) +
    # White background points (for contrast)
    geom_point(data = mean_dat, 
               aes(y = mean_y, shape = !!group_sym), 
               size = 10, position = position_dodge(width = dodge_val), 
               color = "white") +
    # Colored points
    geom_point(data = mean_dat, 
               aes(y = mean_y, shape = !!group_sym), 
               size = 6, position = position_dodge(width = dodge_val))
  
  # Optionally add individual data points
  if (show_individual_points) {
    p <- p + 
      geom_point(aes(group = !!group_sym),
                 position = position_jitterdodge(jitter.width = 0.05, dodge.width = dodge_val),
                 alpha = 0.3, size = 3)
  }
  
  # Add significance stars
  p <- p +
    geom_text(data = star_positions, 
              aes(x = x, y = y_position, label = stars), 
              inherit.aes = FALSE,
              size = 7, vjust = 0, fontface = "bold") +
    # Customization
    scale_shape_manual(values = c(16, 17)) +  # Circle and triangle
    scale_color_manual(values = colors) +
    # scale_color_brewer(palette = "Set2") +  # Nice color palette
    scale_x_discrete(expand = expansion(mult = c(0.2, 0.2))) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
    theme_minimal(base_size = 14) +
    theme(
      panel.background = element_rect(fill = "white", color = NA),
      panel.grid = element_blank(),
      legend.position = legend_position,
      legend.justification = c(1, 1),
      legend.background = element_rect(fill = "white", color = NA),
      axis.text = element_text(size = 12)
    ) +
    labs(
      x = NULL, 
      y = y_axis_title,
      color = group_var_title,
      shape = group_var_title)
  
  return(p)
}

# --- Function for plotting delta values specifically ---
plot_delta_comparison <- function(delta_data, 
                                  delta_var,
                                  y_axis_title,
                                  group_var = "uacr_group",
                                  test_method = "t.test") {
  
  delta_sym <- rlang::ensym(delta_var)
  group_sym <- rlang::ensym(group_var)
  
  # Calculate summary statistics
  summary_stats <- delta_data %>%
    group_by(!!group_sym) %>%
    summarise(
      n = sum(!is.na(!!delta_sym)),
      mean_delta = mean(!!delta_sym, na.rm = TRUE),
      sd_delta = sd(!!delta_sym, na.rm = TRUE),
      se_delta = sd_delta / sqrt(n),
      lower = mean_delta - qt(0.975, n-1) * se_delta,
      upper = mean_delta + qt(0.975, n-1) * se_delta,
      .groups = "drop"
    )
  
  # Perform statistical test
  group_levels <- unique(delta_data[[group_var]])
  group1_data <- delta_data %>% filter(!!group_sym == group_levels[1]) %>% pull(!!delta_sym) %>% na.omit()
  group2_data <- delta_data %>% filter(!!group_sym == group_levels[2]) %>% pull(!!delta_sym) %>% na.omit()
  
  if (test_method == "t.test") {
    test_result <- t.test(group1_data, group2_data)
  } else {
    test_result <- wilcox.test(group1_data, group2_data)
  }
  
  p_value <- test_result$p.value
  sig_label <- case_when(
    p_value < 0.001 ~ "***",
    p_value < 0.01 ~ "**",
    p_value < 0.05 ~ "*",
    p_value < 0.1 ~ ".",
    TRUE ~ "ns"
  )
  
  # Create plot
  p <- ggplot(summary_stats, aes(x = !!group_sym, y = mean_delta, fill = !!group_sym)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, size = 0.5) +
    geom_point(size = 4, shape = 21, color = "black") +
    geom_text(data = data.frame(x = 1.5, y = max(summary_stats$upper) * 1.1),
              aes(x = x, y = y), 
              label = paste0("p = ", round(p_value, 4), " ", sig_label),
              inherit.aes = FALSE, size = 4) +
    scale_fill_brewer(palette = "Set2") +
    theme_minimal(base_size = 14) +
    theme(
      legend.position = "none",
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      axis.title = element_text(face = "bold"),
      plot.title = element_text(face = "bold", hjust = 0.5)
    ) +
    labs(
      x = group_var,
      y = paste0("Δ ", y_axis_title, "\n(HEIRitage - HEIR)"),
      title = paste0("Change in ", y_axis_title, " by ", group_var)
    )
  
  return(p)
}

# --- Helper function to check data structure before plotting ---
check_data_structure <- function(data, y_var, group_var = "uacr_group", 
                                 studys_to_check = c("RENAL HEIR", "RENAL HEIRitage")) {
  
  cat("\n=== Data Structure Check ===\n")
  
  # Check if required columns exist
  required_cols <- c("mrn", "study", group_var, y_var)
  missing_cols <- setdiff(required_cols, names(data))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  
  # Check unique values in study column
  cat("\nUnique studys in data:\n")
  print(unique(data$study))
  
  # Check unique values in group variable
  cat("\nUnique values in", group_var, ":\n")
  print(unique(data[[group_var]]))
  
  # Check for NAs
  cat("\nMissing data summary:\n")
  na_summary <- data %>%
    dplyr::select(all_of(c("study", group_var, y_var))) %>%
    summarise(
      study_NAs = sum(is.na(study)),
      group_NAs = sum(is.na(.data[[group_var]])),
      y_var_NAs = sum(is.na(.data[[y_var]])),
      total_rows = n()
    )
  print(na_summary)
  
  # Sample size by group and study
  cat("\nSample sizes by", group_var, "and study:\n")
  sample_sizes <- data %>%
    filter(study %in% studys_to_check) %>%
    group_by(.data[[group_var]], study) %>%
    summarise(
      n_total = n(),
      n_with_y = sum(!is.na(.data[[y_var]])),
      .groups = "drop"
    ) %>%
    arrange(study, .data[[group_var]])
  print(sample_sizes)
  
  # Check if we have paired data
  cat("\nChecking for paired data (same MRNs at both studys):\n")
  mrn_check <- data %>%
    filter(study %in% studys_to_check, !is.na(.data[[y_var]])) %>%
    group_by(mrn) %>%
    summarise(n_studys = n_distinct(study), .groups = "drop")
  
  paired_count <- sum(mrn_check$n_studys == 2)
  single_study <- sum(mrn_check$n_studys == 1)
  
  cat("  MRNs with both studys:", paired_count, "\n")
  cat("  MRNs with only one study:", single_study, "\n")
  
  return(invisible(NULL))
}


library(tidyverse)
library(ggplot2)
library(broom)

plot_mean_ci_continuous_time <- function(data, 
                                         y_var, 
                                         time_var = "years_from_baseline",  # Continuous time variable
                                         y_axis_title,
                                         x_axis_title = "Years",
                                         group_var = "uacr_group",
                                         group_var_title = "CKD Category",
                                         legend_position = c(0.4, 0.9),
                                         test_at_timepoints = NULL,  # Specific timepoints to test
                                         show_individual_points = T,
                                         show_individual_lines = FALSE,
                                         smooth_method = "loess",
                                         colors = c("#f7b538", "#f28482"),
                                         y_buffer = c(0.1, 0.1)) {  # or "lm" for linear
  
  y_sym <- rlang::ensym(y_var)
  y_name <- rlang::as_name(y_sym)
  group_sym <- rlang::ensym(group_var)
  time_sym <- rlang::ensym(time_var)
  
  # Filter out NA values in grouping variable
  data <- data %>% filter(!is.na(!!group_sym))
  
  # Get unique group levels
  group_levels <- unique(data[[group_var]])
  if (length(group_levels) > 2) {
    warning(paste("Found", length(group_levels), "groups. Using first 2."))
    group_levels <- group_levels[1:2]
    data <- data %>% filter(!!group_sym %in% group_levels)
  }
  
  # Calculate mean and CI at each unique timepoint
  mean_dat <- data %>%
    group_by(!!group_sym, !!time_sym) %>%
    dplyr::summarise(
      n = sum(!is.na(!!y_sym)),
      mean_y = mean(!!y_sym, na.rm = TRUE),
      sd_y = sd(!!y_sym, na.rm = TRUE),
      se_y = sd_y / sqrt(n),
      t_crit = qt(0.975, df = n - 1),
      lower = mean_y - (t_crit * se_y),
      upper = mean_y + (t_crit * se_y),
      .groups = "drop"
    ) %>%
    filter(!is.na(mean_y))  # Remove timepoints with no data
  
  # Perform statistical tests at specified timepoints
  if (!is.null(test_at_timepoints)) {
    test_results <- purrr::map_dfr(test_at_timepoints, function(tp) {
      # Find closest actual timepoint in data
      closest_time <- unique(data[[time_var]])[which.min(abs(unique(data[[time_var]]) - tp))]
      
      test_data <- data %>% 
        filter(abs(!!time_sym - closest_time) < 0.01)  # Small tolerance for floating point
      
      group1_data <- test_data %>% 
        filter(!!group_sym == group_levels[1]) %>% 
        pull(!!y_sym) %>% 
        na.omit()
      
      group2_data <- test_data %>% 
        filter(!!group_sym == group_levels[2]) %>% 
        pull(!!y_sym) %>% 
        na.omit()
      
      if (length(group1_data) >= 2 && length(group2_data) >= 2) {
        test_result <- t.test(group1_data, group2_data)
        
        # Get y position for annotation
        y_pos <- mean_dat %>%
          filter(abs(!!time_sym - closest_time) < 0.01) %>%
          summarise(max_y = max(upper, na.rm = TRUE)) %>%
          pull(max_y)
        
        tibble(
          time = closest_time,
          actual_time = tp,
          p.value = test_result$p.value,
          y_position = y_pos * 1.05,
          stars = case_when(
            p.value < 0.001 ~ "***",
            p.value < 0.01 ~ "**",
            p.value < 0.05 ~ "*",
            p.value < 0.1 ~ ".",
            TRUE ~ ""
          )
        )
      } else {
        tibble(time = closest_time, actual_time = tp, p.value = NA, y_position = NA, stars = "")
      }
    })
    
    cat("\n--- Statistical Test Results ---\n")
    print(test_results %>% select(actual_time, p.value, stars))
    cat("\n")
  } else {
    test_results <- NULL
  }
  
  # Create base plot
  p <- ggplot(data, aes(x = !!time_sym, y = !!y_sym, color = !!group_sym)) +
    geom_hline(yintercept = 0, color = "grey", linetype = "dashed", size = 0.5)
  
  # Add individual trajectories if requested
  if (show_individual_lines) {
    p <- p + 
      geom_line(aes(group = mrn, color = !!group_sym), 
                alpha = 0.5, size = 0.7)
  }
  
  # Add individual points if requested
  if (show_individual_points) {
    p <- p + 
      geom_point(size = 4, aes(shape = !!group_sym))
  }
  
  # # Add smooth trend lines with confidence bands
  # if (!is.null(smooth_method) && smooth_method != "none") {
  #   p <- p +
  #     geom_smooth(aes(fill = !!group_sym), 
  #                 method = smooth_method, 
  #                 alpha = 0.2, 
  #                 size = 0.8,
  #                 se = TRUE)
  # }
  # 

  
  # Add significance annotations if tests were performed
  if (!is.null(test_results)) {
    p <- p +
      geom_text(data = test_results %>% filter(stars != ""),
                aes(x = time, y = y_position, label = stars),
                inherit.aes = FALSE,
                size = 6, fontface = "bold")
  }
  
  # Styling
  p <- p +
    # scale_color_brewer(palette = "Set2") +
    # scale_fill_brewer(palette = "Set2") +
    scale_color_manual(values = colors) +
    scale_fill_manual(values = colors) +
    scale_shape_manual(values = c(16, 17)) +
    theme_minimal(base_size = 14) +
    theme(
      panel.background = element_rect(fill = "white", color = NA),
      panel.grid = element_blank(),
      legend.position = legend_position
    ) +
    labs(
      x = x_axis_title,
      y = y_axis_title,
      color = group_var_title,
      fill = group_var_title,
      shape = group_var_title
    ) +
    scale_y_continuous(expand = expansion(mult = y_buffer))
  
  return(p)
}

create_comparison_boxplot <- function(data, 
                                      y_var, 
                                      x_var = "study",
                                      y_label = NULL,
                                      y_buffer = 10,
                                      p_label_size = 8,
                                      colors = c("#a3b18a", "#588157"),
                                      box_color = "#3a5a40",
                                      comparisons = NULL,
                                      paired = FALSE,
                                      pairing_var = "mrn",
                                      show_lines = TRUE,
                                      comparison_method = "t.test",
                                      p_format = "p.signif",
                                      text_size = 15,
                                      jitter_width = 0,
                                      jitter_size = 2.5,
                                      jitter_alpha = 0.5,
                                      line_color = "gray60",
                                      line_size = 0.25,
                                      line_alpha = 0.7) {
  
  # If no y_label provided, use the variable name
  if (is.null(y_label)) {
    y_label <- y_var
  }
  
  # Set default comparisons if none provided
  if (is.null(comparisons)) {
    # Get unique levels of x_var for automatic comparison
    x_levels <- unique(data[[x_var]])
    if (length(x_levels) == 2) {
      comparisons <- list(x_levels)
    } else {
      warning("More than 2 groups detected. Please specify comparisons manually.")
      comparisons <- list()
    }
  }
  
  # Calculate y_max for scaling
  y_max <- max(data[[y_var]], na.rm = TRUE)
  
  # Create the base plot
  p <- data %>%
    ggplot(aes(x = .data[[x_var]], y = .data[[y_var]], fill = .data[[x_var]])) +
    geom_boxplot(size = 1, width = 0.8, color = box_color, outlier.shape = NA) + 
    geom_jitter(size = jitter_size, alpha = jitter_alpha, width = jitter_width) +
    labs(x = NULL, y = y_label) +
    theme(panel.grid = element_blank(),
          panel.background = element_blank(),
          legend.position = "none",
          text = element_text(size = text_size)) +
    scale_y_continuous(limits = c(NA, y_max + y_buffer))
  
  # Add connecting lines if paired and show_lines is TRUE
  if (paired && show_lines && !is.null(pairing_var)) {
    p <- p + 
      geom_line(aes(group = .data[[pairing_var]]), 
                color = line_color, 
                size = line_size, 
                alpha = line_alpha)
  }
  
  # Add statistical comparisons if specified
  if (length(comparisons) > 0) {
    
    p <- p + 
      stat_compare_means(
        comparisons = comparisons,
        paired = paired, 
        method = comparison_method,
        size = p_label_size,
        vjust = 0,
        label = p_format
      )
  }
  
  # Apply colors if provided
  if (!is.null(colors)) {
    if (length(colors) == 1) {
      # Single color for all groups
      p <- p + scale_fill_manual(values = rep(colors, length(unique(data[[x_var]]))))
    } else {
      # Multiple colors
      p <- p + scale_fill_manual(values = colors)
    }
  }
  
  return(p)
}
