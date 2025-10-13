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

plot_emms_with_brackets <- function(
    emm_table,
    contrast_table,
    outcome_var,
    x_var = "bmi_cat",
    y_label = "Estimated marginal mean",
    x_label = NULL,
    use_stars = TRUE,
    y_step_multiplier = 0.08,
    point_size = 3,
    errorbar_width = 0.2,
    bracket_size = 0.5,
    label_size = 5,
    text_size = 15
) {
  
  # Load required packages
  require(ggplot2)
  require(tidyverse)
  require(ggpubr)
  require(rstatix)
  
  # 1) Subset EMMs and contrasts for this outcome
  emms_sub <- emm_table %>%
    filter(outcome == outcome_var)
  
  contr_sub <- contrast_table %>%
    filter(outcome == outcome_var) %>%
    separate(contrast, into = c("group1", "group2"), sep = " - ", remove = FALSE) %>%
    mutate(
      group1 = trimws(group1),
      group2 = trimws(group2),
      # Label: either stars or formatted p-values
      label = if(use_stars) {
        stars.pval(p.value)
      } else {
        rstatix::p_format(p.value, add.p = TRUE)
      }
    )
  
  # Check if data exists
  if(nrow(emms_sub) == 0) {
    stop(paste("No data found for outcome:", outcome_var))
  }
  
  # 2) Y positions for brackets (stack them above the highest CI)
  y_base <- max(emms_sub$upper.CL, na.rm = TRUE)
  y_step <- y_step_multiplier * diff(range(c(emms_sub$lower.CL, emms_sub$upper.CL), na.rm = TRUE))
  
  contr_sub <- contr_sub %>%
    arrange(p.value) %>%  # optional: nicer stacking
    mutate(y.position = y_base + row_number() * y_step)
  
  # 3) Create the plot
  p <- ggplot(emms_sub, aes_string(x = x_var, y = "emmean")) +
    geom_point(size = point_size) +
    geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = errorbar_width) +
    labs(
      x = x_label,
      y = y_label
    ) +
    theme_minimal() + 
    theme(
      text = element_text(size = text_size),
      panel.grid = element_blank(),
      legend.position = "none"
    )
  
  # Add significance brackets if contrasts exist
  if(nrow(contr_sub) > 0) {
    p <- p + stat_pvalue_manual(
      contr_sub,
      label = "label",
      y.position = "y.position",
      xmin = "group1",
      xmax = "group2",
      tip.length = 0.01,
      bracket.size = bracket_size,
      size = label_size
    )
  }
  
  return(p)
}
