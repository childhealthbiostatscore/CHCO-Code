#Functions

#--------------
#geometric mean
#--------------
render.geometric <- function(x, ...) {
  gm <- exp(mean(log(x[x > 0]), na.rm = TRUE))  # geometric mean
  gsd <- exp(sd(log(x[x > 0]), na.rm = TRUE))    # geometric SD
  c("", 
    " " = sprintf("%.2f (%.2f)", gm, gsd))
}


#---------------------------
#collapsing metadata function
#---------------------------
collapse_visit <- function(data, by_vars) {
 data %>%
  summarise(
   across(
    where(is.character),
   ~ if (all(is.na(.x))) NA_character_ else dplyr::last(na.omit(.x))
  ),
 across(
        where(is.factor),
       ~ if (all(is.na(.x))) NA else dplyr::last(na.omit(.x))
    ),
   across(
    where(is.numeric),
   ~ if (all(is.na(.x))) NA_real_ else mean(.x, na.rm = TRUE)
      ),
     .by = all_of(by_vars)
  )
}

make_race_ethnicity <- function(data) {
 data %>%
  mutate(
   race_ethnicity_condensed = case_when(
    race == "White" & ethnicity == "Not Hispanic or Latino" ~ "Not Hispanic or Latino White",
   race == "Black or African American" & ethnicity == "Not Hispanic or Latino" ~ "Not Hispanic or Latino Black",
  ethnicity == "Hispanic or Latino" ~ "Hispanic or Latino",
 TRUE ~ "Not Hispanic or Latino Other"
)
)
}


#--------------------------------
#functions for pvalues and labels
#--------------------------------
fmt_p <- function(p) {
  case_when(
    is.na(p) ~ NA_character_,
    p < 0.001 ~ "<0.001",
    TRUE ~ formatC(p, format = "f", digits = 3)
  )
}

get_pfas_label <- function(x) {
  dplyr::recode(x, !!!pfas_label_map, .default = x)
}

get_outcome_label <- function(outcome) {
  dplyr::recode(outcome, !!!outcome_label_map, .default = outcome)
}


#-------------------
#imputation function
#-------------------
impute_pfas_lod <- function(pfas_dat, lod, pfas_vars, threshold = 2/3) {
 n_participants <- n_distinct(pfas_dat$record_id)

pfas_to_keep <- purrr::keep(pfas_vars, function(pf) {
 pf %in% names(lod) &&
  (sum(pfas_dat[[pf]] < lod[[pf]], na.rm = TRUE) / n_participants) < threshold
 })

pfas_imputed <- pfas_dat %>%
 dplyr::select(record_id, all_of(pfas_to_keep))

 for (pf in pfas_to_keep) {
  pfas_imputed[[pf]] <- ifelse(
   pfas_imputed[[pf]] < lod[[pf]],
  lod[[pf]] / sqrt(2),
  pfas_imputed[[pf]]
)
 }
  list(
  imputed_data = pfas_imputed,
  pfas_used = pfas_to_keep
)
}

#--------------------------
#PFAS spearman correlations
#--------------------------
cor_mtest <- function(mat, method = "spearman") {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat <- matrix(NA, n, n)
  diag(p.mat) <- 0
  colnames(p.mat) <- colnames(mat)
  rownames(p.mat) <- colnames(mat)
  
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- suppressWarnings(
        cor.test(mat[, i], mat[, j], method = method, exact = FALSE)
      )
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  
  p.mat
}


#------------------------------------------
#baseline clinical data comparisons--no PFAS
#-------------------------------------------
run_baseline_group_model <- function(data, outcome, covars = c("age", "sex")) {
  
  rhs <- c("group", covars)
  m0 <- as.formula(paste0(outcome, " ~ ", paste(rhs, collapse = " + ")))
  
  model_data <- data %>%
    dplyr::select(all.vars(m0)) %>%
    tidyr::drop_na()
  
  m1 <- lm(m0, data = model_data)
  
  # omnibus test for group
  group_rows <- broom::tidy(stats::anova(m1)) %>%
    dplyr::filter(term == "group") %>%
    dplyr::mutate(
      outcome = outcome,
      n = nrow(model_data)
    )
  
  # pairwise Tukey-adjusted comparisons
  emm <- emmeans::emmeans(m1, ~ group)
  pairwise_tbl <- summary(
    pairs(emm, adjust = "tukey"),
    infer = TRUE
  ) %>%
    as.data.frame() %>%
    dplyr::mutate(
      outcome = outcome,
      n = nrow(model_data)
    )
  
  list(
    model = m1,
    omnibus = group_rows,
    emmeans = as.data.frame(summary(emm)),
    pairwise = pairwise_tbl
  )
}



#------------------------------------
#panel boxplots of clinical variables
#------------------------------------
plot_baseline_group_boxplot <- function(data,
                                        outcome,
                                        results_list,
                                        y_label = NULL) {
  
  plot_data <- data %>%
    dplyr::select(group, sex, age, all_of(outcome)) %>%
    tidyr::drop_na()
  
  pair_tbl <- results_list[[outcome]]$pairwise %>%
    dplyr::mutate(
      group1 = stringr::str_trim(stringr::word(contrast, 1, sep = " - ")),
      group2 = stringr::str_trim(stringr::word(contrast, 2, sep = " - "))
    )
  
  # y positions for significance bars
  y_max <- max(plot_data[[outcome]], na.rm = TRUE)
  y_min <- min(plot_data[[outcome]], na.rm = TRUE)
  y_step <- 0.08 * (y_max - y_min)
  
  pair_tbl <- pair_tbl %>%
    dplyr::mutate(
      y.position = y_max + seq_len(n()) * y_step,
      p_label = dplyr::case_when(
        p.value < 0.001 ~ "***",
        p.value < 0.01  ~ "**",
        p.value < 0.05  ~ "*",
        TRUE ~ "ns"
      )
    )
  
  omnibus_p <- results_list[[outcome]]$omnibus$p.value[1]
  
  ggplot(plot_data, aes(x = group, y = .data[[outcome]], fill = group)) + 
    scale_fill_viridis_d(option = "C") +
    geom_boxplot(outlier.shape = NA, alpha = 0.7) +
    geom_jitter(width = 0.12, alpha = 0.6, size = 2) +
    ggpubr::stat_pvalue_manual(
      pair_tbl,
      label = "p_label",
      xmin = "group1",
      xmax = "group2",
      y.position = "y.position",
      tip.length = 0.01,
      hide.ns = TRUE
    ) +
    theme_bw() +
    labs(
      title = paste("Baseline", get_outcome_label(outcome), "by group"),
      subtitle = paste0("Omnibus group p = ", fmt_p(omnibus_p), "; pairwise p adjusted by Tukey"),
      x = "Group",
      y = ifelse(is.null(y_label), get_outcome_label(outcome), y_label)
    ) +
    theme(
      legend.position = "none",
      plot.title = element_text(face = "bold")
    )
}


#------------------------------------------------------
#plotting baseline PFAS across groups--individual plots
#------------------------------------------------------
plot_baseline_pfas_group_boxplot <- function(data,
                                             exposure,
                                             results_list) {
  
  plot_data <- data %>%
    dplyr::select(group, sex, age, all_of(exposure)) %>%
    tidyr::drop_na()
  
  pair_tbl <- results_list[[exposure]]$pairwise %>%
    dplyr::mutate(
      group1 = stringr::str_trim(stringr::word(contrast, 1, sep = " - ")),
      group2 = stringr::str_trim(stringr::word(contrast, 2, sep = " - "))
    )
  
  y_max <- max(plot_data[[exposure]], na.rm = TRUE)
  y_min <- min(plot_data[[exposure]], na.rm = TRUE)
  y_step <- 0.08 * (y_max - y_min)
  
  pair_tbl <- pair_tbl %>%
    dplyr::mutate(
      y.position = y_max + seq_len(n()) * y_step,
      p_label = dplyr::case_when(
        p.value < 0.001 ~ "***",
        p.value < 0.01  ~ "**",
        p.value < 0.05  ~ "*",
        TRUE ~ "ns"
      )
    )
  
  omnibus_p <- results_list[[exposure]]$omnibus$p.value[1]
  
  ggplot(plot_data, aes(x = group, y = .data[[exposure]], fill = group)) +
    scale_fill_viridis_d(option = "C") +
    geom_boxplot(outlier.shape = NA, alpha = 0.7) +
    geom_jitter(width = 0.12, alpha = 0.6, size = 2) +
    ggpubr::stat_pvalue_manual(
      pair_tbl,
      label = "p_label",
      xmin = "group1",
      xmax = "group2",
      y.position = "y.position",
      tip.length = 0.01,
      hide.ns = TRUE
    ) +
    theme_bw() +
    labs(
      title = paste("Baseline", get_pfas_label(exposure), "by group"),
      subtitle = paste0("Omnibus group p = ", fmt_p(omnibus_p), "; pairwise p adjusted by Tukey"),
      x = "Group",
      y = get_pfas_label(exposure)
    ) +
    theme(
      legend.position = "none",
      plot.title = element_text(face = "bold")
    )
}



#---------------------------------------------
#Model 1: Baseline Cross Section PFAS Analysis
#---------------------------------------------
run_baseline_pfas_lm <- function(data, outcome, exposure) {
  
  rhs <- c(exposure, "age", "sex")
  
  m0 <- as.formula(paste0(outcome, " ~ ", paste(rhs, collapse = " + ")))
  
  model_data <- data %>%
    dplyr::select(all.vars(m0)) %>%
    tidyr::drop_na()
  
  m1 <- lm(m0, data = model_data)
  
  list(
    model = m1,
    tidy = broom::tidy(m1) %>%
      dplyr::mutate(
        analysis = "baseline_cross_sectional",
        outcome = outcome,
        exposure = exposure,
        model_formula = paste(deparse(m0), collapse = " ")
      )
  )
}



#-------------------------------------------
#plotting cross-sectional data 
#-------------------------------------------
plot_baseline_pfas_set <- function(data,
                                   outcome,
                                   exposures = pfas_group_vars,
                                   results_tbl = baseline_results_final,
                                   save_pdf = FALSE,
                                   out_dir = dir.results) {
  
  plot_list <- vector("list", length(exposures))
  names(plot_list) <- exposures
  
  for (i in seq_along(exposures)) {
    exposure <- exposures[i]
    
    rhs <- c(exposure, "age", "sex")
    
    m0 <- as.formula(
      paste0(outcome, " ~ ", paste(rhs, collapse = " + "))
    )
    
    plot_data <- data %>%
      dplyr::select(all.vars(m0), group) %>%
      tidyr::drop_na()
    
    if (nrow(plot_data) == 0) {
      plot_list[[i]] <- NULL
      next
    }
    
    m1 <- lm(m0, data = plot_data)
    
    est_row <- broom::tidy(m1) %>%
      dplyr::filter(term == exposure)
    
    fdr_row <- results_tbl %>%
      dplyr::mutate(exposure_clean = stringr::str_replace_all(exposure, "`", "")) %>%
      dplyr::filter(
        outcome == !!outcome,
        exposure %in% c(!!exposure, stringr::str_replace_all(!!exposure, "`", "")) |
          exposure_clean %in% c(!!exposure, stringr::str_replace_all(!!exposure, "`", "")),
        term %in% c(!!exposure, stringr::str_replace_all(!!exposure, "`", ""))
      ) %>%
      dplyr::slice(1)
    
    pred_df <- as.data.frame(ggeffects::ggpredict(m1, terms = exposure))
    
    subtitle_txt <- paste0(
      "Beta = ", round(est_row$estimate, 3),
      "; p = ", fmt_p(est_row$p.value),
      "; FDR p = ", fmt_p(fdr_row$p_fdr)
    )
    
    p <- ggplot(plot_data, aes(x = .data[[exposure]], y = .data[[outcome]], color = group)) +
      geom_point(alpha = 0.7, size = 2) +
      geom_line(
        data = pred_df,
        aes(x = x, y = predicted),
        inherit.aes = FALSE,
        linewidth = 1.1,
        color = "black"
      ) +
      geom_ribbon(
        data = pred_df,
        aes(x = x, ymin = conf.low, ymax = conf.high),
        inherit.aes = FALSE,
        alpha = 0.2
      ) +
      theme_bw() +
      labs(
        title = paste(get_outcome_label(outcome), "vs", get_pfas_label(exposure)),
        subtitle = subtitle_txt,
        x = get_pfas_label(exposure),
        y = get_outcome_label(outcome),
        color = "Group"
      )
    
    plot_list[[i]] <- p
  }
  
  plot_list <- plot_list[!vapply(plot_list, is.null, logical(1))]
  
  if (save_pdf && length(plot_list) > 0) {
    pdf(
      file.path(out_dir, paste0("baseline_", outcome, "_all_pfas.pdf")),
      width = 8,
      height = 6
    )
    for (nm in names(plot_list)) {
      print(plot_list[[nm]])
    }
    dev.off()
  }
  
  return(plot_list)
}


#-----------------------
#Model 2a
#group differences over time
#-------------------------
model_2a <- function(data, outcome) {
  
  rhs <- c("visit_f", "group", "visit_f:group", "age", "sex")
  
  m0 <- as.formula(
    paste0(outcome, " ~ ", paste(rhs, collapse = " + "), " + (1 | record_id)")
  )
  
  model_data <- data %>%
    dplyr::select(record_id, all.vars(m0)) %>%
    tidyr::drop_na()
  
  m1 <- lmerTest::lmer(m0, data = model_data, REML = FALSE)
  
  list(
    model = m1,
    tidy = broom.mixed::tidy(m1, effects = "fixed") %>%
      dplyr::mutate(
        analysis = "time_by_group_base_lmm",
        outcome = outcome,
        model_formula = paste(deparse(m0), collapse = " ")
      )
  )
}


#-------------------
#plotting model 2a
#-------------------
plot_model_2a <- function(data,outcomes = group_time_outcomes,save_pdf = FALSE,out_dir = dir.results) {
  
  plot_list <- vector("list", length(outcomes))
  names(plot_list) <- outcomes
  
  for (i in seq_along(outcomes)) {
    outcome <- outcomes[i]
    
    rhs <- c("visit_f", "group", "visit_f:group", "age", "sex")
    
    m0 <- as.formula(
      paste0(outcome, " ~ ", paste(rhs, collapse = " + "), " + (1 | record_id)")
    )
    
    plot_data <- data %>%
      dplyr::select(record_id, all.vars(m0)) %>%
      tidyr::drop_na()
    
    if (nrow(plot_data) == 0) {
      plot_list[[i]] <- NULL
      next
    }
    
    m1 <- lmerTest::lmer(m0, data = plot_data, REML = FALSE)
    
    pred_df <- as.data.frame(
      ggeffects::ggpredict(m1, terms = c("visit_f", "group"))
    )
    
    int_rows <- broom.mixed::tidy(m1, effects = "fixed") %>%
      dplyr::filter(stringr::str_detect(term, "visit_f.*:group|group.*:visit_f"))
    
    subtitle_text <- if (nrow(int_rows) > 0) {
      paste0(
        "Visit × Group: ",
        paste0(int_rows$term, " p = ", fmt_p(int_rows$p.value), collapse = "; ")
      )
    } else {
      NULL
    }
    
    p <- ggplot(pred_df, aes(x = x, y = predicted, color = group, group = group)) +
      geom_line(linewidth = 1.1) +
      geom_point(size = 2.5) +
      geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.08) +
      scale_color_viridis_d(option = "C") +
      theme_bw() +
      labs(
        title = paste("Unadjusted", get_outcome_label(outcome), "over time by group"),
        subtitle = paste("model 2a: age + sex + group"),
        #x = "Visit",
        #y = get_outcome_label(outcome),
        #color = "Group"
      )
    
    plot_list[[i]] <- p
  }
  
  plot_list <- plot_list[!vapply(plot_list, is.null, logical(1))]

  
  return(plot_list)
}


#-------------------------------------------------
#Model 2b--adjusting for participant baseline uACR
#-------------------------------------------------
model_2b <- function(data, outcome, adjust_acr = TRUE) {
  rhs <- c("visit_f", "group", "visit_f:group", "age", "sex") 
  
  if (adjust_acr && outcome != "acr_u" && outcome != "log_acr_u") {
    rhs <- c(rhs, "log_acr_u_bl")
  }
  
  m0 <- as.formula(
    paste0(outcome, " ~ ", paste(rhs, collapse = " + "), " + (1 | record_id)")
  )
  
  model_data <- data %>%
    dplyr::select(record_id, all.vars(m0)) %>%
    drop_na()
  
  m1 <- lmerTest::lmer(m0, data = model_data, REML = FALSE)
  
  list(
    model = m1,
    tidy = broom.mixed::tidy(m1, effects = "fixed") %>%
      mutate(
        analysis = "time_by_group",
        outcome = outcome,
        model_formula = paste(deparse(m0), collapse = " ")
      )
  )
}



#----------------------------
#plotting model 2b
#----------------------------
plot_model_2b <- function(data,
                          outcomes = group_time_outcomes,
                          adjust_acr = TRUE,
                          save_pdf = FALSE,
                          out_dir = dir.results) {
  
  plot_list <- vector("list", length(outcomes))
  names(plot_list) <- outcomes
  
  for (i in seq_along(outcomes)) {
    outcome <- outcomes[i]
    rhs <- c("visit_f", "group", "visit_f:group", "age", "sex") 
    
    if (adjust_acr && outcome != "acr_u" && outcome != "log_acr_u") {
      rhs <- c(rhs, "log_acr_u_bl")
    }
    m0 <- as.formula(
      paste0(outcome, " ~ ", paste(rhs, collapse = " + "), " + (1 | record_id)")
    )
    
    plot_data <- data %>%
      dplyr::select(record_id, all.vars(m0)) %>%
      drop_na()
    
    if (nrow(plot_data) == 0) {
      plot_list[[i]] <- NULL
      next
    }
    m1 <- lmerTest::lmer(m0, data = plot_data, REML = FALSE)
    
    pred_df <- as.data.frame(
      ggeffects::ggpredict(m1, terms = c("visit_f", "group"))
    )
    
    int_rows <- broom.mixed::tidy(m1, effects = "fixed") %>%
      filter(str_detect(term, "visit_f.*:group|group.*:visit_f"))
    
    subtitle_text <- if (nrow(int_rows) > 0) {
      paste0(
        "Visit × Group: ",
        paste0(int_rows$term, " p = ", fmt_p(int_rows$p.value), collapse = "; ")
      )
    } else {
      NULL
    }
    
    p <- ggplot(pred_df, aes(x = x, y = predicted, color = group, group = group)) +
      geom_line(linewidth = 1.1) +
      geom_point(size = 2.5) +
      geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.08) +
      scale_color_viridis_d(option = "C")+
      theme_bw() +
      labs(
        title = paste("PFAS predicted", get_outcome_label(outcome), "over time by group"),
        subtitle = paste("model 2b: age + sex + group + baseline uACR"),
       # x = "Visit",
      #  y = get_outcome_label(outcome),
      #  color = "Group"
      )
    
    plot_list[[i]] <- p
  }
  
  plot_list <- plot_list[!vapply(plot_list, is.null, logical(1))]
  
  
  return(plot_list)
}

#----------------------------------------------------
# Model 2C: group outcomes adjsuting for tanner stage
#----------------------------------------------------
model_2c <- function(data, outcome) {
  
  rhs <- c("visit_f", "group", "visit_f:group", "age", "sex", "tan_stage")
  
  m0 <- as.formula(
    paste0(outcome, " ~ ", paste(rhs, collapse = " + "), " + (1 | record_id)")
  )
  
  model_data <- data %>%
    dplyr::select(record_id, all.vars(m0)) %>%
    tidyr::drop_na()
  
  m1 <- lmerTest::lmer(m0, data = model_data, REML = FALSE)
  
  list(
    model = m1,
    tidy = broom.mixed::tidy(m1, effects = "fixed") %>%
      dplyr::mutate(
        analysis = "time_by_group_tanner",
        outcome = outcome,
        model_formula = paste(deparse(m0), collapse = " ")
      )
  )
}

#-----------------
#plotting model 2c
#------------------
plot_model_2c <- function(data,
                          outcomes = group_time_outcomes,
                          save_pdf = FALSE,
                          out_dir = dir.results) {
  
  plot_list <- vector("list", length(outcomes))
  names(plot_list) <- outcomes
  
  for (i in seq_along(outcomes)) {
    outcome <- outcomes[i]
    
    rhs <- c("visit_f", "group", "visit_f:group", "age", "sex", "tan_stage")
    
    m0 <- as.formula(
      paste0(outcome, " ~ ", paste(rhs, collapse = " + "), " + (1 | record_id)")
    )
    
    plot_data <- data %>%
      dplyr::select(record_id, all.vars(m0)) %>%
      tidyr::drop_na()
    
    if (nrow(plot_data) == 0) {
      plot_list[[i]] <- NULL
      next
    }
    
    m1 <- lmerTest::lmer(m0, data = plot_data, REML = FALSE)
    
    pred_df <- as.data.frame(
      ggeffects::ggpredict(m1, terms = c("visit_f", "group"))
    )
    
    int_rows <- broom.mixed::tidy(m1, effects = "fixed") %>%
      dplyr::filter(stringr::str_detect(term, "visit_f.*:group|group.*:visit_f"))
    
    subtitle_text <- if (nrow(int_rows) > 0) {
      paste0(
        "Visit × Group: ",
        paste0(int_rows$term, " p = ", fmt_p(int_rows$p.value), collapse = "; ")
      )
    } else {
      NULL
    }
    
    p <- ggplot(pred_df, aes(x = x, y = predicted, color = group, group = group)) +
      geom_line(linewidth = 1.1) +
      geom_point(size = 2.5) +
      geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.08) +
      scale_color_viridis_d(option = "C") +
      theme_bw() +
      labs(
        title = paste("Tanner-adjusted", get_outcome_label(outcome), "over time by group"),
        subtitle = paste("model 2a:", get_outcome_label(outcome), "~ age + sex + group + tanner stage + (1|record_id)"),
       # x = "Visit",
        #y = get_outcome_label(outcome),
        #color = "Group"
      )
    
    plot_list[[i]] <- p
  }
  
  plot_list <- plot_list[!vapply(plot_list, is.null, logical(1))]
  
  return(plot_list)
}

#------------------------------------------------------------------------
#model 2d: group differences over time adjsuint for uACR and tanner stage
#------------------------------------------------------------------------
model_2d <- function(data, outcome, adjust_acr = TRUE) {
  rhs <- c("visit_f", "group", "visit_f:group", "age", "sex", "tan_stage") 
  
  if (adjust_acr && outcome != "acr_u" && outcome != "log_acr_u") {
    rhs <- c(rhs, "log_acr_u_bl")
  }
  
  m0 <- as.formula(
    paste0(outcome, " ~ ", paste(rhs, collapse = " + "), " + (1 | record_id)")
  )
  
  model_data <- data %>%
    dplyr::select(record_id, all.vars(m0)) %>%
    drop_na()
  
  m1 <- lmerTest::lmer(m0, data = model_data, REML = FALSE)
  
  list(
    model = m1,
    tidy = broom.mixed::tidy(m1, effects = "fixed") %>%
      mutate(
        analysis = "time_by_group",
        outcome = outcome,
        model_formula = paste(deparse(m0), collapse = " ")
      )
  )
}

#-----------
#plotting 2d
#-----------
plot_model_2d <- function(data,
                          outcomes = group_time_outcomes,
                          adjust_acr = TRUE,
                          save_pdf = FALSE,
                          out_dir = dir.results) {
  
  plot_list <- vector("list", length(outcomes))
  names(plot_list) <- outcomes
  
  for (i in seq_along(outcomes)) {
    outcome <- outcomes[i]
    rhs <- c("visit_f", "group", "visit_f:group", "age", "sex", "tan_stage") 
    
    if (adjust_acr && outcome != "acr_u" && outcome != "log_acr_u") {
      rhs <- c(rhs, "log_acr_u_bl")
    }
    m0 <- as.formula(
      paste0(outcome, " ~ ", paste(rhs, collapse = " + "), " + (1 | record_id)")
    )
    
    plot_data <- data %>%
      dplyr::select(record_id, all.vars(m0)) %>%
      drop_na()
    
    if (nrow(plot_data) == 0) {
      plot_list[[i]] <- NULL
      next
    }
    m1 <- lmerTest::lmer(m0, data = plot_data, REML = FALSE)
    
    pred_df <- as.data.frame(
      ggeffects::ggpredict(m1, terms = c("visit_f", "group"))
    )
    
    int_rows <- broom.mixed::tidy(m1, effects = "fixed") %>%
      filter(str_detect(term, "visit_f.*:group|group.*:visit_f"))
    
    subtitle_text <- if (nrow(int_rows) > 0) {
      paste0(
        "Visit × Group: ",
        paste0(int_rows$term, " p = ", fmt_p(int_rows$p.value), collapse = "; ")
      )
    } else {
      NULL
    }
    
    p <- ggplot(pred_df, aes(x = x, y = predicted, color = group, group = group)) +
      geom_line(linewidth = 1.1) +
      geom_point(size = 2.5) +
      geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.08) +
      scale_color_viridis_d(option = "C")+
      theme_bw() +
      labs(
        title = paste("PFAS predicted", get_outcome_label(outcome), "over time by group"),
        subtitle = paste("model 2d: age + sex + group + baseline uACR + tanner stage"),
        #x = "Visit",
        #y = get_outcome_label(outcome),
        #color = "Group"
      )
    
    plot_list[[i]] <- p
  }
  
  plot_list <- plot_list[!vapply(plot_list, is.null, logical(1))]
  
  return(plot_list)
}


#----------------------------------------------------
#model 3a: baseline PFAS predicting longitudinal data
#age + sex + group
#taken from hailey's code in archived folder
#----------------------------------------------------
model_3a <- function(data, outcome, exposure) {

  covars <- c("age", "sex", "group")

  needed_vars <- c("record_id", outcome, exposure, "visit_num", covars)
  
  model_data <- data %>%
    dplyr::select(dplyr::all_of(needed_vars)) %>%
    tidyr::drop_na()
  
  # -----------------------------
  # Outlier removal per outcome
  # -----------------------------
  m_out <- mean(model_data[[outcome]], na.rm = TRUE)
  s_out <- sd(model_data[[outcome]], na.rm = TRUE)
  
  removed_ids <- model_data %>%
    dplyr::filter(abs(.data[[outcome]] - m_out) > 3 * s_out) %>%
    dplyr::pull(record_id) %>%
    unique()
  
  model_data <- model_data %>%
    dplyr::filter(abs(.data[[outcome]] - m_out) <= 3 * s_out)
  
  # -----------------------------
  # Sample size info
  # -----------------------------
  n_subj <- dplyr::n_distinct(model_data$record_id)
  n_obs  <- nrow(model_data)
  
  
  #formula
  rhs <- c(exposure, "visit_num", paste0(exposure, ":visit_num"), covars)
  
  m0 <- as.formula(paste0(
    outcome, " ~ ",
    paste(rhs, collapse = " + "),
    " + (1 + visit_num | record_id)"
  ))
  
  # -----------------------------
  # Fit model with fallback
  # -----------------------------
  m1 <- tryCatch(
    lmerTest::lmer(m0, data = model_data, REML = TRUE),
    error = function(e) {
      fallback_formula <- as.formula(
        paste0(
          outcome, " ~ ",
          paste(rhs, collapse = " + "),
          " + (1 | record_id)"
        )
      )
      lmerTest::lmer(fallback_formula, data = model_data, REML = TRUE)
    }
  )
  
  if (lme4::isSingular(m1, tol = 1e-4)) {
    fallback_formula <- as.formula(
      paste0(
        outcome, " ~ ",
        paste(rhs, collapse = " + "),
        " + (1 | record_id)"
      )
    )
    m1 <- lmerTest::lmer(fallback_formula, data = model_data, REML = TRUE)
  }
  
  tidy_tbl <- broom.mixed::tidy(m1, effects = "fixed") %>%
    dplyr::mutate(
      analysis = "baseline_pfas_longitudinal_base",
      outcome = outcome,
      exposure = exposure,
      n_subj = n_subj,
      n_obs = n_obs,
      removed_ids = if (length(removed_ids) > 0) paste(removed_ids, collapse = ", ") else NA_character_,
      model_formula = paste(deparse(formula(m1)), collapse = " ")
    )
  
  list(
    model = m1,
    tidy = tidy_tbl
  )
}


#-----------------------------------
#model 3b adjustin for baseline uACR
#-----------------------------------
model_3b <- function(data, outcome, exposure) {
  covars <- c("age", "sex", "group")
  
  if (!outcome %in% c("acr_u", "log_acr_u")) {
    covars <- c(covars, "log_acr_u_bl")
  }
  
  needed_vars <- c("record_id", outcome, exposure, "visit_num", covars)
  
  model_data <- data %>%
    dplyr::select(dplyr::all_of(needed_vars)) %>%
    tidyr::drop_na()
  
  # -----------------------------
  # Outlier removal per outcome
  # -----------------------------
  m_out <- mean(model_data[[outcome]], na.rm = TRUE)
  s_out <- sd(model_data[[outcome]], na.rm = TRUE)
  
  removed_ids <- model_data %>%
    dplyr::filter(abs(.data[[outcome]] - m_out) > 3 * s_out) %>%
    dplyr::pull(record_id) %>%
    unique()
  
  model_data <- model_data %>%
    dplyr::filter(abs(.data[[outcome]] - m_out) <= 3 * s_out)
  
  # -----------------------------
  # Sample size info
  # -----------------------------
  n_subj <- dplyr::n_distinct(model_data$record_id)
  n_obs  <- nrow(model_data)
  
  # -----------------------------
  # Formula
  # -----------------------------
  rhs <- c(exposure, "visit_num", paste0(exposure, ":visit_num"), covars)
  
  m0 <- as.formula(paste0(
    outcome, " ~ ",
    paste(rhs, collapse = " + "),
    " + (1 + visit_num | record_id)"
  ))
  
  # -----------------------------
  # Fit model with fallback
  # -----------------------------
  m1 <- tryCatch(
    lmerTest::lmer(m0, data = model_data, REML = TRUE),
    error = function(e) {
      fallback_formula <- as.formula(
        paste0(
          outcome, " ~ ",
          paste(rhs, collapse = " + "),
          " + (1 | record_id)"
        )
      )
      lmerTest::lmer(fallback_formula, data = model_data, REML = TRUE)
    }
  )
  
  if (lme4::isSingular(m1, tol = 1e-4)) {
    fallback_formula <- as.formula(
      paste0(
        outcome, " ~ ",
        paste(rhs, collapse = " + "),
        " + (1 | record_id)"
      )
    )
    m1 <- lmerTest::lmer(fallback_formula, data = model_data, REML = TRUE)
  }
  
  tidy_tbl <- broom.mixed::tidy(m1, effects = "fixed") %>%
    dplyr::mutate(
      analysis = "baseline_pfas_longitudinal_acr",
      outcome = outcome,
      exposure = exposure,
      n_subj = n_subj,
      n_obs = n_obs,
      removed_ids = if (length(removed_ids) > 0) paste(removed_ids, collapse = ", ") else NA_character_,
      model_formula = paste(deparse(formula(m1)), collapse = " ")
    )
  
  list(
    model = m1,
    tidy = tidy_tbl
  )
}

#----------------------------------------
#model 3c: age + sex +group + tannerstage
#baseline pfas effect on longitudinal outcomes
#----------------------------------------
model_3c <- function(data, outcome, exposure) {
  
  covars <- c("age", "sex", "group", "tan_stage")
  needed_vars <- c("record_id", outcome, exposure, "visit_num", covars)
  
  model_data <- data %>%
    dplyr::select(dplyr::all_of(needed_vars)) %>%
    tidyr::drop_na()
  
  # -----------------------------
  # Outlier removal per outcome
  # -----------------------------
  m_out <- mean(model_data[[outcome]], na.rm = TRUE)
  s_out <- sd(model_data[[outcome]], na.rm = TRUE)
  
  removed_ids <- model_data %>%
    dplyr::filter(abs(.data[[outcome]] - m_out) > 3 * s_out) %>%
    dplyr::pull(record_id) %>%
    unique()
  
  model_data <- model_data %>%
    dplyr::filter(abs(.data[[outcome]] - m_out) <= 3 * s_out)
  
  # -----------------------------
  # Sample size info
  # -----------------------------
  n_subj <- dplyr::n_distinct(model_data$record_id)
  n_obs  <- nrow(model_data)
  
  # -----------------------------
  # Formula
  # -----------------------------
  rhs <- c(exposure, "visit_num", paste0(exposure, ":visit_num"), covars)
  
  m0 <- as.formula(paste0(
    outcome, " ~ ",
    paste(rhs, collapse = " + "),
    " + (1 + visit_num | record_id)"
  ))
  
  
  # -----------------------------
  # Fit model with fallback
  # -----------------------------
  m1 <- tryCatch(
    lmerTest::lmer(m0, data = model_data, REML = TRUE),
    error = function(e) {
      fallback_formula <- as.formula(
        paste0(
          outcome, " ~ ",
          paste(rhs, collapse = " + "),
          " + (1 | record_id)"
        )
      )
      lmerTest::lmer(fallback_formula, data = model_data, REML = TRUE)
    }
  )
  
  if (lme4::isSingular(m1, tol = 1e-4)) {
    fallback_formula <- as.formula(
      paste0(
        outcome, " ~ ",
        paste(rhs, collapse = " + "),
        " + (1 | record_id)"
      )
    )
    m1 <- lmerTest::lmer(fallback_formula, data = model_data, REML = TRUE)
  }
  
  tidy_tbl <- broom.mixed::tidy(m1, effects = "fixed") %>%
    dplyr::mutate(
      analysis = "baseline_pfas_longitudinal_tanner",
      outcome = outcome,
      exposure = exposure,
      n_subj = n_subj,
      n_obs = n_obs,
      removed_ids = if (length(removed_ids) > 0) paste(removed_ids, collapse = ", ") else NA_character_,
      model_formula = paste(deparse(formula(m1)), collapse = " ")
    )
  
  list(
    model = m1,
    tidy = tidy_tbl
  )
}


#------------------------------------------------------------------
# model 3d: baseline PFAS predicting longitudinal outcomes,
# adjusted for age, sex, group, baseline uACR, and Tanner stage
#------------------------------------------------------------------
model_3d <- function(data, outcome, exposure) {
  
  covars <- c("age", "sex", "group", "tan_stage")
  
  if (!outcome %in% c("acr_u", "log_acr_u")) {
    covars <- c(covars, "log_acr_u_bl")
  }

  needed_vars <- c("record_id", outcome, exposure, "visit_num", covars)
  
  model_data <- data %>%
    dplyr::select(dplyr::all_of(needed_vars)) %>%
    tidyr::drop_na()
  
  # -----------------------------
  # Outlier removal per outcome
  # -----------------------------
  m_out <- mean(model_data[[outcome]], na.rm = TRUE)
  s_out <- sd(model_data[[outcome]], na.rm = TRUE)
  
  removed_ids <- model_data %>%
    dplyr::filter(abs(.data[[outcome]] - m_out) > 3 * s_out) %>%
    dplyr::pull(record_id) %>%
    unique()
  
  model_data <- model_data %>%
    dplyr::filter(abs(.data[[outcome]] - m_out) <= 3 * s_out)
  
  # -----------------------------
  # Sample size info
  # -----------------------------
  n_subj <- dplyr::n_distinct(model_data$record_id)
  n_obs  <- nrow(model_data)
  
  # -----------------------------
  # Formula
  # -----------------------------
  rhs <- c(exposure, "visit_num", paste0(exposure, ":visit_num"), covars)
  
  m0 <- as.formula(paste0(
    outcome, " ~ ",
    paste(rhs, collapse = " + "),
    " + (1 + visit_num | record_id)"
  ))
  

  m1 <- tryCatch(
    lmerTest::lmer(m0, data = model_data, REML = TRUE),
    error = function(e) {
      fallback_formula <- as.formula(
        paste0(
          outcome, " ~ ",
          paste(rhs, collapse = " + "),
          " + (1 | record_id)"
        )
      )
      lmerTest::lmer(fallback_formula, data = model_data, REML = TRUE)
    }
  )
  
  if (lme4::isSingular(m1, tol = 1e-4)) {
    fallback_formula <- as.formula(
      paste0(
        outcome, " ~ ",
        paste(rhs, collapse = " + "),
        " + (1 | record_id)"
      )
    )
    m1 <- lmerTest::lmer(fallback_formula, data = model_data, REML = TRUE)
  }
  
  tidy_tbl <- broom.mixed::tidy(m1, effects = "fixed") %>%
    dplyr::mutate(
      analysis = "baseline_pfas_longitudinal_acr_tanner",
      outcome = outcome,
      exposure = exposure,
      n_subj = n_subj,
      n_obs = n_obs,
      removed_ids = if (length(removed_ids) > 0) paste(removed_ids, collapse = ", ") else NA_character_,
      model_formula = paste(deparse(formula(m1)), collapse = " ")
    )
  
  list(
    model = m1,
    tidy = tidy_tbl
  )
}

