#Functions
render.geometric <- function(x, ...) {
  gm <- exp(mean(log(x[x > 0]), na.rm = TRUE))  # geometric mean
  gsd <- exp(sd(log(x[x > 0]), na.rm = TRUE))    # geometric SD
  c("", 
    " " = sprintf("%.2f (%.2f)", gm, gsd))
}

#collapsing metadata function
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


#fucntions for pvalues and labels
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


#imputation function
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

#PFAS spearman correlations
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

#baseline clinical data comparisons--no PFAS
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

baseline_group_outcome_results <- lapply(
  baseline_model_outcomes,
  function(outcome) {
    run_baseline_group_model(
      data = dat_baseline,
      outcome = outcome,
      covars = c("age", "sex")
    )
  }
)

names(baseline_group_outcome_results) <- baseline_model_outcomes

baseline_group_outcome_omnibus <- bind_rows(
  lapply(baseline_group_outcome_results, `[[`, "omnibus")
) %>%
  dplyr::mutate(
    outcome_label = sapply(outcome, get_outcome_label),
    p_fdr = p.adjust(p.value, method = "fdr"),
    p_fmt = fmt_p(p.value),
    p_fdr_fmt = fmt_p(p_fdr)
  )

baseline_group_outcome_pairwise <- bind_rows(
  lapply(baseline_group_outcome_results, `[[`, "pairwise")
) %>%
  dplyr::mutate(
    outcome_label = sapply(outcome, get_outcome_label),
    p_fmt = fmt_p(p.value)
  )

#panel boxplots of clinical variables
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

#baseline PFAS across groups
baseline_group_pfas_results <- lapply(
  pfas_group_vars,
  function(outcome) {
    run_baseline_group_model(
      data = dat_baseline,
      outcome = outcome,
      covars = c("age", "sex")
    )
  }
)

#plotting baseline PFAS across groups--individual plots
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

#Model 1: Baseline Cross Section PFAS Analysis
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