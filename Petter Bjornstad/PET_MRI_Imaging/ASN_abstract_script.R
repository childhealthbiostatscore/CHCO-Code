##### Fig 1: Renal Metabolic Efficiency: Multi-Modal Imaging


library(tidyverse)
library(ggplot2)
library(patchwork)
library(ggrepel)

your_path_to_onedrive <- "C:/Users/netio/OneDrive - UW/Laura Pyle's files - Biostatistics Core Shared Drive/"
OUTPUT_DIR <- paste0(your_path_to_onedrive, "Dylan Weissenkampen/Projects/Imaging with Shivani/")
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)


## Load and prepare data

harmonized_path <- paste0(your_path_to_onedrive, "Data Harmonization/Data Clean/harmonized_dataset.csv")

harmonized_data <- read.csv(harmonized_path, na = '')

dat <- harmonized_data %>%
  dplyr::select(-dob) %>%
  arrange(date_of_screen) %>%
  dplyr::summarise(
    across(where(negate(is.numeric)), ~ ifelse(all(is.na(.x)), NA_character_, last(na.omit(.x)))),
    across(where(is.numeric),         ~ ifelse(all(is.na(.x)), NA_real_,      mean(na.omit(.x), na.rm = TRUE))),
    .by = c(record_id, visit)
  )

# Calculate bilateral averages for PET kinetic parameters
PET_avg <- function(data) {
  tmp <- data %>% dplyr::select(lc_k2, rc_k2, lm_k2, rm_k2,
                                lc_f,  rc_f,  lm_f,  rm_f)
  avg_c_k2   <- rowMeans(tmp[, c("lc_k2", "rc_k2")], na.rm = TRUE)
  avg_m_k2   <- rowMeans(tmp[, c("lm_k2", "rm_k2")], na.rm = TRUE)
  avg_c_f    <- rowMeans(tmp[, c("lc_f",  "rc_f")],  na.rm = TRUE)
  avg_m_f    <- rowMeans(tmp[, c("lm_f",  "rm_f")],  na.rm = TRUE)
  avg_c_k2_f <- avg_c_k2 / avg_c_f
  avg_m_k2_f <- avg_m_k2 / avg_m_f
  results <- data.frame(avg_c_k2, avg_m_k2, avg_c_f, avg_m_f, avg_c_k2_f, avg_m_k2_f)
  return(results)
}

# Compute bilateral FSOC averages (medullary_fsoc_abs used in Panel A)
dat <- dat %>%
  mutate(
    medullary_fsoc_abs = case_when(
      !is.na(fsoc_l_medulla) & !is.na(fsoc_r_medulla) ~ (fsoc_l_medulla + fsoc_r_medulla) / 2,
      !is.na(fsoc_l_medulla) ~ fsoc_l_medulla,
      !is.na(fsoc_r_medulla) ~ fsoc_r_medulla,
      TRUE ~ NA_real_
    ),
    whole_kidney_fsoc_abs = case_when(
      !is.na(fsoc_l_kidney) & !is.na(fsoc_r_kidney) ~ (fsoc_l_kidney + fsoc_r_kidney) / 2,
      !is.na(fsoc_l_kidney) ~ fsoc_l_kidney,
      !is.na(fsoc_r_kidney) ~ fsoc_r_kidney,
      TRUE ~ NA_real_
    )
  )

dat_results <- dat %>%
  dplyr::select(-any_of(c("avg_c_k2", "avg_m_k2", "avg_c_f",
                          "avg_m_f",  "avg_c_k2_f", "avg_m_k2_f"))) %>%
  bind_cols(PET_avg(dat)) %>%
  filter(!is.na(avg_c_k2)) %>%
  mutate(group = if_else(record_id == "RH2-39-O", "Obese Control", group))

# Filter to 4-group comparison (exclude PKD)
dat_results <- dat_results %>%
  filter(group %in% c("Lean Control", "Obese Control", "Type 1 Diabetes", "Type 2 Diabetes"))


## Shared settings

group_levels_4 <- c("Lean Control", "Obese Control", "Type 1 Diabetes", "Type 2 Diabetes")

group_colors <- c(
  "Lean Control"    = "#87CEEB",
  "Obese Control"   = "#ADD8E6",
  "Type 1 Diabetes" = "#F0E68C",
  "Type 2 Diabetes" = "#CD5C5C"
)

theme_fig <- theme_classic(base_size = 11) +
  theme(
    axis.text.x     = element_text(angle = 45, hjust = 1, size = 10, color = "black"),
    axis.text.y     = element_text(size = 10, color = "black"),
    axis.title      = element_text(size = 11, face = "bold"),
    axis.line       = element_line(color = "black", linewidth = 0.5),
    axis.ticks      = element_line(color = "black", linewidth = 0.4),
    legend.position = "none",
    plot.tag        = element_text(size = 13, face = "bold")
  )


## Build linked dataset and discordant phenotype (Panel E)

create_linked_k2_fsoc_dataset <- function(dat) {
  
  possible_demo_vars <- c("sex", "age", "bmi", "hba1c",
                          "eGFR_CKD_epi", "gfr_fas_cr", "gfr_raw_plasma", "gfr_bsa_plasma",
                          "acr_u", "uacr", "diabetes_duration", "dm_duration",
                          "sbp", "dbp", "weight", "height",
                          "dexa_fat_kg", "dexa_trunk_kg", "dexa_lean_kg",
                          "dexa_body_fat", "creatinine", "cystatin_c")
  available_demo_vars <- possible_demo_vars[possible_demo_vars %in% names(dat)]
  
  k2_vars   <- unique(c("record_id", "mrn", "group", available_demo_vars, "avg_c_k2", "avg_m_k2"))
  k2_vars   <- k2_vars[k2_vars %in% names(dat)]
  fsoc_vars <- c("record_id", "mrn", "group",
                 "fsoc_l_medulla", "fsoc_r_medulla", "fsoc_l_kidney", "fsoc_r_kidney")
  fsoc_vars <- fsoc_vars[fsoc_vars %in% names(dat)]
  
  k2_data <- dat %>%
    filter(!is.na(avg_c_k2) | !is.na(avg_m_k2), !is.na(mrn), mrn != "") %>%
    select(all_of(k2_vars)) %>%
    rename(k2_record_id = record_id, k2_group = group) %>%
    distinct(mrn, .keep_all = TRUE)
  
  fsoc_data <- dat %>%
    filter(!is.na(fsoc_l_medulla) | !is.na(fsoc_r_medulla), !is.na(mrn), mrn != "") %>%
    select(all_of(fsoc_vars)) %>%
    rename(fsoc_record_id = record_id, fsoc_group = group) %>%
    distinct(mrn, .keep_all = TRUE)
  
  linked_data <- inner_join(k2_data, fsoc_data, by = "mrn")
  cat("Linked subjects (K2 + FSOC):", nrow(linked_data), "\n")
  if (nrow(linked_data) == 0) return(NULL)
  
  linked_data <- linked_data %>%
    mutate(
      fsoc_medulla = case_when(
        !is.na(fsoc_l_medulla) & !is.na(fsoc_r_medulla) ~ (fsoc_l_medulla + fsoc_r_medulla) / 2,
        !is.na(fsoc_l_medulla) ~ fsoc_l_medulla,
        !is.na(fsoc_r_medulla) ~ fsoc_r_medulla,
        TRUE ~ NA_real_
      ),
      group = k2_group
    )
  
  return(linked_data)
}

run_lean_control_phenotype <- function(linked_data, k2_var = "avg_c_k2",
                                       fsoc_var = "fsoc_medulla", n_sd = 0.5) {
  
  dat_valid     <- linked_data %>% filter(!is.na(.data[[k2_var]]), !is.na(.data[[fsoc_var]]))
  lean_controls <- dat_valid %>% filter(grepl("Lean", group, ignore.case = TRUE))
  
  if (nrow(lean_controls) < 2) {
    lean_controls <- dat_valid %>% filter(grepl("Control", group, ignore.case = TRUE))
  }
  
  lc_k2_mean   <- mean(lean_controls[[k2_var]],   na.rm = TRUE)
  lc_k2_sd     <- sd(lean_controls[[k2_var]],     na.rm = TRUE)
  lc_fsoc_mean <- mean(lean_controls[[fsoc_var]], na.rm = TRUE)
  lc_fsoc_sd   <- sd(lean_controls[[fsoc_var]],   na.rm = TRUE)
  
  k2_cutoff   <- lc_k2_mean   + n_sd * lc_k2_sd
  fsoc_cutoff <- lc_fsoc_mean - n_sd * lc_fsoc_sd
  
  cat(sprintf("K2 cutoff: %.4f  |  FSOC cutoff: %.2f\n", k2_cutoff, fsoc_cutoff))
  
  dat_valid <- dat_valid %>%
    mutate(
      k2_level    = ifelse(.data[[k2_var]]   > k2_cutoff,   "High K2",    "Normal K2"),
      fsoc_level  = ifelse(.data[[fsoc_var]] < fsoc_cutoff, "Low FSOC",   "Normal FSOC"),
      phenotype_4 = case_when(
        k2_level == "High K2"   & fsoc_level == "Low FSOC"    ~ "High K2 / Low FSOC",
        k2_level == "High K2"   & fsoc_level == "Normal FSOC" ~ "High K2 / Normal FSOC",
        k2_level == "Normal K2" & fsoc_level == "Low FSOC"    ~ "Normal K2 / Low FSOC",
        k2_level == "Normal K2" & fsoc_level == "Normal FSOC" ~ "Normal K2 / Normal FSOC"
      ),
      discordant  = ifelse(k2_level == "High K2" & fsoc_level == "Low FSOC",
                           "Discordant", "Concordant")
    )
  
  cat("Discordant:", sum(dat_valid$discordant == "Discordant"), "/", nrow(dat_valid), "\n")
  
  return(list(
    data               = dat_valid,
    cutoffs            = list(k2 = k2_cutoff, fsoc = fsoc_cutoff),
    lean_control_stats = list(k2_mean   = lc_k2_mean,   k2_sd   = lc_k2_sd,
                              fsoc_mean = lc_fsoc_mean, fsoc_sd = lc_fsoc_sd),
    n_sd               = n_sd
  ))
}

linked_data <- create_linked_k2_fsoc_dataset(dat)
pheno_05sd  <- run_lean_control_phenotype(linked_data, n_sd = 0.5)


## IQR summaries

summarise_iqr <- function(data, var, label, grp_levels = group_levels_4, digits = 3) {
  out <- data %>%
    filter(!is.na(.data[[var]]), group %in% grp_levels) %>%
    mutate(Group = factor(group, levels = grp_levels)) %>%
    group_by(Group) %>%
    summarise(
      n      = n(),
      Median = round(median(.data[[var]], na.rm = TRUE), digits),
      Q1     = round(quantile(.data[[var]], 0.25, na.rm = TRUE), digits),
      Q3     = round(quantile(.data[[var]], 0.75, na.rm = TRUE), digits),
      IQR    = round(IQR(.data[[var]], na.rm = TRUE), digits),
      Mean   = round(mean(.data[[var]], na.rm = TRUE), digits),
      SD     = round(sd(.data[[var]], na.rm = TRUE), digits),
      .groups = "drop"
    )
  cat("\n======================================================\n")
  cat(" IQR SUMMARY:", label, "\n")
  cat("======================================================\n")
  print(out, n = Inf)
  invisible(out)
}

summarise_iqr(dat %>% filter(visit == "baseline", medullary_fsoc_abs >= 0),
              "medullary_fsoc_abs", "Medullary FSOC — Panel A")
summarise_iqr(dat_results, "avg_c_k2",   "Cortical K2 — Panel B")
summarise_iqr(dat_results, "avg_c_k2_f", "Cortical K2/F — Panel C")
summarise_iqr(pheno_05sd$data, "avg_c_k2", "Cortical K2 — dual-imaging subset (Panel E)")


## Boxplot function (Panels A, B, C)
# Kruskal-Wallis overall + pairwise Wilcoxon with BH correction

plot_group <- function(data, value_var, ylab, tag,
                       grp_levels = group_levels_4) {
  
  df <- data %>%
    filter(!is.na(.data[[value_var]]), group %in% grp_levels) %>%
    mutate(grp = factor(group, levels = grp_levels))
  
  ymax <- max(df[[value_var]], na.rm = TRUE)
  yr   <- diff(range(df[[value_var]], na.rm = TRUE))
  
  kw_p   <- kruskal.test(df[[value_var]] ~ df$grp)$p.value
  kw_lab <- ifelse(kw_p < 0.0001, "KW p<0.0001",
                   ifelse(kw_p < 0.001,  "KW p<0.001",
                          paste0("KW p=", formatC(kw_p, format = "f", digits = 3))))
  
  pw <- pairwise.wilcox.test(df[[value_var]], df$grp, p.adjust.method = "BH")$p.value
  
  fmt_p <- function(pv) {
    if (is.na(pv))   return(NA_character_)
    if (pv < 0.0001) return("p<0.0001")
    if (pv < 0.001)  return("p<0.001")
    paste0("p=", formatC(pv, format = "f", digits = 3))
  }
  
  get_pw_p <- function(g1, g2) {
    r <- tryCatch(pw[g2, g1], error = function(e) NA_real_)
    if (!is.na(r)) return(fmt_p(r))
    r <- tryCatch(pw[g1, g2], error = function(e) NA_real_)
    if (!is.na(r)) return(fmt_p(r))
    NA_character_
  }
  
  p_lc_oc  <- get_pw_p("Lean Control",    "Obese Control")
  p_lc_t1d <- get_pw_p("Lean Control",    "Type 1 Diabetes")
  p_lc_t2d <- get_pw_p("Lean Control",    "Type 2 Diabetes")
  p_oc_t1d <- get_pw_p("Obese Control",   "Type 1 Diabetes")
  p_oc_t2d <- get_pw_p("Obese Control",   "Type 2 Diabetes")
  p_t1_t2d <- get_pw_p("Type 1 Diabetes", "Type 2 Diabetes")
  
  b <- ymax + yr * c(0.08, 0.20, 0.32, 0.44, 0.56, 0.70)
  
  add_bracket <- function(p, x1, x2, y, label) {
    if (is.na(label)) return(p)
    p +
      annotate("segment", x = x1, xend = x2,
               y = y, yend = y, linewidth = 0.4, color = "black") +
      annotate("segment", x = x1, xend = x1,
               y = y - yr * 0.01, yend = y, linewidth = 0.4, color = "black") +
      annotate("segment", x = x2, xend = x2,
               y = y - yr * 0.01, yend = y, linewidth = 0.4, color = "black") +
      annotate("text", x = (x1 + x2) / 2, y = y + yr * 0.025,
               label = label, size = 2.8, fontface = "italic")
  }
  
  p <- ggplot(df, aes(x = grp, y = .data[[value_var]], fill = grp)) +
    geom_boxplot(alpha = 0.75, width = 0.55, linewidth = 0.5, outlier.shape = NA) +
    geom_jitter(aes(color = grp), width = 0.15, size = 1.8, alpha = 0.6, shape = 16) +
    scale_fill_manual(values  = group_colors) +
    scale_color_manual(values = group_colors) +
    annotate("text", x = 0.6, y = b[6] + yr * 0.05,
             label = kw_lab, size = 2.8, fontface = "italic", hjust = 0) +
    labs(x = NULL, y = ylab, tag = tag) +
    theme_fig +
    coord_cartesian(ylim = c(NA, b[6] + yr * 0.12))
  
  p <- add_bracket(p, 1, 2, b[1], p_lc_oc)
  p <- add_bracket(p, 2, 3, b[2], p_oc_t1d)
  p <- add_bracket(p, 3, 4, b[3], p_t1_t2d)
  p <- add_bracket(p, 2, 4, b[4], p_oc_t2d)
  p <- add_bracket(p, 1, 3, b[5], p_lc_t1d)
  p <- add_bracket(p, 1, 4, b[6], p_lc_t2d)
  
  p
}


## Panel A

pA <- plot_group(
  data      = dat %>% filter(visit == "baseline", medullary_fsoc_abs >= 0),
  value_var = "medullary_fsoc_abs",
  ylab      = expression(bold("FSOC Value (s"^{-1}*")")),
  tag       = "A"
)


## Panel B

pB <- plot_group(
  data      = dat_results,
  value_var = "avg_c_k2",
  ylab      = expression(bold("Cortical K"[2]*" (s"^{-1}*")")),
  tag       = "B"
)


## Panel C

pC <- plot_group(
  data      = dat_results,
  value_var = "avg_c_k2_f",
  ylab      = expression(bold("Cortical K"[2]*"/F Ratio")),
  tag       = "C"
)


## Panel D
# Load immune proteomics results (saved by FSOC proteomics script)
immune_results <- read.csv("C:/Users/netio/OneDrive - UW/Laura Pyle's files - Biostatistics Core Shared Drive/Dylan Weissenkampen/Projects/Imaging with Shivani/immune_proteins_results.csv")

sig_immune <- immune_results %>% filter(p_fdr < 0.10)

pD <- ggplot(immune_results, aes(x = log2_fc, y = -log10(p_value))) +
  geom_point(
    data  = immune_results %>% filter(p_fdr >= 0.10),
    color = "gray60", size = 1.5, alpha = 0.5
  ) +
  geom_point(
    data  = sig_immune,
    color = "#E41A1C", size = 2.5, alpha = 0.8
  ) +
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed", color = "blue", alpha = 0.5, linewidth = 0.5) +
  geom_vline(xintercept = 0,
             linetype = "dashed", alpha = 0.4, linewidth = 0.5) +
  geom_text_repel(
    data          = sig_immune,
    aes(label     = EntrezGeneSymbol),
    size          = 2.5,
    fontface      = "bold",
    max.overlaps  = 30,
    segment.size  = 0.3,
    segment.color = "gray40",
    box.padding   = 0.4,
    point.padding = 0.2
  ) +
  labs(
    x        = expression("Log"[2]*" Fold Change (Impaired / Normal FSOC)"),
    y        = expression("-log"[10]*"(p-value)"),
    tag      = "D",
    title    = "Immune-Related Proteins: Impaired vs Normal Medullary FSOC",
    subtitle = paste0(
      nrow(immune_results), " immune proteins, ",
      sum(immune_results$p_value < 0.05), " p<0.05, ",
      sum(immune_results$p_fdr   < 0.10), " FDR<0.10"
    )
  ) +
  theme_classic(base_size = 11) +
  theme(
    legend.position  = "top",
    panel.grid.minor = element_blank(),
    plot.tag         = element_text(size = 13, face = "bold"),
    plot.title       = element_text(size = 9,  face = "bold"),
    plot.subtitle    = element_text(size = 8,  color = "#475569")
  )


## Panel E

dat_scat <- pheno_05sd$data
k2_cut   <- pheno_05sd$cutoffs$k2
fsoc_cut <- pheno_05sd$cutoffs$fsoc
lc_stats <- pheno_05sd$lean_control_stats

n_disc  <- sum(dat_scat$discordant == "Discordant", na.rm = TRUE)
n_total <- nrow(dat_scat)

x_range <- range(dat_scat$fsoc_medulla, na.rm = TRUE)
y_range <- range(dat_scat$avg_c_k2,     na.rm = TRUE)
x_pad   <- diff(x_range) * 0.10
y_pad   <- diff(y_range) * 0.10

pE <- ggplot(dat_scat, aes(x = fsoc_medulla, y = avg_c_k2)) +
  annotate("rect",
           xmin = -Inf,    xmax = fsoc_cut,
           ymin = k2_cut,  ymax = Inf,
           fill = "#FF6B6B", alpha = 0.25) +
  annotate("rect",
           xmin = fsoc_cut, xmax = Inf,
           ymin = -Inf,     ymax = k2_cut,
           fill = "#4ECDC4", alpha = 0.25) +
  annotate("rect",
           xmin = -Inf,    xmax = fsoc_cut,
           ymin = -Inf,    ymax = k2_cut,
           fill = "#FFE66D", alpha = 0.20) +
  annotate("rect",
           xmin = fsoc_cut, xmax = Inf,
           ymin = k2_cut,   ymax = Inf,
           fill = "#95E1D3", alpha = 0.20) +
  geom_vline(xintercept = fsoc_cut,
             linetype = "dashed", color = "#E63946", linewidth = 1.0) +
  geom_hline(yintercept = k2_cut,
             linetype = "dashed", color = "#E63946", linewidth = 1.0) +
  geom_vline(xintercept = lc_stats$fsoc_mean,
             linetype = "dotted", color = "blue", linewidth = 0.7) +
  geom_hline(yintercept = lc_stats$k2_mean,
             linetype = "dotted", color = "blue", linewidth = 0.7) +
  geom_point(
    aes(fill = group, shape = discordant),
    size = 5, alpha = 0.85, color = "black", stroke = 0.5
  ) +
  scale_fill_manual(values = group_colors, name = "Group") +
  scale_shape_manual(values = c("Concordant" = 21, "Discordant" = 24), name = "Phenotype") +
  labs(
    x     = expression(bold("Medullary FSOC (s"^{-1}*")")),
    y     = expression(bold("Cortical K"[2]*" \u2013 TCA Metabolism (s"^{-1}*")")),
    tag   = "E",
    title = "K2 vs FSOC Phenotype Classification"
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position  = "right",
    legend.text      = element_text(size = 8.5),
    legend.title     = element_text(size = 9, face = "bold"),
    plot.tag         = element_text(size = 13, face = "bold"),
    plot.title       = element_text(size = 10, face = "bold"),
    plot.subtitle    = element_text(size = 7,  color = "#475569"),
    axis.title       = element_text(size = 10, face = "bold")
  ) +
  guides(
    fill  = guide_legend(override.aes = list(shape = 21, size = 3.5)),
    shape = guide_legend(override.aes = list(fill  = "gray50", size = 3.5))
  )


## Combine

fig1 <- (pA | pB | pC) / (pD | pE) +
  plot_layout(
    heights = c(1, 1.1),
    guides  = "keep"
  )

print(fig1)


## Outputs

ggsave(file.path(OUTPUT_DIR, "Figure1.pdf"),
       fig1, width = 14, height = 11, dpi = 300)

ggsave(file.path(OUTPUT_DIR, "Figure1.tiff"),
       fig1, width = 14, height = 11, dpi = 300, compression = "lzw")

ggsave(file.path(OUTPUT_DIR, "Figure1.png"),
       fig1, width = 14, height = 11, dpi = 300)

cat("\nFigure 1 saved to:", OUTPUT_DIR, "\n")
