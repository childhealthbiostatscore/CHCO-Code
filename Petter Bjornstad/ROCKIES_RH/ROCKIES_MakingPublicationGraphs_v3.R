########################################################################
# ROCKIES Publication - Complete Figure Generation Script
# Figures 1-7 for Cell Metabolism Submission
########################################################################

library(tidyverse)
library(ggpubr)
library(patchwork)
library(ggbeeswarm)
library(data.table)

base_path <- 'C:/Users/netio/Documents/UofW/Rockies/publication_figures/'

########################################################################
# GLOBAL THEME & COLORS
########################################################################

if (.Platform$OS.type == "windows") {
  tryCatch({
    library(extrafont)
    loadfonts(device = "pdf", quiet = TRUE)
    loadfonts(device = "win", quiet = TRUE)
  }, error = function(e) message("extrafont loading skipped: ", e$message))
}

font_family <- "sans"

theme_rockies <- theme_classic(base_size = 11) +
  theme(
    text = element_text(family = font_family),
    axis.title = element_text(size = 10, face = "bold"),
    axis.text = element_text(size = 9, color = "black"),
    plot.title = element_text(size = 11, face = "bold", hjust = 0.5),
    plot.tag = element_text(size = 14, face = "bold"),
    strip.text = element_text(size = 10, face = "bold"),
    strip.background = element_blank(),
    legend.position = "none",
    plot.margin = margin(5, 10, 5, 10)
  )

cols <- list(
  placebo = "#4A90D9", ertu = "#E74C3C", paired_line = "gray60",
  t2d = "#E74C3C", control = "#4A90D9", obese_ctrl = "#F39C12",
  corr_point = "#2C3E50", delta = "#8E44AD",
  gbm_yes = "#E74C3C", gbm_no = "#4A90D9",
  arterio_yes = "#E74C3C", arterio_no = "#4A90D9",
  up = "#E74C3C", down = "#4A90D9", ns = "gray70"
)

save_fig <- function(plot, name, width, height, dpi = 300) {
  ggsave(paste0(base_path, name, ".png"), plot,
         width = width, height = height, dpi = dpi)
  tryCatch({
    ggsave(paste0(base_path, name, ".pdf"), plot,
           width = width, height = height, dpi = dpi, device = cairo_pdf)
  }, error = function(e) {
    message("cairo_pdf failed, trying default pdf...")
    tryCatch({
      ggsave(paste0(base_path, name, ".pdf"), plot,
             width = width, height = height, dpi = dpi)
    }, error = function(e2) {
      message("PDF failed: ", e2$message, " - PNG saved.")
    })
  })
}

########################################################################
# DATA LOADING
########################################################################

# =====================================================================
# ROCKIES trial data — Dutch CSV: semicolon sep, comma decimal
# =====================================================================
rockies_wide <- read.csv2(
  "C:/Users/netio/Downloads/PET data ROCKIES variables(Blad1).csv",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

names(rockies_wide) <- make.names(trimws(names(rockies_wide)))

# Convert all columns: try numeric (replacing comma decimals)
rockies_wide <- rockies_wide %>%
  mutate(across(everything(), ~ {
    conv <- suppressWarnings(as.numeric(gsub(",", ".", as.character(.x))))
    if (sum(!is.na(conv)) >= sum(!is.na(.x)) * 0.3 & sum(!is.na(conv)) > 0) conv else .x
  }))

cat("=== ROCKIES columns ===\n")
print(names(rockies_wide))
cat("\nRows:", nrow(rockies_wide), "\n")
cat("Group 0 (placebo-only):", sum(rockies_wide$Group == 0, na.rm = TRUE), "\n")
cat("Group 1 (crossover):", sum(rockies_wide$Group == 1, na.rm = TRUE), "\n\n")

# =====================================================================
# Reshape to long: one row per participant per treatment
# Group==1 = crossover participants with both arms
# =====================================================================
rockies_crossover <- rockies_wide %>% filter(Group == 1)

rockies_long <- bind_rows(
  # --- Placebo ---
  rockies_crossover %>% transmute(
    id            = Participant,
    treatment     = "Placebo",
    cortical_k2   = K2_cortex_mean_Sherbrook_PLB,
    medullary_k2  = K2_medulla_mean_Sherbrook_PLB,
    cortical_k2f  = K2_F_cortex_mean_Sherbrook_PLB,
    medullary_k2f = K2_F_medulla_mean_Sherbrook_PLB,
    homa_ir       = HOMAIR_Placebo,
    matsuda       = Matsuda_Placebo,
    ogis          = OGIS_placebo,
    gfr           = GFR_Placebo_ml_min,
    erpf          = ERPF_Placebo_ml_min,
    ff            = FF_Placebo,
    tna_sodium    = TNA_sodium_plb_mmolpermin,
    sbp           = SBP_Placebo,
    dbp           = DBP_Placebo,
    weight        = Bodyweight_Placebo,
    bmi           = BMI_Placebo,
    uacr          = Ur_albuminuria_24hr_placebo,
    ur_sodium     = Ur_sodium_24hr_placebo_mmol,
    ur_glucose    = Ur_glucose_24hr_placebo_mmol
  ),
  # --- Ertugliflozin ---
  rockies_crossover %>% transmute(
    id            = Participant,
    treatment     = "Ertugliflozin",
    cortical_k2   = K2_cortex_mean_Sherbrook_SGLT2,
    medullary_k2  = K2_medulla_mean_Sherbrook_SGLT2i,
    cortical_k2f  = K2_F_cortex_mean_Sherbrook_SGLT2,
    medullary_k2f = K2_F_medulla_mean_Sherbrook_SGLT2,
    homa_ir       = HOMAIR_SGLT2i,
    matsuda       = Matsuda_SGLT2i,
    ogis          = OGIS_SGLT2i,
    gfr           = GFR_SGLT2i_ml_min,
    erpf          = ERPF_SGLT2i_ml_min,
    ff            = FF_SGLT2,
    tna_sodium    = TNA_sodium_SLGT2i_mmolpermin,
    sbp           = SBP_verum,
    dbp           = DBP_Verum,
    weight        = Bodyweight_Verum,
    bmi           = BMI_Verum,
    uacr          = Ur_albuminuria_24hr_SGLT2i,
    ur_sodium     = Ur_sodium_24hr_sglt2i_mmol,
    ur_glucose    = Ur_glucose_24hr_sglt2i_mmol
  )
)

rockies_long$treatment <- factor(rockies_long$treatment,
                                 levels = c("Placebo", "Ertugliflozin"))

cat("ROCKIES long:", nrow(rockies_long), "rows,",
    n_distinct(rockies_long$id), "participants\n\n")

rockies_plac <- rockies_long %>% filter(treatment == "Placebo")

# =====================================================================
# Harmonized data (Figures 2-5)
# =====================================================================
harmonized_data <- read.csv(
  "C:/Users/netio/OneDrive - UW/Laura Pyle's files - Biostatistics Core Shared Drive/Data Harmonization/Data Clean/harmonized_dataset.csv",
  na = ''
)

dat <- harmonized_data %>%
  dplyr::select(-dob) %>%
  arrange(date_of_screen) %>%
  dplyr::summarise(
    across(where(negate(is.numeric)), ~ ifelse(all(is.na(.x)), NA_character_, last(na.omit(.x)))),
    across(where(is.numeric), ~ ifelse(all(is.na(.x)), NA_real_, mean(na.omit(.x), na.rm = TRUE))),
    .by = c(record_id, visit)
  )

dat_clinical <- fread("C:/Users/netio/Downloads/UACR_Allparticipants_forGBM.csv")
dat_clinical <- dat_clinical[-str_which(dat_clinical$record_id, '-O')]

########################################################################
# PET AVERAGES — computed once from harmonized data
########################################################################

PET_avg <- function(data) {
  data.frame(
    avg_c_k2   = rowMeans(data[, c("lc_k2", "rc_k2")], na.rm = TRUE),
    avg_m_k2   = rowMeans(data[, c("lm_k2", "rm_k2")], na.rm = TRUE),
    avg_c_f    = rowMeans(data[, c("lc_f", "rc_f")],   na.rm = TRUE),
    avg_m_f    = rowMeans(data[, c("lm_f", "rm_f")],   na.rm = TRUE),
    avg_c_k2_f = rowMeans(data[, c("lc_k2", "rc_k2")], na.rm = TRUE) /
      rowMeans(data[, c("lc_f", "rc_f")],   na.rm = TRUE),
    avg_m_k2_f = rowMeans(data[, c("lm_k2", "rm_k2")], na.rm = TRUE) /
      rowMeans(data[, c("lm_f", "rm_f")],   na.rm = TRUE)
  )
}

pet_col_names <- c('avg_c_k2','avg_m_k2','avg_c_f','avg_m_f','avg_c_k2_f','avg_m_k2_f')
dat_clean <- dat[, !(names(dat) %in% pet_col_names), drop = FALSE]
dat_with_pet <- cbind(dat_clean, PET_avg(dat_clean))

dat_pet_slim <- dat_with_pet[!is.na(dat_with_pet$avg_c_k2) & is.finite(dat_with_pet$avg_c_k2),
                             c("record_id", "avg_c_k2", "avg_c_f", "avg_c_k2_f")]

cat("dat_with_pet:", nrow(dat_with_pet), "rows\n")
cat("dat_pet_slim:", nrow(dat_pet_slim), "rows with PET data\n\n")

########################################################################
# REUSABLE PLOT FUNCTIONS
########################################################################

plot_paired <- function(data, id_col, trt_col, val_col, ylab, tag) {
  df <- data.frame(id = data[[id_col]], trt = data[[trt_col]], val = data[[val_col]])
  df <- df[!is.na(df$val), ]
  df$trt <- factor(df$trt, levels = c("Placebo", "Ertugliflozin"))
  
  summ <- df %>% group_by(trt) %>%
    summarise(m = mean(val), se = sd(val)/sqrt(n()), .groups = "drop")
  
  wide <- df %>% pivot_wider(names_from = trt, values_from = val)
  pv <- wilcox.test(wide$Placebo, wide$Ertugliflozin, paired = TRUE)$p.value
  plab <- ifelse(pv < 0.001, paste0("p=", formatC(pv, format = "e", digits = 1)),
                 paste0("p=", formatC(pv, format = "f", digits = 3)))
  
  ymax <- max(df$val, na.rm = TRUE); ymin <- min(df$val, na.rm = TRUE)
  yr <- ymax - ymin
  
  ggplot(df, aes(trt, val)) +
    geom_line(aes(group = id), color = cols$paired_line, alpha = 0.5, linewidth = 0.4) +
    geom_point(aes(color = trt), size = 2.5, alpha = 0.7) +
    geom_errorbar(data = summ, aes(trt, m, ymin = m - se, ymax = m + se),
                  width = 0.15, linewidth = 0.8, inherit.aes = FALSE) +
    geom_point(data = summ, aes(trt, m), size = 4, shape = 18, inherit.aes = FALSE) +
    annotate("segment", x = 1, xend = 2,
             y = ymax + yr * 0.05, yend = ymax + yr * 0.05, linewidth = 0.4) +
    annotate("text", x = 1.5, y = ymax + yr * 0.1,
             label = plab, size = 3.5, fontface = "italic") +
    scale_color_manual(values = c(cols$placebo, cols$ertu)) +
    coord_cartesian(ylim = c(ymin - yr * 0.05, ymax + yr * 0.15)) +
    labs(x = NULL, y = ylab, tag = tag) +
    theme_rockies
}

plot_corr <- function(data, xvar, yvar, xlab, ylab, tag,
                      point_col = cols$corr_point, method = "spearman") {
  df <- data.frame(x = data[[xvar]], y = data[[yvar]])
  df <- df[!is.na(df$x) & !is.na(df$y), ]
  
  ct <- cor.test(df$x, df$y, method = method)
  rho_label <- ifelse(method == "spearman", "rho", "r")
  lab <- paste0(rho_label, "=", formatC(ct$estimate, format = "f", digits = 2),
                ", p=", ifelse(ct$p.value < 0.001, "<0.001",
                               formatC(ct$p.value, format = "f", digits = 3)))
  
  ggplot(df, aes(x, y)) +
    geom_point(size = 2.5, alpha = 0.7, color = point_col) +
    geom_smooth(method = "lm", se = TRUE, color = "black",
                linewidth = 0.8, fill = "gray85", formula = y ~ x) +
    annotate("text", x = min(df$x) + diff(range(df$x)) * 0.02,
             y = max(df$y) - diff(range(df$y)) * 0.02,
             label = lab, hjust = 0, size = 3, fontface = "italic") +
    labs(x = xlab, y = ylab, tag = tag) +
    theme_rockies
}

plot_delta_corr <- function(data, id_col, trt_col, xvar, yvar,
                            xlab, ylab, tag, method = "spearman") {
  df <- data.frame(id = data[[id_col]], trt = data[[trt_col]],
                   x = data[[xvar]], y = data[[yvar]])
  df <- df[!is.na(df$x) & !is.na(df$y), ]
  
  wx <- df %>% dplyr::select(id, trt, x) %>%
    pivot_wider(names_from = trt, values_from = x)
  wy <- df %>% dplyr::select(id, trt, y) %>%
    pivot_wider(names_from = trt, values_from = y)
  
  deltas <- data.frame(
    dx = wx$Ertugliflozin - wx$Placebo,
    dy = wy$Ertugliflozin - wy$Placebo
  )
  deltas <- deltas[!is.na(deltas$dx) & !is.na(deltas$dy), ]
  
  ct <- cor.test(deltas$dx, deltas$dy, method = method)
  rho_label <- ifelse(method == "spearman", "rho", "r")
  lab <- paste0(rho_label, "=", formatC(ct$estimate, format = "f", digits = 2),
                ", p=", ifelse(ct$p.value < 0.001, "<0.001",
                               formatC(ct$p.value, format = "f", digits = 3)))
  
  ggplot(deltas, aes(dx, dy)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    geom_point(size = 2.5, alpha = 0.7, color = cols$delta) +
    geom_smooth(method = "lm", se = TRUE, color = "black",
                linewidth = 0.8, fill = "gray85", formula = y ~ x) +
    annotate("text", x = min(deltas$dx) + diff(range(deltas$dx)) * 0.02,
             y = max(deltas$dy) - diff(range(deltas$dy)) * 0.02,
             label = lab, hjust = 0, size = 3, fontface = "italic") +
    labs(x = paste0("\u0394 ", xlab), y = paste0("\u0394 ", ylab), tag = tag) +
    theme_rockies
}

plot_two_group <- function(data, group_var, value_var,
                           group_colors, ylab, tag, test = "wilcox.test") {
  df <- data.frame(grp = data[[group_var]], val = data[[value_var]])
  df <- df[!is.na(df$val) & !is.na(df$grp), ]
  df$grp <- factor(df$grp, levels = names(group_colors))
  
  pv <- if (test == "wilcox.test") {
    wilcox.test(val ~ grp, data = df)$p.value
  } else { t.test(val ~ grp, data = df)$p.value }
  
  plab <- ifelse(pv < 0.0001, "p<0.0001",
                 ifelse(pv < 0.001, "p<0.001",
                        paste0("p=", formatC(pv, format = "f", digits = 3))))
  
  ymax <- max(df$val, na.rm = TRUE)
  yr <- diff(range(df$val, na.rm = TRUE))
  
  ggplot(df, aes(grp, val, fill = grp)) +
    geom_boxplot(alpha = 0.3, outlier.shape = NA, width = 0.5) +
    geom_jitter(aes(color = grp), width = 0.12, size = 2, alpha = 0.7) +
    annotate("segment", x = 1, xend = 2,
             y = ymax + yr * 0.05, yend = ymax + yr * 0.05, linewidth = 0.4) +
    annotate("text", x = 1.5, y = ymax + yr * 0.1,
             label = plab, size = 3.5, fontface = "italic") +
    scale_fill_manual(values = group_colors) +
    scale_color_manual(values = group_colors) +
    labs(x = NULL, y = ylab, tag = tag) +
    theme_rockies
}

########################################################################
# =====================================================================
#  FIGURE 1: ROCKIES TRIAL — k2 and k2/F focused
# =====================================================================
########################################################################

# =====================================================================
# Panel A: Study Design Diagram (generated in R)
# =====================================================================

fig1a <- ggplot() +
  xlim(0, 16) + ylim(0, 10) +
  annotate("text", x = 8, y = 9.7, label = "ROCKIES Trial Design (NCT04027530)",
           size = 4.5, fontface = "bold", family = "sans") +
  annotate("text", x = 8, y = 9.25,
           label = "Randomized, Double-Blind, Placebo-Controlled Crossover",
           size = 3.2, fontface = "italic", family = "sans", color = "gray30") +
  # Screening box
  annotate("rect", xmin = 0.3, xmax = 2.5, ymin = 5.7, ymax = 7.3,
           fill = "#F0F0F0", color = "gray40", linewidth = 0.6) +
  annotate("text", x = 1.4, y = 6.8, label = "Screening &\nRandomization",
           size = 2.8, fontface = "bold", lineheight = 0.9) +
  annotate("text", x = 1.4, y = 6.1, label = "n = 20\nT2D adults",
           size = 2.4, color = "gray30", lineheight = 0.9) +
  # Arrow to R
  annotate("segment", x = 2.5, xend = 3.3, y = 6.5, yend = 6.5,
           arrow = arrow(length = unit(0.15, "cm"), type = "closed"), linewidth = 0.7) +
  # Randomization circle
  annotate("point", x = 3.5, y = 6.5, size = 8, shape = 21,
           fill = "#2C3E50", color = "#2C3E50") +
  annotate("text", x = 3.5, y = 6.5, label = "R",
           size = 3.5, fontface = "bold", color = "white") +
  # --- Sequence 1 (top): Ertu -> Washout -> Placebo ---
  annotate("segment", x = 3.5, xend = 4.2, y = 6.9, yend = 8.0, linewidth = 0.6) +
  annotate("rect", xmin = 4.2, xmax = 7.0, ymin = 7.4, ymax = 8.6,
           fill = "#E74C3C", color = "#C0392B", linewidth = 0.6, alpha = 0.15) +
  annotate("rect", xmin = 4.2, xmax = 7.0, ymin = 7.4, ymax = 8.6,
           fill = NA, color = "#C0392B", linewidth = 0.6) +
  annotate("text", x = 5.6, y = 8.25, label = "Ertugliflozin 15 mg",
           size = 2.8, fontface = "bold", color = "#C0392B") +
  annotate("text", x = 5.6, y = 7.75, label = "4 weeks", size = 2.5, color = "#C0392B") +
  annotate("point", x = 6.7, y = 7.55, size = 2, shape = 17, color = "#C0392B") +
  annotate("segment", x = 7.0, xend = 7.8, y = 8.0, yend = 8.0,
           arrow = arrow(length = unit(0.12, "cm"), type = "closed"),
           linewidth = 0.5, color = "gray40") +
  annotate("rect", xmin = 7.8, xmax = 9.8, ymin = 7.4, ymax = 8.6,
           fill = "#F9F9F9", color = "gray50", linewidth = 0.5, linetype = "dashed") +
  annotate("text", x = 8.8, y = 8.25, label = "Washout",
           size = 2.6, fontface = "italic", color = "gray40") +
  annotate("text", x = 8.8, y = 7.75, label = "6 weeks", size = 2.4, color = "gray40") +
  annotate("segment", x = 9.8, xend = 10.6, y = 8.0, yend = 8.0,
           arrow = arrow(length = unit(0.12, "cm"), type = "closed"),
           linewidth = 0.5, color = "gray40") +
  annotate("rect", xmin = 10.6, xmax = 13.4, ymin = 7.4, ymax = 8.6,
           fill = "#4A90D9", color = "#2E6DA4", linewidth = 0.6, alpha = 0.15) +
  annotate("rect", xmin = 10.6, xmax = 13.4, ymin = 7.4, ymax = 8.6,
           fill = NA, color = "#2E6DA4", linewidth = 0.6) +
  annotate("text", x = 12.0, y = 8.25, label = "Matching Placebo",
           size = 2.8, fontface = "bold", color = "#2E6DA4") +
  annotate("text", x = 12.0, y = 7.75, label = "4 weeks", size = 2.5, color = "#2E6DA4") +
  annotate("point", x = 13.1, y = 7.55, size = 2, shape = 17, color = "#2E6DA4") +
  # --- Sequence 2 (bottom): Placebo -> Washout -> Ertu ---
  annotate("segment", x = 3.5, xend = 4.2, y = 6.1, yend = 5.0, linewidth = 0.6) +
  annotate("rect", xmin = 4.2, xmax = 7.0, ymin = 4.4, ymax = 5.6,
           fill = "#4A90D9", color = "#2E6DA4", linewidth = 0.6, alpha = 0.15) +
  annotate("rect", xmin = 4.2, xmax = 7.0, ymin = 4.4, ymax = 5.6,
           fill = NA, color = "#2E6DA4", linewidth = 0.6) +
  annotate("text", x = 5.6, y = 5.25, label = "Matching Placebo",
           size = 2.8, fontface = "bold", color = "#2E6DA4") +
  annotate("text", x = 5.6, y = 4.75, label = "4 weeks", size = 2.5, color = "#2E6DA4") +
  annotate("point", x = 6.7, y = 4.55, size = 2, shape = 17, color = "#2E6DA4") +
  annotate("segment", x = 7.0, xend = 7.8, y = 5.0, yend = 5.0,
           arrow = arrow(length = unit(0.12, "cm"), type = "closed"),
           linewidth = 0.5, color = "gray40") +
  annotate("rect", xmin = 7.8, xmax = 9.8, ymin = 4.4, ymax = 5.6,
           fill = "#F9F9F9", color = "gray50", linewidth = 0.5, linetype = "dashed") +
  annotate("text", x = 8.8, y = 5.25, label = "Washout",
           size = 2.6, fontface = "italic", color = "gray40") +
  annotate("text", x = 8.8, y = 4.75, label = "6 weeks", size = 2.4, color = "gray40") +
  annotate("segment", x = 9.8, xend = 10.6, y = 5.0, yend = 5.0,
           arrow = arrow(length = unit(0.12, "cm"), type = "closed"),
           linewidth = 0.5, color = "gray40") +
  annotate("rect", xmin = 10.6, xmax = 13.4, ymin = 4.4, ymax = 5.6,
           fill = "#E74C3C", color = "#C0392B", linewidth = 0.6, alpha = 0.15) +
  annotate("rect", xmin = 10.6, xmax = 13.4, ymin = 4.4, ymax = 5.6,
           fill = NA, color = "#C0392B", linewidth = 0.6) +
  annotate("text", x = 12.0, y = 5.25, label = "Ertugliflozin 15 mg",
           size = 2.8, fontface = "bold", color = "#C0392B") +
  annotate("text", x = 12.0, y = 4.75, label = "4 weeks", size = 2.5, color = "#C0392B") +
  annotate("point", x = 13.1, y = 4.55, size = 2, shape = 17, color = "#C0392B") +
  # --- End-of-Period Assessments box ---
  annotate("segment", x = 13.4, xend = 14.0, y = 8.0, yend = 6.9, linewidth = 0.5, color = "gray40") +
  annotate("segment", x = 13.4, xend = 14.0, y = 5.0, yend = 6.1, linewidth = 0.5, color = "gray40") +
  annotate("rect", xmin = 13.8, xmax = 15.8, ymin = 5.0, ymax = 8.0,
           fill = "#F5F5DC", color = "#8B7D3C", linewidth = 0.6) +
  annotate("text", x = 14.8, y = 7.65, label = "End-of-Period",
           size = 2.5, fontface = "bold", color = "#5D4E37") +
  annotate("text", x = 14.8, y = 7.3, label = "Assessments",
           size = 2.5, fontface = "bold", color = "#5D4E37") +
  annotate("text", x = 14.8, y = 6.4,
           label = "\u2022 \u00B9\u00B9C-Acetate PET/CT\n\u2022 OGTT\n\u2022 mGFR (iohexol)\n\u2022 ERPF (PAH)\n\u2022 24-hr urine",
           size = 2.1, lineheight = 1.1, hjust = 0.5, color = "#5D4E37") +
  # --- Labels ---
  annotate("text", x = 3.9, y = 8.5, label = "Sequence 1",
           size = 2.3, fontface = "italic", color = "gray50") +
  annotate("text", x = 3.9, y = 4.5, label = "Sequence 2",
           size = 2.3, fontface = "italic", color = "gray50") +
  annotate("text", x = 5.6, y = 3.7, label = "Period 1",
           size = 2.8, fontface = "bold", color = "gray40") +
  annotate("text", x = 8.8, y = 3.7, label = "Washout",
           size = 2.8, fontface = "bold", color = "gray40") +
  annotate("text", x = 12.0, y = 3.7, label = "Period 2",
           size = 2.8, fontface = "bold", color = "gray40") +
  # Timeline
  annotate("segment", x = 4.2, xend = 13.4, y = 3.3, yend = 3.3, linewidth = 0.6, color = "gray50") +
  annotate("segment", x = 4.2, xend = 4.2, y = 3.15, yend = 3.45, linewidth = 0.4, color = "gray50") +
  annotate("segment", x = 7.0, xend = 7.0, y = 3.15, yend = 3.45, linewidth = 0.4, color = "gray50") +
  annotate("segment", x = 7.8, xend = 7.8, y = 3.15, yend = 3.45, linewidth = 0.4, color = "gray50") +
  annotate("segment", x = 9.8, xend = 9.8, y = 3.15, yend = 3.45, linewidth = 0.4, color = "gray50") +
  annotate("segment", x = 10.6, xend = 10.6, y = 3.15, yend = 3.45, linewidth = 0.4, color = "gray50") +
  annotate("segment", x = 13.4, xend = 13.4, y = 3.15, yend = 3.45, linewidth = 0.4, color = "gray50") +
  annotate("text", x = 4.2, y = 2.9, label = "Wk 0", size = 2.2, color = "gray50") +
  annotate("text", x = 7.0, y = 2.9, label = "Wk 4", size = 2.2, color = "gray50") +
  annotate("text", x = 9.8, y = 2.9, label = "Wk 10", size = 2.2, color = "gray50") +
  annotate("text", x = 13.4, y = 2.9, label = "Wk 14", size = 2.2, color = "gray50") +
  annotate("point", x = 5.0, y = 2.3, size = 2, shape = 17, color = "gray40") +
  annotate("text", x = 5.5, y = 2.3, label = "= End-of-period assessment visit",
           size = 2.2, hjust = 0, color = "gray40") +
  theme_void() +
  theme(plot.margin = margin(5, 5, 5, 5))

# =====================================================================
# Panel B: PET Methodology Diagram (external PNG)
# =====================================================================
library(png)

pet_img <- readPNG("C:/Users/netio/Downloads/Kidney C11 Acetate PET Diagram.png")

fig1b <- ggplot() +
  annotation_raster(pet_img, xmin = 0, xmax = 1, ymin = 0, ymax = 1) +
  xlim(0, 1) + ylim(0, 1) +
  theme_void() +
  theme(plot.margin = margin(5, 5, 5, 5))

# --- Panel C: Cortical k2 vs HOMA-IR (placebo) ---
fig1c <- plot_corr(rockies_plac, "homa_ir", "cortical_k2",
                   "HOMA-IR",
                   expression(bold("Cortical k"[2]*" (min"^{-1}*")")),
                   "C", point_col = cols$placebo)

# --- Panel D: Cortical k2/F vs HOMA-IR (placebo) ---
fig1d <- plot_corr(rockies_plac, "homa_ir", "cortical_k2f",
                   "HOMA-IR",
                   expression(bold("Cortical k"[2]*"/F")),
                   "D", point_col = cols$placebo)

# --- Panel E: Medullary k2 vs Tubular Na reabsorption (placebo) ---
fig1e <- plot_corr(rockies_plac, "tna_sodium", "medullary_k2",
                   expression(bold("Tubular Na"^"+"*" Reabsorption (mmol/min)")),
                   expression(bold("Medullary k"[2]*" (min"^{-1}*")")),
                   "E", point_col = cols$placebo)

# --- Panel F: Placeholder (urinary metabolomics, separate data) ---
fig1f <- ggplot() +
  annotate("text", x = 0.5, y = 0.5,
           label = "Panel F\nUrinary Metabolomics\n(separate data needed)",
           size = 4, color = "gray50") +
  theme_void() + labs(tag = "F")

# --- Panel G: Cortical k2 paired ---
fig1g <- plot_paired(rockies_long, "id", "treatment", "cortical_k2",
                     expression(bold("Cortical k"[2]*" (min"^{-1}*")")), "G")

# --- Panel H: Medullary k2 paired ---
fig1h <- plot_paired(rockies_long, "id", "treatment", "medullary_k2",
                     expression(bold("Medullary k"[2]*" (min"^{-1}*")")), "H")

# --- Panel I: Delta OGIS vs Delta cortical k2 ---
fig1i <- plot_delta_corr(rockies_long, "id", "treatment",
                         "ogis", "cortical_k2",
                         "OGIS (ml/min/m\u00B2)",
                         expression("Cortical k"[2]*" (min"^{-1}*")"), "I")

# --- Panel J: Delta OGIS vs Delta medullary k2 ---
fig1j <- plot_delta_corr(rockies_long, "id", "treatment",
                         "ogis", "medullary_k2",
                         "OGIS (ml/min/m\u00B2)",
                         expression("Medullary k"[2]*" (min"^{-1}*")"), "J")

# --- Assemble Figure 1 ---
# Row 1: A (study design) and B (PET diagram)
# Row 2: C, D, E (placebo correlations)
# Row 3: G, H and placeholder for F
# Row 4: I, J (delta correlations)

row1 <- (fig1a + labs(tag = "A")) | (fig1b + labs(tag = "B"))
row2 <- fig1c | fig1d | fig1e
row3 <- fig1f | fig1g | fig1h
row4 <- fig1i | fig1j | plot_spacer()

fig1 <- row1 / row2 / row3 / row4 +
  plot_layout(heights = c(2.5, 2, 2, 2)) +
  plot_annotation(
    title = "Figure 1. SGLT2 Inhibition Reduces Kidney Oxidative Metabolism in the ROCKIES Trial",
    theme = theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5))
  )

save_fig(fig1, "Figure1_ROCKIES", width = 16, height = 18)

# Save individual panels for Illustrator
save_fig(fig1c, "Figure1C_cortical_k2_vs_HOMAIR", width = 5, height = 4.5)
save_fig(fig1d, "Figure1D_cortical_k2F_vs_HOMAIR", width = 5, height = 4.5)
save_fig(fig1e, "Figure1E_medullary_k2_vs_TNa", width = 5, height = 4.5)
save_fig(fig1g, "Figure1G_cortical_k2_paired", width = 4.5, height = 5)
save_fig(fig1h, "Figure1H_medullary_k2_paired", width = 4.5, height = 5)
save_fig(fig1i, "Figure1I_delta_OGIS_vs_delta_cortical_k2", width = 5, height = 4.5)
save_fig(fig1j, "Figure1J_delta_OGIS_vs_delta_medullary_k2", width = 5, height = 4.5)

cat("Figure 1 saved!\n")

########################################################################
# =====================================================================
#  FIGURE 2: T2D vs HEALTHY CONTROLS — k2, F, k2/F
# =====================================================================
########################################################################

dat_fig2 <- dat_with_pet %>%
  filter(!is.na(avg_c_k2) & is.finite(avg_c_k2)) %>%
  filter(group %in% c('Lean Control', 'Type 2 Diabetes')) %>%
  mutate(Cohort = case_when(
    group == "Type 2 Diabetes" ~ "T2D\n(RENAL-HEIRitage)",
    group == "Lean Control"    ~ "Healthy Control\n(CROCODILE)"
  ))

cat("Figure 2 data:", nrow(dat_fig2), "rows\n")

fig2_colors <- c("Healthy Control\n(CROCODILE)" = cols$control,
                 "T2D\n(RENAL-HEIRitage)" = cols$t2d)

fig2b <- plot_two_group(dat_fig2, "Cohort", "avg_c_k2", fig2_colors,
                        ylab = expression(bold("Cortical k"[2]*" (min"^{-1}*")")), tag = "B")
fig2c <- plot_two_group(dat_fig2, "Cohort", "avg_c_f", fig2_colors,
                        ylab = expression(bold("Cortical F (ml/g/min)")), tag = "C")
fig2d <- plot_two_group(dat_fig2, "Cohort", "avg_c_k2_f", fig2_colors,
                        ylab = expression(bold("Cortical k"[2]*"/F")), tag = "D")

fig2 <- (fig2b | fig2c | fig2d) +
  plot_annotation(
    title = "Figure 2. Elevated Kidney Oxidative Metabolism in Type 2 Diabetes",
    theme = theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5))
  )

save_fig(fig2, "Figure2_T2D_vs_Controls", width = 12, height = 5)
cat("Figure 2 saved!\n")

########################################################################
# Figure 3 update: use log10(UACR) on x-axis
# 
# Replace the fig3b/fig3c/fig3d lines in your main script with these.
# The Spearman correlation is computed on the ORIGINAL acr_u values
# (rank-based, so log transform doesn't change rho/p),
# but the x-axis is displayed as log10(UACR).
########################################################################

# --- Add log10 UACR to the data ---
dat_fig3 <- dat_fig3 %>%
  mutate(log10_uacr = log10(acr_u))

# --- Updated plot_corr that uses log10 x but computes rho on raw x ---
plot_corr_log <- function(data, xvar_raw, xvar_log, yvar, xlab, ylab, tag,
                          point_col = cols$corr_point, method = "spearman") {
  df <- data.frame(x_raw = data[[xvar_raw]], 
                   x_log = data[[xvar_log]], 
                   y = data[[yvar]])
  df <- df[!is.na(df$x_raw) & !is.na(df$x_log) & !is.na(df$y) & 
             is.finite(df$x_log), ]
  
  # Spearman on original values (rank-based, same result either way)
  ct <- cor.test(df$x_raw, df$y, method = method)
  lab <- paste0(ifelse(method == "spearman", "rho", "r"), "=",
                formatC(ct$estimate, format = "f", digits = 2),
                ", p=", ifelse(ct$p.value < 0.001, "<0.001",
                               formatC(ct$p.value, format = "f", digits = 3)))
  
  ggplot(df, aes(x_log, y)) +
    geom_point(size = 2.5, alpha = 0.7, color = point_col) +
    geom_smooth(method = "lm", se = TRUE, color = "black",
                linewidth = 0.8, fill = "gray85", formula = y ~ x) +
    annotate("text", x = min(df$x_log) + diff(range(df$x_log)) * 0.02,
             y = max(df$y) - diff(range(df$y)) * 0.02,
             label = lab, hjust = 0, size = 3, fontface = "italic") +
    labs(x = xlab, y = ylab, tag = tag) +
    theme_rockies
}

# --- Panels B-D with log10(UACR) ---
fig3b <- plot_corr_log(dat_fig3, "acr_u", "log10_uacr", "avg_c_k2",
                       expression(bold("log"[10]*"(UACR) (mg/g)")),
                       expression(bold("Cortical k"[2]*" (min"^{-1}*")")),
                       "B", point_col = cols$corr_point)

fig3c <- plot_corr_log(dat_fig3, "acr_u", "log10_uacr", "avg_c_f",
                       expression(bold("log"[10]*"(UACR) (mg/g)")),
                       expression(bold("Cortical F (ml/g/min)")),
                       "C", point_col = cols$corr_point)

fig3d <- plot_corr_log(dat_fig3, "acr_u", "log10_uacr", "avg_c_k2_f",
                       expression(bold("log"[10]*"(UACR) (mg/g)")),
                       expression(bold("Cortical k"[2]*"/F")),
                       "D", point_col = cols$corr_point)

# --- Also update Panel A heatmap to note log10 on x-axis label ---
fig3a <- ggplot(cor_results, aes(x = "UACR", y = var_label, fill = rho)) +
  geom_tile(color = "white", linewidth = 1) +
  geom_text(aes(label = paste0(label, sig)), size = 4, fontface = "bold") +
  scale_fill_gradient2(low = cols$control, mid = "white", high = cols$t2d,
                       midpoint = 0, limits = c(-0.6, 0.6),
                       name = "Spearman\nrho") +
  labs(x = NULL, y = NULL, tag = "A") +
  theme_rockies +
  theme(legend.position = "right",
        axis.text.x = element_text(face = "bold"))

# --- Reassemble Figure 3 ---
fig3 <- (fig3a | fig3b) / (fig3c | fig3d) +
  plot_annotation(
    title = "Figure 3. Kidney Oxidative Metabolism Correlates with Albuminuria",
    subtitle = "Spearman correlations in 40 participants (CROCODILE + RENAL-HEIRitage)",
    theme = theme(
      plot.title    = element_text(size = 14, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 11, hjust = 0.5, color = "gray30")
    )
  )

# Save standalone
save_fig(fig3, "Figure3_UACR_Correlation", width = 12, height = 8)
cat("Figure 3 updated with log10(UACR)\n")

########################################################################
# =====================================================================
#  FIGURE 4: GBM THICKENING — k2, F, k2/F
# =====================================================================
########################################################################

dat_fig4_clin <- dat_clinical %>%
  as.data.frame() %>%
  filter(`GBM thickness` != '' & !is.na(`GBM thickness`)) %>%
  mutate(GBM_Status = ifelse(`GBM thickness` == "yes",
                             "GBM\nThickening", "No GBM\nThickening")) %>%
  dplyr::select(record_id, GBM_Status)

# Drop any existing PET columns from clinical before merging
dat_fig4_clin <- dat_fig4_clin[, !(names(dat_fig4_clin) %in% pet_col_names), drop = FALSE]
dat_fig4 <- merge(dat_fig4_clin, dat_pet_slim, by = "record_id")

cat("Figure 4 data:", nrow(dat_fig4), "rows\n")
cat("  Columns:", paste(names(dat_fig4), collapse = ", "), "\n")
cat("  GBM groups:\n"); print(table(dat_fig4$GBM_Status))

fig4_colors <- c("No GBM\nThickening" = cols$gbm_no, "GBM\nThickening" = cols$gbm_yes)

fig4a <- plot_two_group(dat_fig4, "GBM_Status", "avg_c_k2", fig4_colors,
                        ylab = expression(bold("Cortical k"[2]*" (min"^{-1}*")")), tag = "A")
fig4b <- plot_two_group(dat_fig4, "GBM_Status", "avg_c_f", fig4_colors,
                        ylab = expression(bold("Cortical F (ml/g/min)")), tag = "B")
fig4c <- plot_two_group(dat_fig4, "GBM_Status", "avg_c_k2_f", fig4_colors,
                        ylab = expression(bold("Cortical k"[2]*"/F")), tag = "C")

fig4 <- (fig4a | fig4b | fig4c) +
  plot_annotation(
    title = "Figure 4. Kidney Oxidative Metabolism and GBM Thickening",
    theme = theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5))
  )

save_fig(fig4, "Figure4_GBM", width = 12, height = 5)
cat("Figure 4 saved!\n")

########################################################################
# =====================================================================
#  FIGURE 5: ARTERIOSCLEROSIS — k2, F, k2/F
# =====================================================================
########################################################################

dat_fig5_clin <- dat_clinical %>%
  as.data.frame() %>%
  filter(arteriosclerosis != '' & !is.na(arteriosclerosis)) %>%
  mutate(Arterio_Status = ifelse(arteriosclerosis == "yes",
                                 "Arteriosclerosis", "No\nArteriosclerosis")) %>%
  dplyr::select(record_id, Arterio_Status)

dat_fig5_clin <- dat_fig5_clin[, !(names(dat_fig5_clin) %in% pet_col_names), drop = FALSE]
dat_fig5 <- merge(dat_fig5_clin, dat_pet_slim, by = "record_id")

cat("Figure 5 data:", nrow(dat_fig5), "rows\n")
cat("  Columns:", paste(names(dat_fig5), collapse = ", "), "\n")
cat("  Arterio groups:\n"); print(table(dat_fig5$Arterio_Status))

fig5_colors <- c("No\nArteriosclerosis" = cols$arterio_no,
                 "Arteriosclerosis" = cols$arterio_yes)

fig5a <- plot_two_group(dat_fig5, "Arterio_Status", "avg_c_k2", fig5_colors,
                        ylab = expression(bold("Cortical k"[2]*" (min"^{-1}*")")), tag = "A")
fig5b <- plot_two_group(dat_fig5, "Arterio_Status", "avg_c_f", fig5_colors,
                        ylab = expression(bold("Cortical F (ml/g/min)")), tag = "B")
fig5c <- plot_two_group(dat_fig5, "Arterio_Status", "avg_c_k2_f", fig5_colors,
                        ylab = expression(bold("Cortical k"[2]*"/F")), tag = "C")

fig5 <- (fig5a | fig5b | fig5c) +
  plot_annotation(
    title = "Figure 5. Kidney Oxidative Metabolism and Arteriosclerosis",
    theme = theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5))
  )

save_fig(fig5, "Figure5_Arteriosclerosis", width = 12, height = 5)
cat("Figure 5 saved!\n")

########################################################################
# =====================================================================
#  FIGURES 6-7: scRNAseq — TCA Cycle & OxPhos Gene Expression
#  These are bar plots (log2 FC, colored by direction, * for sig)
#  Generated by existing scripts — load PNGs for combined figure
# =====================================================================
########################################################################

# If you already have Figures 6 and 7 saved from your existing scripts,
# load them here for reference. Update paths as needed:

# fig6_path <- paste0(base_path, "Figure6_TCA_barplots.png")
# fig7_path <- paste0(base_path, "Figure7_OxPhos_barplots.png")

# If you want to regenerate them here in the same style, use this template:
# It expects a data frame with columns: gene, log2FC, padj

plot_gene_bars <- function(de_results, gene_list, title, tag, padj_thresh = 0.05) {
  df <- de_results %>%
    filter(gene %in% gene_list) %>%
    mutate(
      Regulation = ifelse(log2FC >= 0, "Up-regulated", "Down-regulated"),
      Regulation = factor(Regulation, levels = c("Down-regulated", "Up-regulated")),
      sig_label  = ifelse(padj < padj_thresh, "*", ""),
      gene       = factor(gene, levels = gene[order(match(gene, gene_list))])
    )
  
  ggplot(df, aes(x = gene, y = log2FC, fill = Regulation)) +
    geom_col(width = 0.7) +
    geom_text(aes(label = sig_label),
              vjust = ifelse(df$log2FC >= 0, -0.3, 1.3),
              size = 5, fontface = "bold") +
    geom_hline(yintercept = 0, linewidth = 0.4) +
    scale_fill_manual(values = c("Down-regulated" = "#E74C3C",
                                 "Up-regulated" = "#4A90D9")) +
    labs(x = NULL, y = "Log2 Fold Change", title = title, tag = tag) +
    theme_rockies +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.text = element_text(size = 8)
    )
}

tca_genes <- c("ACO1","ACO2","CS","DHTKD1","DLAT","DLD","DLST","FH",
               "IDH1","IDH2","IDH3A","IDH3B","IDH3G","MDH1","MDH2",
               "MPC1","OGDH","OGDHL","PDHA1","PDHB","SDHA","SDHB",
               "SDHC","SDHD","SUCLA2","SUCLG1","SUCLG2")

oxphos_genes <- c("NDUFA1","NDUFA2","NDUFA3","NDUFA4","NDUFA5","NDUFA6",
                  "NDUFA7","NDUFA8","NDUFA9","NDUFA10","NDUFA11","NDUFA12","NDUFA13",
                  "NDUFAB1","NDUFB1","NDUFB2","NDUFB3","NDUFB4","NDUFB5","NDUFB6",
                  "NDUFB7","NDUFB8","NDUFB9","NDUFB10","NDUFB11","NDUFC1","NDUFC2",
                  "NDUFS1","NDUFS2","NDUFS3","NDUFS4","NDUFS5","NDUFS6","NDUFS7",
                  "NDUFS8","NDUFV1","NDUFV2","NDUFV3","SDHA","SDHB","SDHC","SDHD",
                  "UQCRB","UQCRC1","UQCRC2","UQCRFS1","UQCRH","UQCRQ","CYC1",
                  "UQCR10","UQCR11","COX4I1","COX4I2","COX5A","COX5B","COX6A1",
                  "COX6A2","COX6B1","COX6C","COX7A1","COX7A2","COX7B","COX7C","COX8A",
                  "MT-CO1","MT-CO2","MT-CO3","ATP5F1A","ATP5F1B","ATP5F1C","ATP5F1D",
                  "ATP5F1E","ATP5MC1","ATP5MC2","ATP5MC3","ATP5MF","ATP5MG","ATP5PB",
                  "ATP5PD","ATP5PF","ATP5PO")

# Uncomment when NEBULA results are available:
# nebula_pt_all  <- read.csv("path/to/nebula_total_PT.csv")
# nebula_pt_s1s2 <- read.csv("path/to/nebula_PT_S1S2.csv")
# nebula_pt_s3   <- read.csv("path/to/nebula_PT_S3.csv")
# nebula_apt     <- read.csv("path/to/nebula_aPT.csv")
#
# # Figure 6: TCA Cycle
# fig6a <- plot_gene_bars(nebula_pt_all,  tca_genes, "PT",      "A")
# fig6b <- plot_gene_bars(nebula_pt_s1s2, tca_genes, "PT-S1/S2","B")
# fig6c <- plot_gene_bars(nebula_pt_s3,   tca_genes, "PT-S3",   "C")
# fig6d <- plot_gene_bars(nebula_apt,     tca_genes, "aPT",     "D")
#
# fig6 <- (fig6a) / (fig6b | fig6c | fig6d) +
#   plot_layout(heights = c(1.2, 1)) +
#   plot_annotation(
#     title = "Figure 6. TCA Cycle Gene Expression in Proximal Tubule Cells",
#     subtitle = "TCA cycle gene expression in proximal tubule cells",
#     theme = theme(plot.title = element_text(size = 13, face = "bold"),
#                   plot.subtitle = element_text(size = 10)))
# save_fig(fig6, "Figure6_TCA_barplots", width = 14, height = 10)
#
# # Figure 7: OxPhos
# fig7a <- plot_gene_bars(nebula_pt_all,  oxphos_genes, "PT",      "A")
# fig7b <- plot_gene_bars(nebula_pt_s1s2, oxphos_genes, "PT-S1/S2","B")
# fig7c <- plot_gene_bars(nebula_pt_s3,   oxphos_genes, "PT-S3",   "C")
# fig7d <- plot_gene_bars(nebula_apt,     oxphos_genes, "aPT",     "D")
#
# fig7 <- (fig7a) / (fig7b | fig7c | fig7d) +
#   plot_layout(heights = c(1.2, 1)) +
#   plot_annotation(
#     title = "Figure 7. OxPhos Gene Expression in Proximal Tubule Cells",
#     subtitle = "Oxidative phosphorylation gene expression in proximal tubule cells",
#     theme = theme(plot.title = element_text(size = 13, face = "bold"),
#                   plot.subtitle = element_text(size = 10)))
# save_fig(fig7, "Figure7_OxPhos_barplots", width = 16, height = 10)

########################################################################
# STATUS
########################################################################

cat("\n================================================================\n")
cat("  ROCKIES PUBLICATION FIGURES - STATUS\n")
cat("================================================================\n")
cat("  FIGURE 1 (ROCKIES):            --> GENERATED\n")
cat("    Panels C: Cortical k2 vs HOMA-IR\n")
cat("    Panels D: Cortical k2/F vs HOMA-IR\n")
cat("    Panels E: Medullary k2 vs TNa reabsorption\n")
cat("    Panels G: Cortical k2 paired\n")
cat("    Panels H: Medullary k2 paired\n")
cat("    Panels I: Delta OGIS vs Delta cortical k2\n")
cat("    Panels J: Delta OGIS vs Delta medullary k2\n")
cat("    Panel A,B: Study design       --> Illustrator\n")
cat("    Panel F: Metabolomics         --> Separate data needed\n")
cat("  FIGURE 2 (T2D vs Controls):    --> GENERATED\n")
cat("  FIGURE 3 (UACR Correlations):  --> GENERATED\n")
cat("  FIGURE 4 (GBM Thickening):     --> GENERATED\n")
cat("  FIGURE 5 (Arteriosclerosis):   --> GENERATED\n")
cat("  FIGURE 6 (TCA Barplots):       --> Use existing script / plot_gene_bars()\n")
cat("  FIGURE 7 (OxPhos Barplots):     --> Use existing script / plot_gene_bars()\n")
cat("================================================================\n")






###combine the files 


########################################################################
# ROCKIES Publication - Combined PDF (All Figures, One Page Each)
# Run AFTER the main figure generation script
########################################################################

library(tidyverse)
library(patchwork)
library(png)

base_path <- 'C:/Users/netio/Documents/UofW/Rockies/publication_figures/'

# Page dimensions (letter-ish, landscape for wider figures)
pg_w <- 16
pg_h <- 18

########################################################################
# Open multi-page PDF
########################################################################

tryCatch({
  cairo_pdf(paste0(base_path, "ROCKIES_All_Figures_Combined.pdf"),
            width = pg_w, height = pg_h, onefile = TRUE)
  pdf_device <- "cairo"
}, error = function(e) {
  pdf(paste0(base_path, "ROCKIES_All_Figures_Combined.pdf"),
      width = pg_w, height = pg_h, onefile = TRUE)
  pdf_device <- "standard"
})

########################################################################
# PAGE 1: Figure 1 — ROCKIES Trial
########################################################################

# Row 1: A (study design) + B (PET diagram)
row1 <- (fig1a + labs(tag = "A")) | (fig1b + labs(tag = "B"))

# Row 2: Placebo correlations
row2 <- fig1c | fig1d | fig1e

# Row 3: Paired plots (+ placeholder for F)
row3 <- fig1f | fig1g | fig1h

# Row 4: Delta correlations
row4 <- fig1i | fig1j | plot_spacer()

fig1_full <- row1 / row2 / row3 / row4 +
  plot_layout(heights = c(2.5, 2, 2, 2)) +
  plot_annotation(
    title = "Figure 1. SGLT2 Inhibition Reduces Kidney Oxidative Metabolism in the ROCKIES Trial",
    theme = theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5))
  )

print(fig1_full)

########################################################################
# PAGE 2: Figure 2 — T2D vs Healthy Controls
########################################################################

fig2_full <- (fig2b | fig2c | fig2d) +
  plot_annotation(
    title = "Figure 2. Elevated Kidney Oxidative Metabolism in Type 2 Diabetes",
    subtitle = "18 T2D (RENAL-HEIRitage) vs 11 Healthy Controls (CROCODILE)",
    theme = theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 11, hjust = 0.5, color = "gray30")
    )
  )

print(fig2_full)

########################################################################
# PAGE 3: Figure 3 — UACR Correlations
########################################################################

fig3_full <- (fig3a | fig3b) / (fig3c | fig3d) +
  plot_annotation(
    title = "Figure 3. Kidney Oxidative Metabolism Correlates with Albuminuria",
    subtitle = "Spearman correlations in 40 participants (CROCODILE + RENAL-HEIRitage)",
    theme = theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 11, hjust = 0.5, color = "gray30")
    )
  )

print(fig3_full)

########################################################################
# PAGE 4: Figure 4 — GBM Thickening
########################################################################

fig4_full <- (fig4a | fig4b | fig4c) +
  plot_annotation(
    title = "Figure 4. Kidney Oxidative Metabolism and Glomerular Basement Membrane Thickening",
    subtitle = paste0("GBM Thickening (n=", sum(dat_fig4$GBM_Status == "GBM\nThickening"),
                      ") vs No GBM Thickening (n=", sum(dat_fig4$GBM_Status == "No GBM\nThickening"), ")"),
    theme = theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 11, hjust = 0.5, color = "gray30")
    )
  )

print(fig4_full)

########################################################################
# PAGE 5: Figure 5 — Arteriosclerosis
########################################################################

fig5_full <- (fig5a | fig5b | fig5c) +
  plot_annotation(
    title = "Figure 5. Kidney Oxidative Metabolism and Arteriosclerosis",
    subtitle = paste0("Arteriosclerosis (n=", sum(dat_fig5$Arterio_Status == "Arteriosclerosis"),
                      ") vs No Arteriosclerosis (n=", sum(dat_fig5$Arterio_Status == "No\nArteriosclerosis"), ")"),
    theme = theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 11, hjust = 0.5, color = "gray30")
    )
  )

print(fig5_full)

########################################################################
# PAGES 6-7: Figures 6 & 7 — scRNAseq (if available)
# If these exist as saved PNGs from your existing script, load them.
# Otherwise these pages are skipped.
########################################################################

# Try loading existing Figure 6 PNG
fig6_png_path <- paste0(base_path, "Figure6_TCA_Cycle.png")
if (file.exists(fig6_png_path)) {
  fig6_img <- readPNG(fig6_png_path)
  fig6_page <- ggplot() +
    annotation_raster(fig6_img, xmin = 0, xmax = 1, ymin = 0, ymax = 1) +
    xlim(0, 1) + ylim(0, 1) +
    theme_void() +
    ggtitle("Figure 6. TCA Cycle Gene Expression in Proximal Tubule Cells") +
    theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5))
  print(fig6_page)
  cat("Figure 6 page added from PNG\n")
} else {
  # If fig6a-d objects exist (generated by plot_gene_bars), use them
  if (exists("fig6a") & exists("fig6b") & exists("fig6c") & exists("fig6d")) {
    fig6_page <- (fig6a) / (fig6b | fig6c | fig6d) +
      plot_layout(heights = c(1.2, 1)) +
      plot_annotation(
        title = "Figure 6. TCA Cycle Gene Expression in Proximal Tubule Cells",
        theme = theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5))
      )
    print(fig6_page)
    cat("Figure 6 page added from plot objects\n")
  } else {
    cat("Figure 6 skipped — no PNG or plot objects found\n")
  }
}

# Try loading existing Figure 7 PNG
fig7_png_path <- paste0(base_path, "Figure7_OxPhos.png")
if (file.exists(fig7_png_path)) {
  fig7_img <- readPNG(fig7_png_path)
  fig7_page <- ggplot() +
    annotation_raster(fig7_img, xmin = 0, xmax = 1, ymin = 0, ymax = 1) +
    xlim(0, 1) + ylim(0, 1) +
    theme_void() +
    ggtitle("Figure 7. OxPhos Gene Expression in Proximal Tubule Cells") +
    theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5))
  print(fig7_page)
  cat("Figure 7 page added from PNG\n")
} else {
  if (exists("fig7a") & exists("fig7b") & exists("fig7c") & exists("fig7d")) {
    fig7_page <- (fig7a) / (fig7b | fig7c | fig7d) +
      plot_layout(heights = c(1.2, 1)) +
      plot_annotation(
        title = "Figure 7. OxPhos Gene Expression in Proximal Tubule Cells",
        theme = theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5))
      )
    print(fig7_page)
    cat("Figure 7 page added from plot objects\n")
  } else {
    cat("Figure 7 skipped — no PNG or plot objects found\n")
  }
}

########################################################################
# Close PDF
########################################################################

dev.off()

cat("\n========================================================\n")
cat("  Combined PDF saved:\n")
cat("  ", paste0(base_path, "ROCKIES_All_Figures_Combined.pdf"), "\n")
cat("========================================================\n")












###Figures 6 and 7 updated
########################################################################
# ROCKIES Publication - Combined PDF (All Figures, One Page Each)
# Regenerates Figures 6 & 7 from CSV data
# Run AFTER the main figure generation script (Figures 1-5)
########################################################################

library(tidyverse)
library(patchwork)
library(cowplot)
library(data.table)
library(png)
library(grid)
library(gridExtra)

base_path <- 'C:/Users/netio/Documents/UofW/Rockies/publication_figures/'
scrna_path <- 'C:/Users/netio/Documents/UofW/Rockies/Hailey_Dotplots/'

pg_w <- 16
pg_h <- 18

########################################################################
# FUNCTION: Build a bar plot for one cell type / pathway
########################################################################

plot_gene_bars <- function(data_file, celltype, panel_tag) {
  
  dt <- fread(data_file)
  
  if ("logFC_groupType_2_Diabetes" %in% names(dt))
    setnames(dt, "logFC_groupType_2_Diabetes", "logFC")
  if ("p_groupType_2_Diabetes" %in% names(dt))
    setnames(dt, "p_groupType_2_Diabetes", "pvalue")
  
  dt <- dt %>%
    mutate(
      direction = factor(
        ifelse(logFC < 0, "Down-regulated", "Up-regulated"),
        levels = c("Up-regulated", "Down-regulated")
      ),
      Significance = ifelse(pvalue < 0.05, "*", "")
    )
  
  max_abs <- max(abs(dt$logFC), na.rm = TRUE)
  dt <- dt %>%
    mutate(ast_y     = ifelse(logFC >= 0, logFC + max_abs * 0.08,
                              logFC - max_abs * 0.08),
           ast_vjust = ifelse(logFC >= 0, 0, 1))
  
  p <- ggplot(dt, aes(x = gene, y = logFC, fill = direction)) +
    geom_col(color = "black", linewidth = 0.3, width = 0.7) +
    geom_text(aes(y = ast_y, label = Significance, vjust = ast_vjust),
              size = 5, show.legend = FALSE) +
    geom_hline(yintercept = 0, linewidth = 0.4) +
    scale_fill_manual(
      name   = "Regulation",
      values = c("Up-regulated" = "#4575b4", "Down-regulated" = "#d73027"),
      drop   = FALSE
    ) +
    scale_y_continuous(expand = expansion(mult = c(0.12, 0.15))) +
    labs(x = NULL, y = "Log2 Fold Change", subtitle = celltype, tag = panel_tag) +
    theme_classic(base_size = 11) +
    theme(
      text              = element_text(color = "black"),
      axis.text.x       = element_text(angle = 90, hjust = 1, vjust = 0.5,
                                       size = 9, color = "black"),
      axis.text.y       = element_text(size = 10, color = "black"),
      axis.title        = element_text(size = 11, face = "bold"),
      axis.title.y      = element_text(margin = margin(r = 8)),
      axis.line         = element_line(linewidth = 0.5),
      plot.subtitle     = element_text(size = 12, face = "bold", hjust = 0.5,
                                       margin = margin(b = 8)),
      plot.tag          = element_text(size = 14, face = "bold"),
      plot.tag.position = c(0.01, 0.98),
      plot.margin       = margin(10, 10, 15, 10),
      legend.position   = "none",
      panel.background  = element_rect(fill = "white", color = NA),
      plot.background   = element_rect(fill = "white", color = NA)
    )
  
  return(p)
}

########################################################################
# HELPER: find the right CSV for a given pathway + cell type
########################################################################

find_csv <- function(pathway, celltype) {
  celltype_safe <- str_replace_all(str_replace_all(celltype, "/", "_"), "-", "_")
  all_csv <- list.files(scrna_path, pattern = "\\.csv$", full.names = TRUE)
  hits <- all_csv[str_detect(all_csv, fixed(pathway)) &
                    str_detect(all_csv, fixed(celltype_safe))]
  if (length(hits) == 0) {
    hits <- all_csv[str_detect(tolower(all_csv), tolower(pathway)) &
                      str_detect(all_csv, fixed(celltype_safe))]
  }
  if (length(hits) == 0) {
    warning("No CSV found for ", pathway, " / ", celltype)
    return(NULL)
  }
  return(hits[1])
}

########################################################################
# Standalone legend as a ggplot (not extracted grob)
# This is placed as a full patchwork panel so boxes render properly
########################################################################

legend_panel <- ggplot(data.frame(
  x = c("Up-regulated", "Down-regulated"),
  y = c(1, 1),
  fill = factor(c("Up-regulated", "Down-regulated"),
                levels = c("Up-regulated", "Down-regulated"))
), aes(x, y, fill = fill)) +
  geom_col(width = 0.5, color = "black", linewidth = 0.5) +
  scale_fill_manual(
    name   = "Regulation",
    values = c("Up-regulated" = "#4575b4", "Down-regulated" = "#d73027")
  ) +
  theme_void() +
  theme(
    legend.position   = "bottom",
    legend.direction  = "horizontal",
    legend.text       = element_text(size = 12, color = "black", margin = margin(l = 4, r = 12)),
    legend.title      = element_text(size = 12, face = "bold", color = "black"),
    legend.key.size   = unit(0.7, "cm"),
    legend.key        = element_rect(colour = "black", linewidth = 0.5),
    legend.margin     = margin(t = 5, b = 5),
    plot.background   = element_rect(fill = "white", color = NA)
  ) +
  guides(fill = guide_legend(
    override.aes = list(colour = "black", linewidth = 0.5)
  ))

# Extract just the legend grob
legend_grob <- get_legend(legend_panel)

# Wrap it as a patchwork-compatible element
legend_wrap <- wrap_elements(full = legend_grob)

########################################################################
# HELPER: assemble 4-panel figure with shared bottom legend
########################################################################

build_figure <- function(plots, title_text, save_name) {
  
  combined <- (plots[["PT"]]) /
    (plots[["PT-S1/S2"]] | plots[["PT-S3"]] | plots[["aPT"]]) /
    legend_wrap +
    plot_layout(heights = c(5, 4, 0.5)) +
    plot_annotation(
      title = title_text,
      theme = theme(
        plot.title = element_text(size = 14, face = "bold",
                                  hjust = 0.5, color = "black")
      )
    )
  
  ggsave(paste0(base_path, save_name, ".pdf"), plot = combined,
         width = 14, height = 11, units = "in", device = cairo_pdf, dpi = 300)
  ggsave(paste0(base_path, save_name, ".png"), plot = combined,
         width = 14, height = 11, units = "in", dpi = 300)
  cat(save_name, "saved (PDF + PNG)\n")
  
  return(combined)
}

########################################################################
# BUILD FIGURE 6: TCA Cycle
########################################################################

celltypes <- c("PT", "PT-S1/S2", "PT-S3", "aPT")
tags       <- c("A", "B", "C", "D")

tca_plots <- list()
for (i in seq_along(celltypes)) {
  csv <- find_csv("TCA", celltypes[i])
  if (!is.null(csv)) {
    tca_plots[[celltypes[i]]] <- plot_gene_bars(csv, celltypes[i], tags[i])
  }
}

fig6_final <- build_figure(
  tca_plots,
  title_text = "Figure 6. TCA Cycle Gene Expression in Proximal Tubule Cells",
  save_name  = "Figure6_TCA_Cycle"
)

########################################################################
# BUILD FIGURE 7: OxPhos
########################################################################

oxphos_plots <- list()
for (i in seq_along(celltypes)) {
  csv <- find_csv("OX_PHOS", celltypes[i])
  if (!is.null(csv)) {
    oxphos_plots[[celltypes[i]]] <- plot_gene_bars(csv, celltypes[i], tags[i])
  }
}

fig7_final <- build_figure(
  oxphos_plots,
  title_text = "Figure 7. OxPhos Gene Expression in Proximal Tubule Cells",
  save_name  = "Figure7_OxPhos"
)

########################################################################
# OPEN COMBINED MULTI-PAGE PDF
########################################################################

out_file <- paste0(base_path, "ROCKIES_All_Figures_Combined.pdf")

tryCatch({
  cairo_pdf(out_file, width = pg_w, height = pg_h, onefile = TRUE)
  cat("Using cairo_pdf device\n")
}, error = function(e) {
  pdf(out_file, width = pg_w, height = pg_h, onefile = TRUE)
  cat("Using standard pdf device\n")
})

########################################################################
# PAGE 1: Figure 1
########################################################################
row1 <- (fig1a + labs(tag = "A")) | (fig1b + labs(tag = "B"))
row2 <- fig1c | fig1d | fig1e
row3 <- fig1f | fig1g | fig1h
row4 <- fig1i | fig1j | plot_spacer()

fig1_full <- row1 / row2 / row3 / row4 +
  plot_layout(heights = c(4, 2, 2, 2)) +
  plot_annotation(
    title = "Figure 1. SGLT2 Inhibition Reduces Kidney Oxidative Metabolism in the ROCKIES Trial",
    theme = theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5))
  )
print(fig1_full)

########################################################################
# PAGE 2: Figure 2
########################################################################
fig2_full <- (fig2b | fig2c | fig2d) +
  plot_annotation(
    title = "Figure 2. Elevated Kidney Oxidative Metabolism in Type 2 Diabetes",
    subtitle = "18 T2D (RENAL-HEIRitage) vs 11 Healthy Controls (CROCODILE)",
    theme = theme(
      plot.title    = element_text(size = 14, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 11, hjust = 0.5, color = "gray30")
    )
  )
print(fig2_full)

########################################################################
# PAGE 3: Figure 3
########################################################################
fig3_full <- (fig3a | fig3b) / (fig3c | fig3d) +
  plot_annotation(
    title = "Figure 3. Kidney Oxidative Metabolism Correlates with Albuminuria",
    subtitle = "Spearman correlations in 40 participants (CROCODILE + RENAL-HEIRitage)",
    theme = theme(
      plot.title    = element_text(size = 14, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 11, hjust = 0.5, color = "gray30")
    )
  )
print(fig3_full)

########################################################################
# PAGE 4: Figure 4
########################################################################
fig4_full <- (fig4a | fig4b | fig4c) +
  plot_annotation(
    title = "Figure 4. Kidney Oxidative Metabolism and Glomerular Basement Membrane Thickening",
    subtitle = paste0(
      "GBM Thickening (n=", sum(dat_fig4$GBM_Status == "GBM\nThickening"),
      ") vs No GBM Thickening (n=", sum(dat_fig4$GBM_Status == "No GBM\nThickening"), ")"
    ),
    theme = theme(
      plot.title    = element_text(size = 14, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 11, hjust = 0.5, color = "gray30")
    )
  )
print(fig4_full)

########################################################################
# PAGE 5: Figure 5
########################################################################
fig5_full <- (fig5a | fig5b | fig5c) +
  plot_annotation(
    title = "Figure 5. Kidney Oxidative Metabolism and Arteriosclerosis",
    subtitle = paste0(
      "Arteriosclerosis (n=", sum(dat_fig5$Arterio_Status == "Arteriosclerosis"),
      ") vs No Arteriosclerosis (n=", sum(dat_fig5$Arterio_Status == "No\nArteriosclerosis"), ")"
    ),
    theme = theme(
      plot.title    = element_text(size = 14, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 11, hjust = 0.5, color = "gray30")
    )
  )
print(fig5_full)



########################################################################
# PAGE 6: Figure 6 — TCA Cycle
########################################################################
print(fig6_final)

########################################################################
# PAGE 7: Figure 7 — OxPhos
########################################################################
print(fig7_final)

########################################################################
# Close PDF
########################################################################
dev.off()

cat("\n========================================================\n")
cat("  Combined PDF saved:\n")
cat("  ", out_file, "\n")
cat("========================================================\n")
