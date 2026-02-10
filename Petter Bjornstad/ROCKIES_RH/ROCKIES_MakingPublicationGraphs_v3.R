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
fig1 <- (fig1c | fig1d | fig1e) /
  (fig1f | fig1g | fig1h) /
  (fig1i | fig1j | plot_spacer()) +
  plot_annotation(
    title = "Figure 1. SGLT2 Inhibition Reduces Kidney Oxidative Metabolism in the ROCKIES Trial",
    theme = theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5))
  )

save_fig(fig1, "Figure1_ROCKIES", width = 14, height = 13)

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
# =====================================================================
#  FIGURE 3: UACR CORRELATIONS — k2, F, k2/F
# =====================================================================
########################################################################

dat_fig3 <- dat_with_pet %>%
  filter(!is.na(avg_c_k2) & is.finite(avg_c_k2)) %>%
  filter(group %in% c('Lean Control', 'Obese Control', 'Type 2 Diabetes')) %>%
  filter(record_id != 'CRC-55')

cat("Figure 3 data:", nrow(dat_fig3), "rows\n")

cor_vars <- c("avg_c_k2", "avg_c_f", "avg_c_k2_f")
cor_results <- map_dfr(cor_vars, function(v) {
  ct <- cor.test(dat_fig3$acr_u, dat_fig3[[v]], method = "spearman", use = "complete.obs")
  tibble(var = v, rho = ct$estimate, pval = ct$p.value)
}) %>% mutate(
  label = formatC(rho, format = "f", digits = 2),
  sig = case_when(pval < 0.01 ~ "**", pval < 0.05 ~ "*", TRUE ~ ""),
  var_label = factor(case_when(
    var == "avg_c_k2"   ~ "Cortical k2",
    var == "avg_c_f"    ~ "Cortical F",
    var == "avg_c_k2_f" ~ "k2/F"
  ), levels = c("Cortical k2", "Cortical F", "k2/F"))
)

fig3a <- ggplot(cor_results, aes(x = "UACR", y = var_label, fill = rho)) +
  geom_tile(color = "white", linewidth = 1) +
  geom_text(aes(label = paste0(label, sig)), size = 4, fontface = "bold") +
  scale_fill_gradient2(low = cols$control, mid = "white", high = cols$t2d,
                       midpoint = 0, limits = c(-0.6, 0.6), name = "Spearman\nrho") +
  labs(x = NULL, y = NULL, tag = "A") +
  theme_rockies + theme(legend.position = "right")

fig3b <- plot_corr(dat_fig3, "acr_u", "avg_c_k2", "UACR (mg/g)",
                   expression(bold("Cortical k"[2]*" (min"^{-1}*")")), "B")
fig3c <- plot_corr(dat_fig3, "acr_u", "avg_c_f", "UACR (mg/g)",
                   expression(bold("Cortical F (ml/g/min)")), "C")
fig3d <- plot_corr(dat_fig3, "acr_u", "avg_c_k2_f", "UACR (mg/g)",
                   expression(bold("Cortical k"[2]*"/F")), "D")

fig3 <- (fig3a | fig3b) / (fig3c | fig3d) +
  plot_annotation(
    title = "Figure 3. Kidney Oxidative Metabolism Correlates with Albuminuria",
    theme = theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5))
  )

save_fig(fig3, "Figure3_UACR_Correlations", width = 11, height = 9)
cat("Figure 3 saved!\n")

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
#  FIGURES 6-7: scRNAseq VOLCANO PLOTS (need NEBULA results)
# =====================================================================
########################################################################

plot_volcano <- function(de_results, gene_list, title, tag, padj_thresh = 0.05) {
  df <- de_results %>%
    filter(gene %in% gene_list) %>%
    mutate(
      sig = case_when(
        padj < padj_thresh & log2FC > 0  ~ "Up",
        padj < padj_thresh & log2FC < 0  ~ "Down",
        TRUE ~ "NS"),
      label = ifelse(padj < padj_thresh, gene, NA_character_))
  
  ggplot(df, aes(log2FC, -log10(pvalue))) +
    geom_point(aes(color = sig), size = 2, alpha = 0.7) +
    ggrepel::geom_text_repel(aes(label = label), size = 2.5, max.overlaps = 15) +
    geom_hline(yintercept = -log10(padj_thresh), linetype = "dashed", color = "gray50") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    scale_color_manual(values = c("Up" = cols$up, "Down" = cols$down, "NS" = cols$ns)) +
    labs(x = expression(bold("log"[2]*" FC")),
         y = expression(bold("-log"[10]*" p")), title = title, tag = tag) +
    theme_rockies + theme(legend.position = "bottom")
}

tca_genes <- c("CS","ACO1","ACO2","IDH1","IDH2","IDH3A","IDH3B","IDH3G",
               "OGDH","OGDHL","DLST","DLD","SUCLG1","SUCLG2","SUCLA2",
               "SDHA","SDHB","SDHC","SDHD","FH","MDH1","MDH2","PCK1",
               "PCK2","PC","ACLY","DLAT","PDHA1","PDHB")

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

# Uncomment when NEBULA results available:
# nebula_pt_all  <- read.csv("path/to/nebula_total_PT.csv")
# nebula_pt_s1s2 <- read.csv("path/to/nebula_PT_S1S2.csv")
# nebula_pt_s3   <- read.csv("path/to/nebula_PT_S3.csv")
# nebula_apt     <- read.csv("path/to/nebula_aPT.csv")
#
# fig6 <- (plot_volcano(nebula_pt_all, tca_genes, "Total PT", "A") |
#           plot_volcano(nebula_pt_s1s2, tca_genes, "PT-S1/S2", "B")) /
#          (plot_volcano(nebula_pt_s3, tca_genes, "PT-S3", "C") |
#           plot_volcano(nebula_apt, tca_genes, "aPT", "D")) +
#   plot_annotation(title = "Figure 6. TCA Cycle Genes in PT Cells")
# save_fig(fig6, "Figure6_TCA_Volcanos", width = 12, height = 10)
#
# fig7 <- (plot_volcano(nebula_pt_all, oxphos_genes, "Total PT", "A") |
#           plot_volcano(nebula_pt_s1s2, oxphos_genes, "PT-S1/S2", "B")) /
#          (plot_volcano(nebula_pt_s3, oxphos_genes, "PT-S3", "C") |
#           plot_volcano(nebula_apt, oxphos_genes, "aPT", "D")) +
#   plot_annotation(title = "Figure 7. OxPhos Genes in PT Cells")
# save_fig(fig7, "Figure7_OxPhos_Volcanos", width = 12, height = 10)

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
cat("  FIGURE 6 (TCA Volcanos):       --> NEED NEBULA results\n")
cat("  FIGURE 7 (OxPhos Volcanos):    --> NEED NEBULA results\n")
cat("================================================================\n")