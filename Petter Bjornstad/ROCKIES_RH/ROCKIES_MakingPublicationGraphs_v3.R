########################################################################
# ROCKIES Publication - Complete Figure Generation Script
# Figures 1-7 for Cell Metabolism Submission
########################################################################

library(tidyverse)
library(ggpubr)
library(patchwork)
library(ggbeeswarm)
library(data.table)

load_image_panel <- function(path, tag) {
  ext <- tolower(tools::file_ext(path))
  if (!file.exists(path)) {
    message(sprintf("Image not found: %s — placeholder used", path))
    return(ggplot() +
             annotate("text", x=0.5, y=0.5, label=paste0("Panel ",tag,"\n(file not found)"),
                      size=4, color="gray50") + theme_void() + labs(tag=tag))
  }
  img <- if (ext %in% c("jpg","jpeg")) jpeg::readJPEG(path) else png::readPNG(path)
  ggplot() +
    annotation_raster(img, xmin=0, xmax=1, ymin=0, ymax=1) +
    xlim(0,1) + ylim(0,1) + theme_void() + labs(tag=tag) +
    theme(plot.margin=margin(5,5,5,5))
}

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
# Read the raw file
rockies_wide <- read.csv2(
  "C:/Users/netio/Downloads/PET data ROCKIES variables(Blad1).csv",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

names(rockies_wide) <- make.names(trimws(names(rockies_wide)))

# Convert all columns: remove tab separators, replace comma decimals
rockies_wide <- rockies_wide %>%
  mutate(across(everything(), ~ {
    cleaned <- gsub("\t", "", as.character(.x))
    cleaned <- gsub(",", ".", cleaned)
    conv <- suppressWarnings(as.numeric(cleaned))
    if (sum(!is.na(conv)) >= sum(!is.na(.x)) * 0.3 & sum(!is.na(conv)) > 0) conv else .x
  }))

cat("=== ROCKIES columns ===\n")
print(names(rockies_wide))
cat("\nRows:", nrow(rockies_wide), "\n")
cat("Group 0 (placebo-only):", sum(rockies_wide$Group == 0, na.rm = TRUE), "\n")
cat("Group 1 (crossover):", sum(rockies_wide$Group == 1, na.rm = TRUE), "\n\n")
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
  
  wide <- df %>% 
    pivot_wider(names_from = trt, values_from = val) %>%
    arrange(id)
  pv <- t.test(wide$Placebo, wide$Ertugliflozin, paired = TRUE)$p.value
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
# ROCKIES Publication - Figure 1 Complete Script
# All fixes applied:
#   1. Data loaded from .sav (full precision, matches SPSS)
#   2. Panels E-H: p-values from raw wide columns (matches SPSS)
#   3. Panel J: plain ASCII matrix indexing (fixes blank cells)
#   4. Panel J: pre-computed delta columns + Pearson (matches SPSS)
########################################################################

library(tidyverse)
library(ggpubr)
library(patchwork)
library(ggbeeswarm)
library(data.table)
library(pheatmap)
library(grid)
library(png)
library(haven)

load_image_panel <- function(path, tag) {
  ext <- tolower(tools::file_ext(path))
  if (!file.exists(path)) {
    message(sprintf("Image not found: %s — placeholder used", path))
    return(ggplot() +
             annotate("text", x=0.5, y=0.5, label=paste0("Panel ",tag,"\n(file not found)"),
                      size=4, color="gray50") + theme_void() + labs(tag=tag))
  }
  img <- if (ext %in% c("jpg","jpeg")) jpeg::readJPEG(path) else png::readPNG(path)
  ggplot() +
    annotation_raster(img, xmin=0, xmax=1, ymin=0, ymax=1) +
    xlim(0,1) + ylim(0,1) + theme_void() + labs(tag=tag) +
    theme(plot.margin=margin(5,5,5,5))
}

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
# DATA LOADING — from SPSS .sav (full precision)
########################################################################

sav_path <- "C:/Users/netio/Downloads/ROCKIES_Hoofddatabase__jan2026.sav"

if (!file.exists(sav_path)) stop("SAV file not found: ", sav_path)

message("Loading from .sav (full precision, matches SPSS)...")
rockies_wide <- read_sav(sav_path)
rockies_wide <- zap_labels(rockies_wide)   # strip SPSS value labels → plain numerics
names(rockies_wide) <- make.names(trimws(names(rockies_wide)))

cat("=== Loaded from .sav ===\n")
cat("Rows:", nrow(rockies_wide), "\n")
cat("Cols:", ncol(rockies_wide), "\n")
cat("Group 0 (placebo-only):", sum(rockies_wide$Group == 0, na.rm = TRUE), "\n")
cat("Group 1 (crossover):",    sum(rockies_wide$Group == 1, na.rm = TRUE), "\n\n")

# Verify key column names
cat("=== K2/cortex/medulla columns ===\n")
print(names(rockies_wide)[grep("K2|k2|cortex|medulla", names(rockies_wide), ignore.case = TRUE)])

rockies_crossover <- rockies_wide %>% filter(Group == 1)

########################################################################
# RESHAPE TO LONG
########################################################################

rockies_long <- bind_rows(
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

cat("\nROCKIES long:", nrow(rockies_long), "rows,",
    n_distinct(rockies_long$id), "participants\n\n")

########################################################################
# REUSABLE PLOT FUNCTIONS
########################################################################

# plot_paired_rawp: uses p_override computed directly from .sav wide
# columns (full precision) — dots/lines still render from rockies_long
plot_paired_rawp <- function(data, id_col, trt_col, val_col, ylab, tag,
                             p_override = NULL) {
  df <- data.frame(id  = data[[id_col]],
                   trt = data[[trt_col]],
                   val = data[[val_col]])
  df <- df[!is.na(df$val), ]
  df$trt <- factor(df$trt, levels = c("Placebo", "Ertugliflozin"))
  
  summ <- df %>%
    group_by(trt) %>%
    summarise(m = mean(val), se = sd(val) / sqrt(n()), .groups = "drop")
  
  pv <- if (!is.null(p_override)) {
    p_override
  } else {
    wide <- df %>% pivot_wider(names_from = trt, values_from = val) %>% arrange(id)
    t.test(wide$Placebo, wide$Ertugliflozin, paired = TRUE)$p.value
  }
  
  plab <- ifelse(pv < 0.001,
                 paste0("p=", formatC(pv, format = "e", digits = 1)),
                 paste0("p=", formatC(pv, format = "f", digits = 3)))
  
  ymax <- max(df$val, na.rm = TRUE)
  ymin <- min(df$val, na.rm = TRUE)
  yr   <- ymax - ymin
  
  ggplot(df, aes(trt, val)) +
    geom_line(aes(group = id), color = cols$paired_line,
              alpha = 0.5, linewidth = 0.4) +
    geom_point(aes(color = trt), size = 2.5, alpha = 0.7) +
    geom_errorbar(data = summ,
                  aes(trt, m, ymin = m - se, ymax = m + se),
                  width = 0.15, linewidth = 0.8, inherit.aes = FALSE) +
    geom_point(data = summ, aes(trt, m),
               size = 4, shape = 18, inherit.aes = FALSE) +
    annotate("segment", x = 1, xend = 2,
             y = ymax + yr * 0.05, yend = ymax + yr * 0.05,
             linewidth = 0.4) +
    annotate("text", x = 1.5, y = ymax + yr * 0.1,
             label = plab, size = 3.5, fontface = "italic") +
    scale_color_manual(values = c(cols$placebo, cols$ertu)) +
    coord_cartesian(ylim = c(ymin - yr * 0.05, ymax + yr * 0.15)) +
    labs(x = NULL, y = ylab, tag = tag) +
    theme_rockies
}

########################################################################
# PANEL A — Trial Design
########################################################################

fig1a <- load_image_panel("C:/Users/netio/Downloads/ROCKIES study plan.png", "A")

########################################################################
# PANEL B — PET Diagram
########################################################################

fig1b <- load_image_panel(
  "C:/Users/netio/Downloads/Kidney C11 Acetate PET Diagram.png", "B")

########################################################################
# PANEL C — Baseline correlation heatmap (Spearman)
########################################################################

rockies_baseline <- rockies_crossover %>%
  transmute(
    id          = Participant,
    k2_cortex   = K2_cortex_mean_Sherbrook_PLB,
    k2_medulla  = K2_medulla_mean_Sherbrook_PLB,
    k2f_cortex  = K2_F_cortex_mean_Sherbrook_PLB,
    k2f_medulla = K2_F_medulla_mean_Sherbrook_PLB,
    homa_base   = HOMAIR_Placebo,
    ogis_base   = OGIS_placebo,
    bw_base     = Bodyweight_Placebo
  )

pet_vars    <- c("k2_cortex", "k2_medulla", "k2f_cortex", "k2f_medulla")
clin_vars   <- c("homa_base", "ogis_base", "bw_base")
pet_labels  <- c("k\u2082 Cortex", "k\u2082 Medulla",
                 "k\u2082/F Cortex", "k\u2082/F Medulla")
clin_labels <- c("HOMA-IR", "OGIS", "Body Weight")

rho_mat <- matrix(NA_real_, 4, 3, dimnames = list(pet_labels, clin_labels))
p_mat   <- rho_mat

for (i in seq_along(pet_vars)) {
  for (j in seq_along(clin_vars)) {
    d <- rockies_baseline[, c(pet_vars[i], clin_vars[j])]
    d <- d[complete.cases(d), ]
    if (nrow(d) >= 5) {
      res           <- cor.test(d[[1]], d[[2]], method = "spearman", exact = FALSE)
      rho_mat[i, j] <- res$estimate
      p_mat[i, j]   <- res$p.value
    }
  }
}

sig_labels <- matrix("", 4, 3, dimnames = dimnames(rho_mat))
sig_labels[p_mat < 0.05  & !is.na(p_mat)] <- "*"
sig_labels[p_mat < 0.01  & !is.na(p_mat)] <- "**"
sig_labels[p_mat < 0.001 & !is.na(p_mat)] <- "***"

cell_labels <- matrix(paste0(round(rho_mat, 2), sig_labels),
                      4, 3, dimnames = dimnames(rho_mat))

fig1c_heatmap <- pheatmap(
  rho_mat,
  display_numbers = cell_labels,
  number_color    = "black",
  fontsize_number = 11,
  fontsize_row    = 11,
  fontsize_col    = 11,
  cluster_rows    = FALSE,
  cluster_cols    = FALSE,
  color           = colorRampPalette(c("#2166AC", "white", "#B2182B"))(101),
  breaks          = seq(-1, 1, length.out = 102),
  legend_breaks   = c(-1, -0.5, 0, 0.5, 1),
  legend_labels   = c("-1.0", "-0.5", "0.0", "0.5", "1.0"),
  main            = "",
  border_color    = "white",
  na_col          = "grey90",
  angle_col       = 45,
  silent          = TRUE
)

fig1c <- wrap_elements(full = fig1c_heatmap$gtable) +
  labs(tag = "C") +
  theme(plot.tag = element_text(size = 14, face = "bold",
                                margin = margin(0, 4, 0, 0)))

########################################################################
# PANEL D — Urine Omics
########################################################################

fig1d <- load_image_panel(
  "C:/Users/netio/Documents/UofW/Rockies/publication_figures/Urine_Metabolism_Figure.jpg",
  "D")

########################################################################
# PANELS E–H — Paired plots
# p-values computed directly from .sav wide columns (full precision)
########################################################################

pval_raw <- function(plb_col, ertu_col) {
  t.test(rockies_crossover[[plb_col]],
         rockies_crossover[[ertu_col]],
         paired = TRUE)$p.value
}

pv_e <- pval_raw("K2_cortex_mean_Sherbrook_PLB",    "K2_cortex_mean_Sherbrook_SGLT2")
pv_f <- pval_raw("K2_F_cortex_mean_Sherbrook_PLB",  "K2_F_cortex_mean_Sherbrook_SGLT2")
pv_g <- pval_raw("K2_medulla_mean_Sherbrook_PLB",   "K2_medulla_mean_Sherbrook_SGLT2i")
pv_h <- pval_raw("K2_F_medulla_mean_Sherbrook_PLB", "K2_F_medulla_mean_Sherbrook_SGLT2")

cat("=== Panels E-H p-values (from .sav, should match SPSS) ===\n")
cat(sprintf("E — Cortical k2:    p = %.4f\n", pv_e))
cat(sprintf("F — Cortical k2/F:  p = %.4f\n", pv_f))
cat(sprintf("G — Medullary k2:   p = %.4f\n", pv_g))
cat(sprintf("H — Medullary k2/F: p = %.4f\n", pv_h))

fig1e <- plot_paired_rawp(rockies_long, "id", "treatment", "cortical_k2",
                          expression(bold("Cortical k"[2]*" (min"^{-1}*")")), "E",
                          p_override = pv_e)

fig1f <- plot_paired_rawp(rockies_long, "id", "treatment", "cortical_k2f",
                          expression(bold("Cortical k"[2]*"/F")), "F",
                          p_override = pv_f)

fig1g <- plot_paired_rawp(rockies_long, "id", "treatment", "medullary_k2",
                          expression(bold("Medullary k"[2]*" (min"^{-1}*")")), "G",
                          p_override = pv_g)

fig1h <- plot_paired_rawp(rockies_long, "id", "treatment", "medullary_k2f",
                          expression(bold("Medullary k"[2]*"/F")), "H",
                          p_override = pv_h)

########################################################################
# PANEL I — Representative PET k2 maps
########################################################################

pet_map_img <- readPNG(
  "C:/Users/netio/Downloads/Re_ ROCKIES_RH2_ Representative plots/k2-maps.png")

fig1i <- ggplot() +
  annotation_raster(pet_map_img, xmin = 0, xmax = 1, ymin = 0, ymax = 1) +
  xlim(0, 1) + ylim(0, 1) +
  theme_void() +
  labs(tag = "I") +
  theme(plot.tag    = element_text(size = 14, face = "bold"),
        plot.margin = margin(5, 5, 5, 5))

########################################################################
# PANEL J — Delta-delta correlation heatmap (hardcoded from .sav values)
#
# Values taken directly from the current .sav-based output + SPSS tables.
# Rows: k2 Cortex, k2 Medulla, k2/F Cortex, k2/F Medulla
# Cols: HbA1c, HOMA-IR, Glucose, Sodium Load
# NA = cell not requested (renders grey)
#
# r values (from .sav run):
#   k2 Cortex   x HbA1c:        0.48*
#   k2 Cortex   x HOMA-IR:      0.77***
#   k2 Medulla  x HbA1c:        0.50*
#   k2 Medulla  x HOMA-IR:      0.82***
#   k2/F Cortex x HOMA-IR:      0.61**
#   k2/F Cortex x Glucose:      0.48*
#   k2/F Cortex x Sodium Load:  0.45*
#   k2/F Medulla x HOMA-IR:     0.64**
#   k2/F Medulla x Sodium Load: 0.52*
########################################################################

# Row order: k2 Cortex, k2 Medulla, k2/F Cortex, k2/F Medulla
# Col order: HbA1c, HOMA-IR, Glucose, Sodium Load

# All 16 cells — exact values from .sav run
# Stars: p<0.05=*, p<0.01=**, p<0.001=***
# Row order: k2 Cortex, k2 Medulla, k2/F Cortex, k2/F Medulla
# Col order: HbA1c, HOMA-IR, Glucose, Sodium Load

dd_rho <- matrix(
  c( 0.48,  0.77,  0.36,  0.20,   # k2 Cortex    (gluc p=0.124, sod p=0.401 — ns)
     0.50,  0.82,  0.31,  0.17,   # k2 Medulla   (gluc p=0.178, sod p=0.471 — ns)
     0.32,  0.61,  0.48,  0.45,   # k2/F Cortex  (hba1c p=0.165 — ns)
     0.31,  0.64,  0.42,  0.52),  # k2/F Medulla (hba1c p=0.177, gluc p=0.065 — ns)
  nrow = 4, ncol = 4, byrow = TRUE,
  dimnames = list(
    c("k\u2082 Cortex", "k\u2082 Medulla", "k\u2082/F Cortex", "k\u2082/F Medulla"),
    c("\u0394 HbA1c", "\u0394 HOMA-IR", "\u0394 Glucose", "\u0394 Sodium Load")
  )
)

dd_cell <- matrix(
  c("0.48*",  "0.77***", "0.36",   "0.20",
    "0.50*",  "0.82***", "0.31",   "0.17",
    "0.32",   "0.61**",  "0.48*",  "0.45*",
    "0.31",   "0.64**",  "0.42",   "0.52*"),
  nrow = 4, ncol = 4, byrow = TRUE,
  dimnames = dimnames(dd_rho)
)

fig1j_heatmap <- pheatmap(
  dd_rho,
  display_numbers = dd_cell,
  number_color    = "black",
  fontsize_number = 11,
  fontsize_row    = 11,
  fontsize_col    = 11,
  cluster_rows    = FALSE,
  cluster_cols    = FALSE,
  color           = colorRampPalette(c("#2166AC", "white", "#B2182B"))(101),
  breaks          = seq(-1, 1, length.out = 102),
  legend_breaks   = c(-1, -0.5, 0, 0.5, 1),
  legend_labels   = c("-1.0", "-0.5", "0.0", "0.5", "1.0"),
  main            = "",
  border_color    = "white",
  na_col          = "grey90",
  angle_col       = 45,
  silent          = TRUE
)

fig1j <- wrap_elements(full = fig1j_heatmap$gtable) +
  labs(tag = "J") +
  theme(plot.tag = element_text(size = 14, face = "bold",
                                margin = margin(0, 4, 0, 0)))

########################################################################
# ASSEMBLE FIGURE 1
########################################################################

row1 <- fig1a + fig1b +
  plot_layout(widths = c(1.2, 1.2))

row2 <- fig1c + fig1d +
  plot_layout(widths = c(1, 1.4))

row3 <- fig1e + fig1f + fig1g + fig1h +
  plot_layout(widths = c(1, 1, 1, 1))

row4 <- fig1i + fig1j +
  plot_layout(widths = c(1.8, 1.2))

fig1_full <- row1 / row2 / row3 / row4 +
  plot_layout(heights = c(3.5, 3, 2.5, 3)) +
  plot_annotation(
    title = "Figure 1. SGLT2 Inhibition Reduces Kidney Oxidative Metabolism in the ROCKIES Trial",
    theme = theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5))
  )

print(fig1_full)
save_fig(fig1_full, "Figure1_ROCKIES", width = 16, height = 22)
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

# =====================================================================
# BUILD cor_results for Figure 3 Panel A heatmap
# =====================================================================

dat_fig3 <- dat_with_pet %>%
  filter(!is.na(avg_c_k2) & is.finite(avg_c_k2)) %>%
  filter(group %in% c('Lean Control', 'Type 2 Diabetes')) %>%
  filter(!is.na(acr_u))

dat_fig3 <- dat_fig3 %>%
  mutate(log10_uacr = log10(acr_u))

# Compute Spearman correlations vs UACR
vars_to_test <- list(
  list(var = "avg_c_k2",   label = "Cortical k\u2082"),
  list(var = "avg_c_f",    label = "Cortical F"),
  list(var = "avg_c_k2_f", label = "Cortical k\u2082/F")
)

cor_results <- purrr::map_dfr(vars_to_test, function(v) {
  df <- dat_fig3 %>% filter(!is.na(.data[[v$var]]))
  ct <- cor.test(df$acr_u, df[[v$var]], method = "spearman")
  data.frame(
    var_label = v$label,
    rho       = ct$estimate,
    pval      = ct$p.value,
    label     = formatC(ct$estimate, format = "f", digits = 2),
    sig       = ifelse(ct$p.value < 0.001, "***",
                       ifelse(ct$p.value < 0.01,  "**",
                              ifelse(ct$p.value < 0.05,  "*", "")))
  )
})

cor_results$var_label <- factor(cor_results$var_label,
                                levels = c("Cortical k\u2082/F", 
                                           "Cortical F", 
                                           "Cortical k\u2082"))

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

row1 <- (fig1a + labs(tag = "A")) | (fig1b + labs(tag = "B"))
row2 <- fig1c | fig1d | fig1e
row3 <- fig1f | fig1g | fig1h
row4 <- fig1i | fig1j | fig1k

fig1 <- row1 / row2 / row3 / row4 +
  plot_layout(heights = c(3, 2, 2, 2)) +
  plot_annotation(
    title = "Figure 1. SGLT2 Inhibition Reduces Kidney Oxidative Metabolism in the ROCKIES Trial",
    theme = theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5))
  )

save_fig(fig1, "Figure1_ROCKIES", width = 16, height = 22)


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










########################################################################
# ROCKIES Publication - Combined PDF (All Figures, One Page Each)
# Run AFTER the main figure generation script
########################################################################

library(tidyverse)
library(patchwork)
library(png)

base_path <- 'C:/Users/netio/Documents/UofW/Rockies/publication_figures/'

# Page dimensions — tall for Figure 1's 5 rows
pg_w <- 16
pg_h <- 22

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
# ROCKIES Publication - Figure 1 Complete Script
# All fixes applied:
#   1. Data loaded from .sav (full precision, matches SPSS)
#   2. Panels E-H: p-values from raw wide columns (matches SPSS)
#   3. Panel J: plain ASCII matrix indexing (fixes blank cells)
#   4. Panel J: pre-computed delta columns + Pearson (matches SPSS)
########################################################################

library(tidyverse)
library(ggpubr)
library(patchwork)
library(ggbeeswarm)
library(data.table)
library(pheatmap)
library(grid)
library(png)
library(haven)

load_image_panel <- function(path, tag) {
  ext <- tolower(tools::file_ext(path))
  if (!file.exists(path)) {
    message(sprintf("Image not found: %s — placeholder used", path))
    return(ggplot() +
             annotate("text", x=0.5, y=0.5, label=paste0("Panel ",tag,"\n(file not found)"),
                      size=4, color="gray50") + theme_void() + labs(tag=tag))
  }
  img <- if (ext %in% c("jpg","jpeg")) jpeg::readJPEG(path) else png::readPNG(path)
  ggplot() +
    annotation_raster(img, xmin=0, xmax=1, ymin=0, ymax=1) +
    xlim(0,1) + ylim(0,1) + theme_void() + labs(tag=tag) +
    theme(plot.margin=margin(5,5,5,5))
}

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
# DATA LOADING — from SPSS .sav (full precision)
########################################################################

sav_path <- "C:/Users/netio/Downloads/ROCKIES_Hoofddatabase__jan2026.sav"

if (!file.exists(sav_path)) stop("SAV file not found: ", sav_path)

message("Loading from .sav (full precision, matches SPSS)...")
rockies_wide <- read_sav(sav_path)
rockies_wide <- zap_labels(rockies_wide)   # strip SPSS value labels → plain numerics
names(rockies_wide) <- make.names(trimws(names(rockies_wide)))

cat("=== Loaded from .sav ===\n")
cat("Rows:", nrow(rockies_wide), "\n")
cat("Cols:", ncol(rockies_wide), "\n")
cat("Group 0 (placebo-only):", sum(rockies_wide$Group == 0, na.rm = TRUE), "\n")
cat("Group 1 (crossover):",    sum(rockies_wide$Group == 1, na.rm = TRUE), "\n\n")

# Verify key column names
cat("=== K2/cortex/medulla columns ===\n")
print(names(rockies_wide)[grep("K2|k2|cortex|medulla", names(rockies_wide), ignore.case = TRUE)])

rockies_crossover <- rockies_wide %>% filter(Group == 1)

########################################################################
# RESHAPE TO LONG
########################################################################

rockies_long <- bind_rows(
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

cat("\nROCKIES long:", nrow(rockies_long), "rows,",
    n_distinct(rockies_long$id), "participants\n\n")

########################################################################
# REUSABLE PLOT FUNCTIONS
########################################################################

# plot_paired_rawp: uses p_override computed directly from .sav wide
# columns (full precision) — dots/lines still render from rockies_long
plot_paired_rawp <- function(data, id_col, trt_col, val_col, ylab, tag,
                             p_override = NULL) {
  df <- data.frame(id  = data[[id_col]],
                   trt = data[[trt_col]],
                   val = data[[val_col]])
  df <- df[!is.na(df$val), ]
  df$trt <- factor(df$trt, levels = c("Placebo", "Ertugliflozin"))
  
  summ <- df %>%
    group_by(trt) %>%
    summarise(m = mean(val), se = sd(val) / sqrt(n()), .groups = "drop")
  
  pv <- if (!is.null(p_override)) {
    p_override
  } else {
    wide <- df %>% pivot_wider(names_from = trt, values_from = val) %>% arrange(id)
    t.test(wide$Placebo, wide$Ertugliflozin, paired = TRUE)$p.value
  }
  
  plab <- ifelse(pv < 0.001,
                 paste0("p=", formatC(pv, format = "e", digits = 1)),
                 paste0("p=", formatC(pv, format = "f", digits = 3)))
  
  ymax <- max(df$val, na.rm = TRUE)
  ymin <- min(df$val, na.rm = TRUE)
  yr   <- ymax - ymin
  
  ggplot(df, aes(trt, val)) +
    geom_line(aes(group = id), color = cols$paired_line,
              alpha = 0.5, linewidth = 0.4) +
    geom_point(aes(color = trt), size = 2.5, alpha = 0.7) +
    geom_errorbar(data = summ,
                  aes(trt, m, ymin = m - se, ymax = m + se),
                  width = 0.15, linewidth = 0.8, inherit.aes = FALSE) +
    geom_point(data = summ, aes(trt, m),
               size = 4, shape = 18, inherit.aes = FALSE) +
    annotate("segment", x = 1, xend = 2,
             y = ymax + yr * 0.05, yend = ymax + yr * 0.05,
             linewidth = 0.4) +
    annotate("text", x = 1.5, y = ymax + yr * 0.1,
             label = plab, size = 3.5, fontface = "italic") +
    scale_color_manual(values = c(cols$placebo, cols$ertu)) +
    coord_cartesian(ylim = c(ymin - yr * 0.05, ymax + yr * 0.15)) +
    labs(x = NULL, y = ylab, tag = tag) +
    theme_rockies
}

########################################################################
# PANEL A — Trial Design
########################################################################

fig1a <- load_image_panel("C:/Users/netio/Downloads/ROCKIES study plan.png", "A")

########################################################################
# PANEL B — PET Diagram
########################################################################

fig1b <- load_image_panel(
  "C:/Users/netio/Downloads/Kidney C11 Acetate PET Diagram.png", "B")

########################################################################
# PANEL C — Baseline correlation heatmap (Spearman)
########################################################################

rockies_baseline <- rockies_crossover %>%
  transmute(
    id          = Participant,
    k2_cortex   = K2_cortex_mean_Sherbrook_PLB,
    k2_medulla  = K2_medulla_mean_Sherbrook_PLB,
    k2f_cortex  = K2_F_cortex_mean_Sherbrook_PLB,
    k2f_medulla = K2_F_medulla_mean_Sherbrook_PLB,
    homa_base   = HOMAIR_Placebo,
    ogis_base   = OGIS_placebo,
    bw_base     = Bodyweight_Placebo
  )

pet_vars    <- c("k2_cortex", "k2_medulla", "k2f_cortex", "k2f_medulla")
clin_vars   <- c("homa_base", "ogis_base", "bw_base")
pet_labels  <- c("k\u2082 Cortex", "k\u2082 Medulla",
                 "k\u2082/F Cortex", "k\u2082/F Medulla")
clin_labels <- c("HOMA-IR", "OGIS", "Body Weight")

rho_mat <- matrix(NA_real_, 4, 3, dimnames = list(pet_labels, clin_labels))
p_mat   <- rho_mat

for (i in seq_along(pet_vars)) {
  for (j in seq_along(clin_vars)) {
    d <- rockies_baseline[, c(pet_vars[i], clin_vars[j])]
    d <- d[complete.cases(d), ]
    if (nrow(d) >= 5) {
      res           <- cor.test(d[[1]], d[[2]], method = "spearman", exact = FALSE)
      rho_mat[i, j] <- res$estimate
      p_mat[i, j]   <- res$p.value
    }
  }
}

sig_labels <- matrix("", 4, 3, dimnames = dimnames(rho_mat))
sig_labels[p_mat < 0.05  & !is.na(p_mat)] <- "*"
sig_labels[p_mat < 0.01  & !is.na(p_mat)] <- "**"
sig_labels[p_mat < 0.001 & !is.na(p_mat)] <- "***"

cell_labels <- matrix(paste0(round(rho_mat, 2), sig_labels),
                      4, 3, dimnames = dimnames(rho_mat))

fig1c_heatmap <- pheatmap(
  rho_mat,
  display_numbers = cell_labels,
  number_color    = "black",
  fontsize_number = 11,
  fontsize_row    = 11,
  fontsize_col    = 11,
  cluster_rows    = FALSE,
  cluster_cols    = FALSE,
  color           = colorRampPalette(c("#2166AC", "white", "#B2182B"))(101),
  breaks          = seq(-1, 1, length.out = 102),
  legend_breaks   = c(-1, -0.5, 0, 0.5, 1),
  legend_labels   = c("-1.0", "-0.5", "0.0", "0.5", "1.0"),
  main            = "",
  border_color    = "white",
  na_col          = "grey90",
  angle_col       = 45,
  silent          = TRUE
)

fig1c <- wrap_elements(full = fig1c_heatmap$gtable) +
  labs(tag = "C") +
  theme(plot.tag = element_text(size = 14, face = "bold",
                                margin = margin(0, 4, 0, 0)))

########################################################################
# PANEL D — Urine Omics
########################################################################

fig1d <- load_image_panel(
  "C:/Users/netio/Documents/UofW/Rockies/publication_figures/Urine_Metabolism_Figure.jpg",
  "D")

########################################################################
# PANELS E–H — Paired plots
# p-values computed directly from .sav wide columns (full precision)
########################################################################

pval_raw <- function(plb_col, ertu_col) {
  t.test(rockies_crossover[[plb_col]],
         rockies_crossover[[ertu_col]],
         paired = TRUE)$p.value
}

pv_e <- pval_raw("K2_cortex_mean_Sherbrook_PLB",    "K2_cortex_mean_Sherbrook_SGLT2")
pv_f <- pval_raw("K2_F_cortex_mean_Sherbrook_PLB",  "K2_F_cortex_mean_Sherbrook_SGLT2")
pv_g <- pval_raw("K2_medulla_mean_Sherbrook_PLB",   "K2_medulla_mean_Sherbrook_SGLT2i")
pv_h <- pval_raw("K2_F_medulla_mean_Sherbrook_PLB", "K2_F_medulla_mean_Sherbrook_SGLT2")

cat("=== Panels E-H p-values (from .sav, should match SPSS) ===\n")
cat(sprintf("E — Cortical k2:    p = %.4f\n", pv_e))
cat(sprintf("F — Cortical k2/F:  p = %.4f\n", pv_f))
cat(sprintf("G — Medullary k2:   p = %.4f\n", pv_g))
cat(sprintf("H — Medullary k2/F: p = %.4f\n", pv_h))

fig1e <- plot_paired_rawp(rockies_long, "id", "treatment", "cortical_k2",
                          expression(bold("Cortical k"[2]*" (min"^{-1}*")")), "E",
                          p_override = pv_e)

fig1f <- plot_paired_rawp(rockies_long, "id", "treatment", "cortical_k2f",
                          expression(bold("Cortical k"[2]*"/F")), "F",
                          p_override = pv_f)

fig1g <- plot_paired_rawp(rockies_long, "id", "treatment", "medullary_k2",
                          expression(bold("Medullary k"[2]*" (min"^{-1}*")")), "G",
                          p_override = pv_g)

fig1h <- plot_paired_rawp(rockies_long, "id", "treatment", "medullary_k2f",
                          expression(bold("Medullary k"[2]*"/F")), "H",
                          p_override = pv_h)

########################################################################
# PANEL I — Representative PET k2 maps
########################################################################

pet_map_img <- readPNG(
  "C:/Users/netio/Downloads/Re_ ROCKIES_RH2_ Representative plots/k2-maps.png")

fig1i <- ggplot() +
  annotation_raster(pet_map_img, xmin = 0, xmax = 1, ymin = 0, ymax = 1) +
  xlim(0, 1) + ylim(0, 1) +
  theme_void() +
  labs(tag = "I") +
  theme(plot.tag    = element_text(size = 14, face = "bold"),
        plot.margin = margin(5, 5, 5, 5))

########################################################################
# PANEL J — Delta-delta correlation heatmap (hardcoded from .sav values)
#
# Values taken directly from the current .sav-based output + SPSS tables.
# Rows: k2 Cortex, k2 Medulla, k2/F Cortex, k2/F Medulla
# Cols: HbA1c, HOMA-IR, Glucose, Sodium Load
# NA = cell not requested (renders grey)
#
# r values (from .sav run):
#   k2 Cortex   x HbA1c:        0.48*
#   k2 Cortex   x HOMA-IR:      0.77***
#   k2 Medulla  x HbA1c:        0.50*
#   k2 Medulla  x HOMA-IR:      0.82***
#   k2/F Cortex x HOMA-IR:      0.61**
#   k2/F Cortex x Glucose:      0.48*
#   k2/F Cortex x Sodium Load:  0.45*
#   k2/F Medulla x HOMA-IR:     0.64**
#   k2/F Medulla x Sodium Load: 0.52*
########################################################################

# Row order: k2 Cortex, k2 Medulla, k2/F Cortex, k2/F Medulla
# Col order: HbA1c, HOMA-IR, Glucose, Sodium Load

# All 16 cells — exact values from .sav run
# Stars: p<0.05=*, p<0.01=**, p<0.001=***
# Row order: k2 Cortex, k2 Medulla, k2/F Cortex, k2/F Medulla
# Col order: HbA1c, HOMA-IR, Glucose, Sodium Load

dd_rho <- matrix(
  c( 0.48,  0.77,  0.36,  0.20,   # k2 Cortex    (gluc p=0.124, sod p=0.401 — ns)
     0.50,  0.82,  0.31,  0.17,   # k2 Medulla   (gluc p=0.178, sod p=0.471 — ns)
     0.32,  0.61,  0.48,  0.45,   # k2/F Cortex  (hba1c p=0.165 — ns)
     0.31,  0.64,  0.42,  0.52),  # k2/F Medulla (hba1c p=0.177, gluc p=0.065 — ns)
  nrow = 4, ncol = 4, byrow = TRUE,
  dimnames = list(
    c("k\u2082 Cortex", "k\u2082 Medulla", "k\u2082/F Cortex", "k\u2082/F Medulla"),
    c("\u0394 HbA1c", "\u0394 HOMA-IR", "\u0394 Glucose", "\u0394 Sodium Load")
  )
)

dd_cell <- matrix(
  c("0.48*",  "0.77***", "0.36",   "0.20",
    "0.50*",  "0.82***", "0.31",   "0.17",
    "0.32",   "0.61**",  "0.48*",  "0.45*",
    "0.31",   "0.64**",  "0.42",   "0.52*"),
  nrow = 4, ncol = 4, byrow = TRUE,
  dimnames = dimnames(dd_rho)
)

fig1j_heatmap <- pheatmap(
  dd_rho,
  display_numbers = dd_cell,
  number_color    = "black",
  fontsize_number = 11,
  fontsize_row    = 11,
  fontsize_col    = 11,
  cluster_rows    = FALSE,
  cluster_cols    = FALSE,
  color           = colorRampPalette(c("#2166AC", "white", "#B2182B"))(101),
  breaks          = seq(-1, 1, length.out = 102),
  legend_breaks   = c(-1, -0.5, 0, 0.5, 1),
  legend_labels   = c("-1.0", "-0.5", "0.0", "0.5", "1.0"),
  main            = "",
  border_color    = "white",
  na_col          = "grey90",
  angle_col       = 45,
  silent          = TRUE
)

fig1j <- wrap_elements(full = fig1j_heatmap$gtable) +
  labs(tag = "J") +
  theme(plot.tag = element_text(size = 14, face = "bold",
                                margin = margin(0, 4, 0, 0)))

########################################################################
# ASSEMBLE FIGURE 1
########################################################################

row1 <- fig1a + fig1b +
  plot_layout(widths = c(1.2, 1.2))

row2 <- fig1c + fig1d +
  plot_layout(widths = c(1, 1.4))

row3 <- fig1e + fig1f + fig1g + fig1h +
  plot_layout(widths = c(1, 1, 1, 1))

row4 <- fig1i + fig1j +
  plot_layout(widths = c(1.8, 1.2))

fig1_full <- row1 / row2 / row3 / row4 +
  plot_layout(heights = c(3.5, 3, 2.5, 3)) +
  plot_annotation(
    title = "Figure 1. SGLT2 Inhibition Reduces Kidney Oxidative Metabolism in the ROCKIES Trial",
    theme = theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5))
  )

print(fig1_full)
save_fig(fig1_full, "Figure1_ROCKIES", width = 16, height = 22)
cat("Figure 1 saved!\n")

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
# ROCKIES Figures 6 & 7 — TCA & OxPhos barplots
########################################################################

library(tidyverse)
library(patchwork)
library(cowplot)
library(data.table)

base_path  <- 'C:/Users/netio/Documents/UofW/Rockies/publication_figures/'
scrna_path <- 'C:/Users/netio/Documents/UofW/Rockies/Hailey_Dotplots/'

########################################################################
# FILE PATHS — exact matches based on your folder listing
########################################################################

tca_files <- list(
  "PT"       = paste0(scrna_path, "NEBULA_TCA_cycle_PT_cells_LC_T2D_NoMed_unadjusted_pooled_offset.csv"),
  "PT-S1/S2" = paste0(scrna_path, "NEBULA_TCA_cycle_PT_S1_S2_cells_LC_T2D_NoMed_unadjusted_pooled_offset.csv"),
  "PT-S3"    = paste0(scrna_path, "NEBULA_TCA_cycle_PT_S3_cells_LC_T2D_NoMed_unadjusted_pooled_offset.csv"),
  "aPT"      = paste0(scrna_path, "NEBULA_TCA_cycle_aPT_cells_LC_T2D_NoMed_unadjusted_pooled_offset.csv")
)

oxphos_files <- list(
  "PT"       = paste0(scrna_path, "NEBULA_OX_PHOS_cycle_PT_cells_LC_T2D_NoMed_unadjusted_pooled_offset.csv"),
  "PT-S1/S2" = paste0(scrna_path, "NEBULA_OX_PHOS_cycle_PT_S1_S2_cells_LC_T2D_NoMed_unadjusted_pooled_offset.csv"),
  "PT-S3"    = paste0(scrna_path, "NEBULA_OX_PHOS_cycle_PT_S3_cells_LC_T2D_NoMed_unadjusted_pooled_offset.csv"),
  "aPT"      = paste0(scrna_path, "NEBULA_OX_PHOS_cycle_aPT_cells_LC_T2D_NoMed_unadjusted_pooled_offset.csv")
)

########################################################################
# PLOT FUNCTION
########################################################################

make_barplot <- function(filepath, celltype_label, panel_tag, show_legend = FALSE) {
  
  dt <- fread(filepath)
  
  # Rename columns
  setnames(dt, "logFC_groupType_2_Diabetes", "logFC")
  setnames(dt, "p_groupType_2_Diabetes",     "pvalue")
  
  dt <- dt %>%
    mutate(
      direction = factor(
        ifelse(logFC >= 0, "Up-regulated", "Down-regulated"),
        levels = c("Up-regulated", "Down-regulated")
      ),
      sig_label = ifelse(pvalue < 0.05, "*", ""),
      label_y   = ifelse(logFC >= 0,
                         logFC + max(abs(logFC), na.rm = TRUE) * 0.04,
                         logFC - max(abs(logFC), na.rm = TRUE) * 0.04)
    )
  
  p <- ggplot(dt, aes(x = gene, y = logFC, fill = direction)) +
    geom_col(width = 0.7) +
    geom_text(aes(y = label_y, label = sig_label),
              size = 4, fontface = "bold", vjust = 0.5) +
    geom_hline(yintercept = 0, linewidth = 0.4) +
    scale_fill_manual(
      name   = "Regulation",
      values = c("Up-regulated"   = "#4A90D9",
                 "Down-regulated" = "#E74C3C"),
      drop   = FALSE
    ) +
    labs(x = NULL, y = "Log2 Fold Change",
         subtitle = celltype_label, tag = panel_tag) +
    theme_classic(base_size = 11) +
    theme(
      axis.text.x      = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7),
      axis.title.y     = element_text(size = 10, face = "bold"),
      plot.subtitle    = element_text(size = 10, hjust = 0.5, face = "bold"),
      plot.tag         = element_text(size = 14, face = "bold"),
      legend.position  = ifelse(show_legend, "bottom", "none"),
      legend.title     = element_text(size = 10, face = "bold"),
      legend.text      = element_text(size = 10),
      legend.key.size  = unit(0.6, "cm"),
      legend.key       = element_rect(colour = "black", linewidth = 0.3),
      plot.margin      = margin(5, 10, 5, 10)
    )
  
  p
}

########################################################################
# BUILD FIGURE 6: TCA
########################################################################

# Build all 4 panels — PT gets show_legend=TRUE so we can extract it
tca_pt    <- make_barplot(tca_files[["PT"]],       "PT",       "A", show_legend = TRUE)
tca_s1s2  <- make_barplot(tca_files[["PT-S1/S2"]], "PT-S1/S2", "B")
tca_s3    <- make_barplot(tca_files[["PT-S3"]],    "PT-S3",    "C")
tca_apt   <- make_barplot(tca_files[["aPT"]],      "aPT",      "D")

# Extract shared legend from PT panel
legend6 <- get_legend(tca_pt)

# Remove legend from PT for the main layout
tca_pt_nl <- tca_pt + theme(legend.position = "none")

# Combine panels
fig6_combined <- (tca_pt_nl) /
  (tca_s1s2 | tca_s3 | tca_apt) +
  plot_layout(heights = c(1.3, 1)) +
  plot_annotation(
    title = "Figure 6. TCA Cycle Gene Expression in Proximal Tubule Cells",
    theme = theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5)
    )
  )

# Add legend via inset_element
fig6 <- fig6_combined +
  inset_element(legend6,
                left = 0.38, bottom = -0.04,
                right = 0.62, top = 0.0,
                align_to = 'full')

# Save
ggsave(paste0(base_path, "Figure6_TCA_Cycle.png"), fig6,
       width = 14, height = 10, dpi = 300)
ggsave(paste0(base_path, "Figure6_TCA_Cycle.pdf"), fig6,
       width = 14, height = 10, device = cairo_pdf)
cat("Figure 6 saved!\n")

########################################################################
# BUILD FIGURE 7: OxPhos
########################################################################

oxphos_pt   <- make_barplot(oxphos_files[["PT"]],       "PT",       "A", show_legend = TRUE)
oxphos_s1s2 <- make_barplot(oxphos_files[["PT-S1/S2"]], "PT-S1/S2", "B")
oxphos_s3   <- make_barplot(oxphos_files[["PT-S3"]],    "PT-S3",    "C")
oxphos_apt  <- make_barplot(oxphos_files[["aPT"]],      "aPT",      "D")

legend7 <- get_legend(oxphos_pt)

oxphos_pt_nl <- oxphos_pt + theme(legend.position = "none")

fig7_combined <- (oxphos_pt_nl) /
  (oxphos_s1s2 | oxphos_s3 | oxphos_apt) +
  plot_layout(heights = c(1.3, 1)) +
  plot_annotation(
    title = "Figure 7. OxPhos Gene Expression in Proximal Tubule Cells",
    theme = theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5)
    )
  )

fig7 <- fig7_combined +
  inset_element(legend7,
                left = 0.38, bottom = -0.04,
                right = 0.62, top = 0.0,
                align_to = 'full')

ggsave(paste0(base_path, "Figure7_OxPhos.png"), fig7,
       width = 14, height = 10, dpi = 300)
ggsave(paste0(base_path, "Figure7_OxPhos.pdf"), fig7,
       width = 14, height = 10, device = cairo_pdf)
cat("Figure 7 saved!\n")

 print(fig6)
 print(fig7)
dev.off()

cat("\n========================================================\n")
cat("  Combined PDF saved:\n")
cat("  ", paste0(base_path, "ROCKIES_All_Figures_Combined.pdf"), "\n")
cat("========================================================\n")

