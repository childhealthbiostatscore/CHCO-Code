library(tidyverse)
library(ggpubr)
library(patchwork)
library(ggbeeswarm)
library(data.table)
library(ggplotify)
library(grid)
library(png)
library(haven)
library(cowplot)
library(qpdf)
library(ggtext)
library(conflicted)
conflicts_prefer(dplyr::last)
conflicts_prefer(ggpubr::get_legend)
conflicts_prefer(dplyr::select)
conflicts_prefer(dplyr::filter)

if (.Platform$OS.type == "windows") {
  tryCatch({
    library(extrafont)
    loadfonts(device = "pdf", quiet = TRUE)
    loadfonts(device = "win", quiet = TRUE)
  }, error = function(e) message("extrafont loading skipped"))
}

base_path  <- "C:/Users/netio/Documents/UofW/Rockies/publication_figures/"
scrna_path <- "C:/Users/netio/Documents/UofW/Rockies/Hailey_Dotplots/"

########################################################################
# THEME & COLORS
########################################################################

theme_rockies <- theme_classic(base_size = 11) +
  theme(
    text              = element_text(family = "sans"),
    axis.title        = element_text(size = 10, face = "bold"),
    axis.text         = element_text(size = 9, color = "black"),
    plot.title        = element_text(size = 11, face = "bold", hjust = 0.5),
    plot.tag          = element_text(size = 14, face = "bold"),
    plot.tag.position = "topleft",
    strip.text        = element_text(size = 10, face = "bold"),
    strip.background  = element_blank(),
    legend.position   = "none",
    plot.margin       = margin(5, 10, 5, 10)
  )

cols <- list(
  placebo     = "#3f78c1",
  ertu        = "#ce9a28",
  paired_line = "gray60",
  control     = "#00A99D",
  obese_ctrl  = "#FFC000",
  t2d         = "#C55A11",
  corr_point  = "#2C3E50",
  delta       = "#8E44AD",
  gbm_yes     = "#E74C3C", gbm_no     = "#4A90D9",
  arterio_yes = "#E74C3C", arterio_no = "#4A90D9",
  down        = "#E74C3C", up         = "#4A90D9", ns = "gray70"
)

# Markdown labels for heatmap y-axes (HTML subscript rendered by element_markdown)
pet_labels <- c("k2 Cortex", "k2 Medulla", "k2/F Cortex", "k2/F Medulla")

pet_md_labels <- c(
  "k2 Cortex"    = "k<sub>2</sub> Cortex",
  "k2 Medulla"   = "k<sub>2</sub> Medulla",
  "k2/F Cortex"  = "k<sub>2</sub>/F Cortex",
  "k2/F Medulla" = "k<sub>2</sub>/F Medulla"
)

pet_delta_md_labels <- c(
  "k2 Cortex"    = "\u0394 k<sub>2</sub> Cortex",
  "k2 Medulla"   = "\u0394 k<sub>2</sub> Medulla",
  "k2/F Cortex"  = "\u0394 k<sub>2</sub>/F Cortex",
  "k2/F Medulla" = "\u0394 k<sub>2</sub>/F Medulla"
)

########################################################################
# HELPERS
########################################################################

load_image_panel <- function(path, tag) {
  ext <- tolower(tools::file_ext(path))
  if (!file.exists(path)) {
    message(sprintf("Image not found: %s — placeholder used", path))
    return(ggplot() +
             annotate("text", x = 0.5, y = 0.5,
                      label = paste0("Panel ", tag, "\n(file not found)"),
                      size = 4, color = "gray50") +
             theme_void() + labs(tag = tag) +
             theme(plot.tag = element_text(size = 14, face = "bold"),
                   plot.tag.position = "topleft"))
  }
  img <- if (ext %in% c("jpg", "jpeg")) jpeg::readJPEG(path) else png::readPNG(path)
  ggplot() +
    annotation_raster(img, xmin = 0, xmax = 1, ymin = 0, ymax = 1) +
    xlim(0, 1) + ylim(0, 1) + theme_void() + labs(tag = tag) +
    theme(plot.margin = margin(5, 5, 5, 5),
          plot.tag = element_text(size = 14, face = "bold"),
          plot.tag.position = "topleft")
}

save_fig <- function(plot, name, width, height, dpi = 300) {
  ggsave(paste0(base_path, name, ".png"), plot, width = width, height = height, dpi = dpi)
  tryCatch(
    ggsave(paste0(base_path, name, ".pdf"), plot,
           width = width, height = height, dpi = dpi, device = cairo_pdf),
    error = function(e) {
      tryCatch(
        ggsave(paste0(base_path, name, ".pdf"), plot, width = width, height = height, dpi = dpi),
        error = function(e2) message("PDF failed: ", e2$message, " — PNG saved.")
      )
    }
  )
}

save_page <- function(plot_obj, filename, w, h) {
  path <- paste0(base_path, filename)
  tryCatch(
    ggsave(path, plot_obj, width = w, height = h, device = cairo_pdf),
    error = function(e) ggsave(path, plot_obj, width = w, height = h)
  )
  path
}

########################################################################
# PLOT FUNCTIONS
########################################################################

plot_paired_rawp <- function(data, id_col, trt_col, val_col, ylab, tag,
                             p_override = NULL) {
  df <- data.frame(id = data[[id_col]], trt = data[[trt_col]], val = data[[val_col]])
  df <- df[!is.na(df$val), ]
  df$trt <- factor(df$trt, levels = c("Placebo", "Ertugliflozin"))
  
  summ <- df %>%
    group_by(trt) %>%
    summarise(m = mean(val), se = sd(val) / sqrt(n()), .groups = "drop")
  
  pv <- if (!is.null(p_override)) p_override else {
    wide <- df %>% pivot_wider(names_from = trt, values_from = val) %>% arrange(id)
    t.test(wide$Placebo, wide$Ertugliflozin, paired = TRUE)$p.value
  }
  
  plab <- ifelse(pv < 0.001, "p<0.001", paste0("p=", sprintf("%.3f", pv)))
  
  ymax <- max(df$val, na.rm = TRUE)
  ymin <- min(df$val, na.rm = TRUE)
  yr   <- ymax - ymin
  
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

plot_three_group <- function(data, group_var, value_var,
                             group_colors, ylab, tag, group_order = NULL) {
  df <- data.frame(grp = data[[group_var]], val = data[[value_var]])
  df <- df[!is.na(df$val) & !is.na(df$grp), ]
  df$grp <- factor(df$grp, levels = if (!is.null(group_order)) group_order else levels(factor(df$grp)))
  lvls   <- levels(df$grp)
  
  kw_p   <- kruskal.test(val ~ grp, data = df)$p.value
  kw_lab <- ifelse(kw_p < 0.0001, "KW p<0.0001",
                   ifelse(kw_p < 0.001,  "KW p<0.001",
                          paste0("KW p=", sprintf("%.3f", kw_p))))
  
  pairs    <- combn(lvls, 2, simplify = FALSE)
  pair_pvs <- sapply(pairs, function(pr)
    wilcox.test(val ~ grp, data = df[df$grp %in% pr, ], exact = FALSE)$p.value)
  
  ymax    <- max(df$val, na.rm = TRUE)
  ymin    <- min(df$val, na.rm = TRUE)
  yr      <- ymax - ymin
  n_pairs <- length(pairs)
  step    <- yr * 0.13
  y_brack <- ymax + yr * 0.05 + (seq_len(n_pairs) - 1) * step
  
  p <- ggplot(df, aes(grp, val, fill = grp)) +
    geom_boxplot(alpha = 0.3, outlier.shape = NA, width = 0.5) +
    geom_jitter(aes(color = grp), width = 0.12, size = 2, alpha = 0.7) +
    annotate("text", x = 2,
             y = ymax + yr * 0.05 + n_pairs * step + yr * 0.04,
             label = kw_lab, size = 3.5, fontface = "italic") +
    scale_fill_manual(values  = group_colors) +
    scale_color_manual(values = group_colors) +
    coord_cartesian(ylim = c(ymin - yr * 0.05,
                             ymax + yr * 0.05 + n_pairs * step + yr * 0.14)) +
    labs(x = NULL, y = ylab, tag = tag) +
    theme_rockies
  
  for (k in seq_along(pairs)) {
    pr   <- pairs[[k]]
    pv   <- pair_pvs[k]
    plab <- ifelse(pv < 0.0001, "p<0.0001",
                   ifelse(pv < 0.001,  "p<0.001",
                          paste0("p=", sprintf("%.3f", pv))))
    x1 <- which(lvls == pr[1])
    x2 <- which(lvls == pr[2])
    ys <- y_brack[k]
    p  <- p +
      annotate("segment", x = x1, xend = x2, y = ys, yend = ys, linewidth = 0.4) +
      annotate("text", x = (x1 + x2) / 2, y = ys + yr * 0.03,
               label = plab, size = 3, fontface = "italic")
  }
  p
}

plot_two_group <- function(data, group_var, value_var,
                           group_colors, ylab, tag, test = "wilcox.test") {
  df <- data.frame(grp = data[[group_var]], val = data[[value_var]])
  df <- df[!is.na(df$val) & !is.na(df$grp), ]
  df$grp <- factor(df$grp, levels = names(group_colors))
  
  pv   <- if (test == "wilcox.test") wilcox.test(val ~ grp, data = df)$p.value
  else t.test(val ~ grp, data = df)$p.value
  plab <- ifelse(pv < 0.0001, "p<0.0001",
                 ifelse(pv < 0.001,  "p<0.001",
                        paste0("p=", sprintf("%.3f", pv))))
  
  ymax <- max(df$val, na.rm = TRUE)
  yr   <- diff(range(df$val, na.rm = TRUE))
  
  ggplot(df, aes(grp, val, fill = grp)) +
    geom_boxplot(alpha = 0.3, outlier.shape = NA, width = 0.5) +
    geom_jitter(aes(color = grp), width = 0.12, size = 2, alpha = 0.7) +
    annotate("segment", x = 1, xend = 2,
             y = ymax + yr * 0.05, yend = ymax + yr * 0.05, linewidth = 0.4) +
    annotate("text", x = 1.5, y = ymax + yr * 0.1,
             label = plab, size = 3.5, fontface = "italic") +
    scale_fill_manual(values  = group_colors) +
    scale_color_manual(values = group_colors) +
    labs(x = NULL, y = ylab, tag = tag) +
    theme_rockies
}

plot_corr_log <- function(data, xvar_raw, xvar_log, yvar, xlab, ylab, tag,
                          point_col = cols$corr_point, method = "spearman") {
  df <- data.frame(x_raw = data[[xvar_raw]], x_log = data[[xvar_log]], y = data[[yvar]])
  df <- df[!is.na(df$x_raw) & !is.na(df$x_log) & !is.na(df$y) & is.finite(df$x_log), ]
  ct  <- cor.test(df$x_raw, df$y, method = method)
  lab <- paste0(ifelse(method == "spearman", "rho", "r"), "=",
                sprintf("%.2f", ct$estimate), ", p=",
                ifelse(ct$p.value < 0.001, "<0.001", sprintf("%.3f", ct$p.value)))
  
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

plot_delta_corr <- function(data, id_col, trt_col, xvar, yvar,
                            xlab, ylab, tag, method = "spearman") {
  df <- data.frame(id = data[[id_col]], trt = data[[trt_col]],
                   x  = data[[xvar]],  y   = data[[yvar]])
  df <- df[!is.na(df$x) & !is.na(df$y), ]
  
  wx <- df %>% dplyr::select(id, trt, x) %>% pivot_wider(names_from = trt, values_from = x)
  wy <- df %>% dplyr::select(id, trt, y) %>% pivot_wider(names_from = trt, values_from = y)
  
  deltas <- data.frame(dx = wx$Ertugliflozin - wx$Placebo,
                       dy = wy$Ertugliflozin - wy$Placebo)
  deltas <- deltas[!is.na(deltas$dx) & !is.na(deltas$dy), ]
  
  ct  <- cor.test(deltas$dx, deltas$dy, method = method)
  lab <- paste0(ifelse(method == "spearman", "rho", "r"), "=",
                sprintf("%.2f", ct$estimate), ", p=",
                ifelse(ct$p.value < 0.001, "<0.001", sprintf("%.3f", ct$p.value)))
  
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

make_barplot <- function(filepath, celltype_label, panel_tag) {
  dt <- fread(filepath)
  setnames(dt, "logFC_groupType_2_Diabetes", "logFC")
  setnames(dt, "p_groupType_2_Diabetes",     "pvalue")
  dt <- dt %>%
    mutate(
      direction = factor(ifelse(logFC >= 0, "Up-regulated", "Down-regulated"),
                         levels = c("Up-regulated", "Down-regulated")),
      sig_label = ifelse(pvalue < 0.05, "*", ""),
      label_y   = ifelse(logFC >= 0,
                         logFC + max(abs(logFC), na.rm = TRUE) * 0.04,
                         logFC - max(abs(logFC), na.rm = TRUE) * 0.04)
    )
  ggplot(dt, aes(x = gene, y = logFC, fill = direction)) +
    geom_col(width = 0.7) +
    geom_text(aes(y = label_y, label = sig_label), size = 4, fontface = "bold", vjust = 0.5) +
    geom_hline(yintercept = 0, linewidth = 0.4) +
    scale_fill_manual(name   = "Regulation",
                      values = c("Up-regulated" = cols$up, "Down-regulated" = cols$down),
                      drop   = FALSE) +
    labs(x = NULL, y = "Log2 Fold Change", subtitle = celltype_label, tag = panel_tag) +
    theme_classic(base_size = 11) +
    theme(axis.text.x      = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7),
          axis.title.y     = element_text(size = 10, face = "bold"),
          plot.subtitle    = element_text(size = 10, hjust = 0.5, face = "bold"),
          plot.tag         = element_text(size = 14, face = "bold"),
          plot.tag.position = "topleft",
          legend.position  = "none",
          plot.margin      = margin(5, 10, 5, 10))
}

########################################################################
# DATA LOADING
########################################################################

sav_path <- "C:/Users/netio/Downloads/ROCKIES_Hoofddatabase__jan2026.sav"
if (!file.exists(sav_path)) stop("SAV file not found: ", sav_path)

rockies_wide <- read_sav(sav_path) %>% zap_labels()
names(rockies_wide) <- make.names(trimws(names(rockies_wide)))
rockies_crossover <- rockies_wide %>% filter(Group == 1)

rockies_long <- bind_rows(
  rockies_crossover %>% transmute(
    id = Participant, treatment = "Placebo",
    cortical_k2  = K2_cortex_mean_Sherbrook_PLB,    medullary_k2  = K2_medulla_mean_Sherbrook_PLB,
    cortical_k2f = K2_F_cortex_mean_Sherbrook_PLB,  medullary_k2f = K2_F_medulla_mean_Sherbrook_PLB,
    homa_ir = HOMAIR_Placebo, matsuda = Matsuda_Placebo, ogis = OGIS_placebo,
    gfr = GFR_Placebo_ml_min, erpf = ERPF_Placebo_ml_min, ff = FF_Placebo,
    tna_sodium = TNA_sodium_plb_mmolpermin,
    sbp = SBP_Placebo, dbp = DBP_Placebo, weight = Bodyweight_Placebo, bmi = BMI_Placebo,
    uacr = Ur_albuminuria_24hr_placebo, ur_sodium = Ur_sodium_24hr_placebo_mmol,
    ur_glucose = Ur_glucose_24hr_placebo_mmol
  ),
  rockies_crossover %>% transmute(
    id = Participant, treatment = "Ertugliflozin",
    cortical_k2  = K2_cortex_mean_Sherbrook_SGLT2,   medullary_k2  = K2_medulla_mean_Sherbrook_SGLT2i,
    cortical_k2f = K2_F_cortex_mean_Sherbrook_SGLT2, medullary_k2f = K2_F_medulla_mean_Sherbrook_SGLT2,
    homa_ir = HOMAIR_SGLT2i, matsuda = Matsuda_SGLT2i, ogis = OGIS_SGLT2i,
    gfr = GFR_SGLT2i_ml_min, erpf = ERPF_SGLT2i_ml_min, ff = FF_SGLT2,
    tna_sodium = TNA_sodium_SLGT2i_mmolpermin,
    sbp = SBP_verum, dbp = DBP_Verum, weight = Bodyweight_Verum, bmi = BMI_Verum,
    uacr = Ur_albuminuria_24hr_SGLT2i, ur_sodium = Ur_sodium_24hr_sglt2i_mmol,
    ur_glucose = Ur_glucose_24hr_sglt2i_mmol
  )
)
rockies_long$treatment <- factor(rockies_long$treatment, levels = c("Placebo", "Ertugliflozin"))

harmonized_data <- read.csv(
  "C:/Users/netio/OneDrive - UW/Laura Pyle's files - Biostatistics Core Shared Drive/Data Harmonization/Data Clean/harmonized_dataset.csv",
  na = ""
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
dat_clinical  <- dat_clinical[-str_which(dat_clinical$record_id, "-O")]

PET_avg <- function(data) {
  data.frame(
    avg_c_k2   = rowMeans(data[, c("lc_k2", "rc_k2")], na.rm = TRUE),
    avg_m_k2   = rowMeans(data[, c("lm_k2", "rm_k2")], na.rm = TRUE),
    avg_c_f    = rowMeans(data[, c("lc_f",  "rc_f")],  na.rm = TRUE),
    avg_m_f    = rowMeans(data[, c("lm_f",  "rm_f")],  na.rm = TRUE),
    avg_c_k2_f = rowMeans(data[, c("lc_k2", "rc_k2")], na.rm = TRUE) /
      rowMeans(data[, c("lc_f",  "rc_f")],  na.rm = TRUE),
    avg_m_k2_f = rowMeans(data[, c("lm_k2", "rm_k2")], na.rm = TRUE) /
      rowMeans(data[, c("lm_f",  "rm_f")],  na.rm = TRUE)
  )
}

pet_col_names <- c("avg_c_k2", "avg_m_k2", "avg_c_f", "avg_m_f", "avg_c_k2_f", "avg_m_k2_f")
dat_clean    <- dat[, !(names(dat) %in% pet_col_names), drop = FALSE]
dat_with_pet <- cbind(dat_clean, PET_avg(dat_clean))
dat_pet_slim <- dat_with_pet[!is.na(dat_with_pet$avg_c_k2) & is.finite(dat_with_pet$avg_c_k2),
                             c("record_id", "avg_c_k2", "avg_c_f", "avg_c_k2_f")]

########################################################################
# FIGURE 1
########################################################################

fig1a <- load_image_panel("C:/Users/netio/Downloads/ROCKIES study plan_v5.png", "A")
fig1b <- load_image_panel("C:/Users/netio/Downloads/Kidney C11 Acetate PET Diagram (3).png", "B")

# Panel C: PET k2 maps (coord_fixed preserved)
pet_map_img <- readPNG("C:/Users/netio/Downloads/Re_ ROCKIES_RH2_ Representative plots/k2-maps.png")
fig1c <- ggplot() +
  annotation_raster(pet_map_img, xmin = 0, xmax = 1, ymin = 0, ymax = 1) +
  xlim(0, 1) + ylim(0, 1) +
  coord_fixed(ratio = dim(pet_map_img)[1] / dim(pet_map_img)[2]) +
  theme_void() + labs(tag = "C") +
  theme(plot.tag = element_text(size = 14, face = "bold"),
        plot.tag.position = "topleft",
        plot.margin = margin(5, 5, 5, 5))

# Panel D: Baseline correlation heatmap
clin_labels <- c("HOMA-IR", "OGIS", "Body Weight")

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

pet_vars <- c("k2_cortex", "k2_medulla", "k2f_cortex", "k2f_medulla")
clin_vars <- c("homa_base", "ogis_base", "bw_base")

rho_mat <- matrix(NA_real_, 4, 3, dimnames = list(pet_labels, clin_labels))
p_mat   <- rho_mat

for (i in seq_along(pet_vars)) {
  for (j in seq_along(clin_vars)) {
    d <- rockies_baseline[, c(pet_vars[i], clin_vars[j])]
    d <- d[complete.cases(d), ]
    if (nrow(d) >= 5) {
      res           <- cor.test(d[[1]], d[[2]], method = "pearson", exact = FALSE)
      rho_mat[i, j] <- res$estimate
      p_mat[i, j]   <- res$p.value
    }
  }
}

sig_labels <- matrix("", 4, 3, dimnames = dimnames(rho_mat))
sig_labels[p_mat < 0.05  & !is.na(p_mat)] <- "*"
sig_labels[p_mat < 0.01  & !is.na(p_mat)] <- "**"
sig_labels[p_mat < 0.001 & !is.na(p_mat)] <- "***"

fig1d_long <- data.frame(
  pet  = factor(rep(pet_labels, times = length(clin_labels)), levels = rev(pet_labels)),
  clin = factor(rep(clin_labels, each  = length(pet_labels)), levels = clin_labels),
  rho  = as.vector(rho_mat),
  sig  = as.vector(sig_labels)
) %>% mutate(label = ifelse(is.na(rho), "", paste0(sprintf("%.2f", rho), sig)))

fig1d <- ggplot(fig1d_long, aes(x = clin, y = pet, fill = rho)) +
  geom_tile(color = "white", linewidth = 1) +
  geom_text(aes(label = label), size = 3.5, fontface = "bold") +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                       midpoint = 0, limits = c(-1, 1), name = "Pearson") +
  scale_y_discrete(labels = pet_md_labels) +
  labs(x = NULL, y = NULL, tag = "D") +
  theme_rockies +
  theme(legend.position = "right",
        axis.text.x = element_text(angle = 45, hjust = 1, size = 9, color = "black"),
        axis.text.y = element_markdown(size = 9))

# Panel E: Urine omics image
urine_img <- jpeg::readJPEG(
  "C:/Users/netio/Documents/UofW/Rockies/publication_figures/Urine_Metabolism_Figure.jpg")
fig1e <- ggplot() +
  annotation_raster(urine_img, xmin = 0, xmax = 1, ymin = 0, ymax = 1) +
  xlim(0, 1) + ylim(0, 1) +
  coord_fixed(ratio = dim(urine_img)[1] / dim(urine_img)[2]) +
  theme_void() + labs(tag = "E") +
  theme(plot.tag = element_text(size = 14, face = "bold"),
        plot.tag.position = "topleft",
        plot.margin = margin(5, 5, 5, 5))

# Panels F–I: Paired treatment plots
pval_raw <- function(plb_col, ertu_col)
  t.test(rockies_crossover[[plb_col]], rockies_crossover[[ertu_col]], paired = TRUE)$p.value

pv_f <- pval_raw("K2_cortex_mean_Sherbrook_PLB",    "K2_cortex_mean_Sherbrook_SGLT2")
pv_g <- pval_raw("K2_F_cortex_mean_Sherbrook_PLB",  "K2_F_cortex_mean_Sherbrook_SGLT2")
pv_h <- pval_raw("K2_medulla_mean_Sherbrook_PLB",   "K2_medulla_mean_Sherbrook_SGLT2i")
pv_i <- pval_raw("K2_F_medulla_mean_Sherbrook_PLB", "K2_F_medulla_mean_Sherbrook_SGLT2")

fig1f <- plot_paired_rawp(rockies_long, "id", "treatment", "cortical_k2",
                          expression(bold("Cortical k"[2]*" (min"^{-1}*")")), "G", p_override = pv_f)
fig1g <- plot_paired_rawp(rockies_long, "id", "treatment", "cortical_k2f",
                          expression(bold("Cortical k"[2]*"/F")), "H", p_override = pv_g)
fig1h <- plot_paired_rawp(rockies_long, "id", "treatment", "medullary_k2",
                          expression(bold("Medullary k"[2]*" (min"^{-1}*")")), "I", p_override = pv_h)
fig1i <- plot_paired_rawp(rockies_long, "id", "treatment", "medullary_k2f",
                          expression(bold("Medullary k"[2]*"/F")), "J", p_override = pv_i)

# Panel J: Delta-delta heatmap — Δ on both axes, k2 subscripts via element_markdown
dd_rho <- matrix(
  c(0.48, 0.77, 0.36, 0.20,
    0.50, 0.82, 0.31, 0.17,
    0.32, 0.61, 0.48, 0.45,
    0.31, 0.64, 0.42, 0.52),
  nrow = 4, ncol = 4, byrow = TRUE,
  dimnames = list(pet_labels,
                  c("Delta HbA1c", "Delta HOMA-IR", "Delta Glucose", "Delta Sodium Load"))
)
dd_cell <- matrix(
  c("0.48*", "0.77***", "0.36",  "0.20",
    "0.50*", "0.82***", "0.31",  "0.17",
    "0.32",  "0.61**",  "0.48*", "0.45*",
    "0.31",  "0.64**",  "0.42",  "0.52*"),
  nrow = 4, ncol = 4, byrow = TRUE, dimnames = dimnames(dd_rho)
)

delta_col_labels <- c(
  "Delta HbA1c"       = "\u0394 HbA1c",
  "Delta HOMA-IR"     = "\u0394 HOMA-IR",
  "Delta Glucose"     = "\u0394 Glucose",
  "Delta Sodium Load" = "\u0394 Sodium Load"
)

fig1j_long <- data.frame(
  pet   = factor(rep(pet_labels, times = ncol(dd_rho)), levels = rev(pet_labels)),
  delta = factor(rep(names(delta_col_labels), each = nrow(dd_rho)),
                 levels = names(delta_col_labels)),
  rho   = as.vector(dd_rho),
  cell  = as.vector(dd_cell)
)

fig1j <- ggplot(fig1j_long, aes(x = delta, y = pet, fill = rho)) +
  geom_tile(color = "white", linewidth = 1) +
  geom_text(aes(label = cell), size = 3.5, fontface = "bold") +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                       midpoint = 0, limits = c(-1, 1), name = "Pearson") +
  scale_x_discrete(labels = delta_col_labels) +
  scale_y_discrete(labels = pet_delta_md_labels) +
  labs(x = NULL, y = NULL, tag = "F") +
  theme_rockies +
  theme(legend.position = "right",
        axis.text.x = element_text(angle = 45, hjust = 1, size = 9, color = "black"),
        axis.text.y = element_markdown(size = 9))

# Assemble Figure 1
fig1_full <-
  (fig1a + fig1b + plot_layout(widths = c(1, 1))) /
  (fig1c + fig1d + plot_layout(widths = c(1.8, 1.2))) /
  (fig1e + fig1j + plot_layout(widths = c(1.4, 1.0))) /
  (fig1f + fig1g + fig1h + fig1i + plot_layout(widths = c(1, 1, 1, 1))) +
  plot_layout(heights = c(3.5, 3.5, 3.0, 2.5)) +
  plot_annotation(
    title = "Figure 1. SGLT2 Inhibition Reduces Kidney Oxidative Metabolism in the ROCKIES Trial",
    theme = theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5))
  )

save_fig(fig1_full, "Figure1_ROCKIES", width = 16, height = 22)
cat("Figure 1 saved!\n")

########################################################################
# FIGURE 2
########################################################################

dat_fig2_base <- dat_with_pet %>%
  filter(!is.na(avg_c_k2) & is.finite(avg_c_k2)) %>%
  filter(group %in% c("Lean Control", "Obese Control", "Type 2 Diabetes")) %>%
  mutate(Cohort = case_when(
    group == "Lean Control"    ~ "Healthy Control",
    group == "Obese Control"   ~ "Obese Control",
    group == "Type 2 Diabetes" ~ "T2D"
  )) %>%
  filter(record_id != 'CRC-55') %>% 
  filter(Cohort == "Healthy Control" | is.na(epic_sglti2_1) | epic_sglti2_1 != 'Yes')

cohort_ns <- dat_fig2_base %>% filter(!is.na(avg_c_k2)) %>% count(Cohort)
get_n     <- function(grp) cohort_ns$n[cohort_ns$Cohort == grp]

lbl_hc  <- paste0("Healthy Control\n(n=", get_n("Healthy Control"), ")")
lbl_oc  <- paste0("Obese Control\n(n=",   get_n("Obese Control"),   ")")
lbl_t2d <- paste0("T2D\n(n=",             get_n("T2D"),             ")")

dat_fig2 <- dat_fig2_base %>%
  mutate(
    Cohort_label = case_when(
      Cohort == "Healthy Control" ~ lbl_hc,
      Cohort == "Obese Control"   ~ lbl_oc,
      Cohort == "T2D"             ~ lbl_t2d
    ),
    Cohort_label = factor(Cohort_label, levels = c(lbl_hc, lbl_oc, lbl_t2d))
  )

grp_order   <- c(lbl_hc, lbl_oc, lbl_t2d)
fig2_colors <- setNames(c(cols$control, cols$obese_ctrl, cols$t2d), grp_order)

fig2a <- plot_three_group(dat_fig2, "Cohort_label", "avg_c_k2", fig2_colors,
                          ylab = expression(bold("Cortical k"[2]*" (min"^{-1}*")")), tag = "A", group_order = grp_order)
fig2b <- plot_three_group(dat_fig2, "Cohort_label", "avg_c_f", fig2_colors,
                          ylab = expression(bold("Cortical F (ml/g/min)")), tag = "B", group_order = grp_order)
fig2c <- plot_three_group(dat_fig2, "Cohort_label", "avg_c_k2_f", fig2_colors,
                          ylab = expression(bold("Cortical k"[2]*"/F")), tag = "C", group_order = grp_order)

fig2 <- (fig2a | fig2b | fig2c) +
  plot_annotation(
    title    = "Figure 2. Elevated Kidney Oxidative Metabolism in Type 2 Diabetes",
    subtitle = paste0("Healthy Control (n=", get_n("Healthy Control"),
                      "), Obese Control (n=", get_n("Obese Control"),
                      "), T2D (n=",           get_n("T2D"), ")"),
    theme = theme(
      plot.title    = element_text(size = 12, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 10, hjust = 0.5, color = "gray30")
    )
  )

save_fig(fig2, "Figure2_T2D_vs_Controls", width = 13, height = 6)
cat("Figure 2 saved!\n")

########################################################################
# FIGURE 3
########################################################################

dat_fig3 <- dat_with_pet %>%
  filter(!is.na(avg_c_k2) & is.finite(avg_c_k2)) %>%
  filter(group %in% c("Lean Control", "Type 2 Diabetes")) %>%
  filter(!is.na(acr_u)) %>%
  mutate(log10_uacr = log10(acr_u)) %>% 
  filter(group == "Lean Control" | is.na(epic_sglti2_1) | epic_sglti2_1 != 'Yes')

vars_to_test <- list(
  list(var = "avg_c_k2",   label = "Cortical k2"),
  list(var = "avg_c_f",    label = "Cortical F"),
  list(var = "avg_c_k2_f", label = "Cortical k2/F")
)

cor_results <- purrr::map_dfr(vars_to_test, function(v) {
  df <- dat_fig3 %>% filter(!is.na(.data[[v$var]]))
  ct <- cor.test(df$acr_u, df[[v$var]], method = "spearman")
  data.frame(
    var_label = v$label, rho = ct$estimate, pval = ct$p.value,
    label = sprintf("%.2f", ct$estimate),
    sig = ifelse(ct$p.value < 0.001, "***", ifelse(ct$p.value < 0.01, "**",
                                                   ifelse(ct$p.value < 0.05, "*", "")))
  )
})
cor_results$var_label <- factor(cor_results$var_label,
                                levels = c("Cortical k2/F", "Cortical F", "Cortical k2"))

fig3a_md_labels <- c(
  "Cortical k2"   = "Cortical k<sub>2</sub>",
  "Cortical F"    = "Cortical F",
  "Cortical k2/F" = "Cortical k<sub>2</sub>/F"
)

fig3a <- ggplot(cor_results, aes(x = "UACR", y = var_label, fill = rho)) +
  geom_tile(color = "white", linewidth = 1) +
  geom_text(aes(label = paste0(label, sig)), size = 4, fontface = "bold") +
  scale_fill_gradient2(low = cols$control, mid = "white", high = cols$t2d,
                       midpoint = 0, limits = c(-1, 1), name = "Spearman\nrho") +
  scale_y_discrete(labels = fig3a_md_labels) +
  labs(x = NULL, y = NULL, tag = "A") +
  theme_rockies +
  theme(legend.position = "right",
        axis.text.x = element_text(face = "bold"),
        axis.text.y = element_markdown(size = 9))

fig3b <- plot_corr_log(dat_fig3, "acr_u", "log10_uacr", "avg_c_k2",
                       expression(bold("log"[10]*"(UACR) (mg/g)")), expression(bold("Cortical k"[2]*" (min"^{-1}*")")), "B")
fig3c <- plot_corr_log(dat_fig3, "acr_u", "log10_uacr", "avg_c_f",
                       expression(bold("log"[10]*"(UACR) (mg/g)")), expression(bold("Cortical F (ml/g/min)")), "C")
fig3d <- plot_corr_log(dat_fig3, "acr_u", "log10_uacr", "avg_c_k2_f",
                       expression(bold("log"[10]*"(UACR) (mg/g)")), expression(bold("Cortical k"[2]*"/F")), "D")

fig3 <- (fig3a | fig3b) / (fig3c | fig3d) +
  plot_annotation(
    title    = "Figure 3. Kidney Oxidative Metabolism Correlates with Albuminuria",
    theme = theme(
      plot.title    = element_text(size = 14, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 11, hjust = 0.5, color = "gray30")
    )
  )

save_fig(fig3, "Figure3_UACR_Correlation", width = 12, height = 8)
cat("Figure 3 saved!\n")

########################################################################
# FIGURE 4
########################################################################

dat_fig4_clin <- dat_clinical %>%
  as.data.frame() %>%
  filter(`GBM thickness` != "" & !is.na(`GBM thickness`)) %>%
  mutate(GBM_Status = ifelse(`GBM thickness` == "yes", "GBM\nThickening", "No GBM\nThickening")) %>%
  dplyr::select(record_id, GBM_Status)
dat_fig4_clin <- dat_fig4_clin[, !(names(dat_fig4_clin) %in% pet_col_names), drop = FALSE]
dat_fig4 <- merge(dat_fig4_clin, dat_pet_slim, by = "record_id")

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
# FIGURE 5
########################################################################

dat_fig5_clin <- dat_clinical %>%
  as.data.frame() %>%
  filter(arteriosclerosis != "" & !is.na(arteriosclerosis)) %>%
  mutate(Arterio_Status = ifelse(arteriosclerosis == "yes", "Arteriosclerosis", "No\nArteriosclerosis")) %>%
  dplyr::select(record_id, Arterio_Status)
dat_fig5_clin <- dat_fig5_clin[, !(names(dat_fig5_clin) %in% pet_col_names), drop = FALSE]
dat_fig5 <- merge(dat_fig5_clin, dat_pet_slim, by = "record_id")

fig5_colors <- c("No\nArteriosclerosis" = cols$arterio_no, "Arteriosclerosis" = cols$arterio_yes)

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
# FIGURES 6 & 7
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

legend_dummy <- ggplot(
  data.frame(
    x         = 1:2,
    direction = factor(c("Up-regulated", "Down-regulated"),
                       levels = c("Up-regulated", "Down-regulated"))
  ),
  aes(x = x, y = 1, fill = direction)) +
  geom_tile() +
  scale_fill_manual(
    name   = "Regulation",
    values = c("Up-regulated" = cols$up, "Down-regulated" = cols$down)
  ) +
  theme_void() +
  theme(
    legend.position = "bottom",
    legend.title    = element_text(size = 10, face = "bold"),
    legend.text     = element_text(size = 10),
    legend.key.size = unit(0.6, "cm")
  )
shared_legend <- get_legend(legend_dummy)

fig6_panels <- make_barplot(tca_files[["PT"]], "PT", "A") /
  (make_barplot(tca_files[["PT-S1/S2"]], "PT-S1/S2", "B") |
     make_barplot(tca_files[["PT-S3"]],   "PT-S3",     "C") |
     make_barplot(tca_files[["aPT"]],     "aPT",       "D")) +
  plot_layout(heights = c(1.3, 1)) +
  plot_annotation(
    title = "Figure 6. TCA Cycle Gene Expression in Proximal Tubule Cells",
    theme = theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5))
  )
fig6 <- plot_grid(as.grob(fig6_panels), shared_legend, ncol = 1, rel_heights = c(1, 0.05))

fig7_panels <- make_barplot(oxphos_files[["PT"]], "PT", "A") /
  (make_barplot(oxphos_files[["PT-S1/S2"]], "PT-S1/S2", "B") |
     make_barplot(oxphos_files[["PT-S3"]],   "PT-S3",     "C") |
     make_barplot(oxphos_files[["aPT"]],     "aPT",       "D")) +
  plot_layout(heights = c(1.3, 1)) +
  plot_annotation(
    title = "Figure 7. OxPhos Gene Expression in Proximal Tubule Cells",
    theme = theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5))
  )
fig7 <- plot_grid(as.grob(fig7_panels), shared_legend, ncol = 1, rel_heights = c(1, 0.05))

save_fig(fig6, "Figure6_TCA_Cycle", width = 14, height = 10.5)
save_fig(fig7, "Figure7_OxPhos",    width = 14, height = 10.5)
cat("Figures 6 & 7 saved!\n")

########################################################################
# COMBINED PDF
########################################################################

pages <- c(
  save_page(fig1_full, "_tmp_fig1.pdf", w = 16,  h = 22),
  save_page(fig2,      "_tmp_fig2.pdf", w = 13,  h = 6),
  save_page(fig3,      "_tmp_fig3.pdf", w = 12,  h = 9),
  save_page(fig4,      "_tmp_fig4.pdf", w = 13,  h = 5.5),
  save_page(fig5,      "_tmp_fig5.pdf", w = 13,  h = 5.5)
)

for (fig_obj in list(fig6, fig7)) {
  tmp <- paste0(base_path, "_tmp_fig", length(pages) + 1, ".pdf")
  pdf(tmp, width = 14, height = 10.5)
  print(fig_obj)
  dev.off()
  pages <- c(pages, tmp)
}

combined_path <- paste0(base_path, "ROCKIES_All_Figures_Combined_noSGLT2i.pdf")
qpdf::pdf_combine(pages, output = combined_path)
file.remove(pages)

cat("\nCombined PDF saved:", combined_path, "\n")

