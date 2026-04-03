########################################################################
# Figure 2 Standalone Script
# T2D vs Healthy Controls vs Obese Controls — k2, F, k2/F
########################################################################

library(tidyverse)
library(patchwork)
library(data.table)

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
    text             = element_text(family = font_family),
    axis.title       = element_text(size = 10, face = "bold"),
    axis.text        = element_text(size = 9, color = "black"),
    plot.title       = element_text(size = 11, face = "bold", hjust = 0.5),
    plot.tag         = element_text(size = 14, face = "bold"),
    strip.text       = element_text(size = 10, face = "bold"),
    strip.background = element_blank(),
    legend.position  = "none",
    plot.margin      = margin(5, 10, 5, 10)
  )

cols <- list(
  placebo      = "#4A90D9", ertu = "#E74C3C", paired_line = "gray60",
  t2d          = "#E74C3C", control = "#4A90D9", obese_ctrl = "#F39C12",
  corr_point   = "#2C3E50", delta = "#8E44AD",
  gbm_yes      = "#E74C3C", gbm_no = "#4A90D9",
  arterio_yes  = "#E74C3C", arterio_no = "#4A90D9",
  up           = "#E74C3C", down = "#4A90D9", ns = "gray70"
)

base_path <- 'C:/Users/netio/Documents/UofW/Rockies/publication_figures/'

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

harmonized_data <- read.csv(
  "C:/Users/netio/OneDrive - UW/Laura Pyle's files - Biostatistics Core Shared Drive/Data Harmonization/Data Clean/harmonized_dataset.csv",
  na = ''
)

dat <- harmonized_data %>%
  dplyr::select(-dob) %>%
  arrange(date_of_screen) %>%
  dplyr::summarise(
    across(where(negate(is.numeric)), ~ ifelse(all(is.na(.x)), NA_character_, last(na.omit(.x)))),
    across(where(is.numeric),         ~ ifelse(all(is.na(.x)), NA_real_,      mean(na.omit(.x), na.rm = TRUE))),
    .by = c(record_id, visit)
  )

########################################################################
# PET AVERAGES
########################################################################

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

pet_col_names <- c('avg_c_k2','avg_m_k2','avg_c_f','avg_m_f','avg_c_k2_f','avg_m_k2_f')
dat_clean    <- dat[, !(names(dat) %in% pet_col_names), drop = FALSE]
dat_with_pet <- cbind(dat_clean, PET_avg(dat_clean))

cat("dat_with_pet:", nrow(dat_with_pet), "rows\n")

########################################################################
# PLOT FUNCTION
########################################################################

plot_two_group <- function(data, group_var, value_var,
                           group_colors, ylab, tag) {
  df <- data.frame(grp = data[[group_var]], val = data[[value_var]])
  df <- df[!is.na(df$val) & !is.na(df$grp), ]
  df$grp <- factor(df$grp, levels = names(group_colors))
  
  n_groups <- nlevels(df$grp)
  ymax <- max(df$val, na.rm = TRUE)
  yr   <- diff(range(df$val, na.rm = TRUE))
  
  p <- ggplot(df, aes(grp, val, fill = grp)) +
    geom_boxplot(alpha = 0.3, outlier.shape = NA, width = 0.5) +
    geom_jitter(aes(color = grp), width = 0.12, size = 2, alpha = 0.7) +
    scale_fill_manual(values  = group_colors) +
    scale_color_manual(values = group_colors) +
    labs(x = NULL, y = ylab, tag = tag) +
    theme_rockies
  
  if (n_groups == 2) {
    # Simple Wilcoxon for two groups
    pv   <- wilcox.test(val ~ grp, data = df)$p.value
    plab <- ifelse(pv < 0.0001, "p<0.0001",
                   ifelse(pv < 0.001,  "p<0.001",
                          paste0("p=", formatC(pv, format = "f", digits = 3))))
    p <- p +
      annotate("segment", x = 1, xend = 2,
               y = ymax + yr * 0.05, yend = ymax + yr * 0.05, linewidth = 0.4) +
      annotate("text", x = 1.5, y = ymax + yr * 0.1,
               label = plab, size = 3, fontface = "italic") +
      coord_cartesian(ylim = c(NA, ymax + yr * 0.18))
    
  } else {
    # Kruskal-Wallis overall + pairwise Wilcoxon with BH correction
    kw_p  <- kruskal.test(val ~ grp, data = df)$p.value
    kw_lab <- ifelse(kw_p < 0.0001, "KW p<0.0001",
                     ifelse(kw_p < 0.001,  "KW p<0.001",
                            paste0("KW p=", formatC(kw_p, format = "f", digits = 3))))
    
    pw <- pairwise.wilcox.test(df$val, df$grp, p.adjust.method = "BH")$p.value
    
    fmt_p <- function(pv) {
      if (is.na(pv))        return(NA_character_)
      if (pv < 0.0001)      return("p<0.0001")
      if (pv < 0.001)       return("p<0.001")
      paste0("p=", formatC(pv, format = "f", digits = 3))
    }
    
    # bracket heights: comparisons 1-2, 1-3, 2-3
    b1 <- ymax + yr * 0.08   # groups 1 vs 2
    b2 <- ymax + yr * 0.20   # groups 1 vs 3
    b3 <- ymax + yr * 0.34   # groups 2 vs 3
    
    lvls <- levels(df$grp)
    
    get_pw_p <- function(g1, g2) {
      # pw matrix: rows = later groups, cols = earlier groups
      r <- match(g2, rownames(pw)); c <- match(g1, colnames(pw))
      if (!is.na(r) && !is.na(c)) return(pw[r, c])
      r <- match(g1, rownames(pw)); c <- match(g2, colnames(pw))
      if (!is.na(r) && !is.na(c)) return(pw[r, c])
      return(NA_real_)
    }
    
    p12 <- fmt_p(get_pw_p(lvls[1], lvls[2]))
    p13 <- fmt_p(get_pw_p(lvls[1], lvls[3]))
    p23 <- fmt_p(get_pw_p(lvls[2], lvls[3]))
    
    p <- p +
      # bracket 1–2
      annotate("segment", x = 1, xend = 2, y = b1, yend = b1, linewidth = 0.4) +
      annotate("text", x = 1.5, y = b1 + yr * 0.03, label = p12, size = 2.8, fontface = "italic") +
      # bracket 1–3
      annotate("segment", x = 1, xend = 3, y = b2, yend = b2, linewidth = 0.4) +
      annotate("text", x = 2,   y = b2 + yr * 0.03, label = p13, size = 2.8, fontface = "italic") +
      # bracket 2–3
      annotate("segment", x = 2, xend = 3, y = b3, yend = b3, linewidth = 0.4) +
      annotate("text", x = 2.5, y = b3 + yr * 0.03, label = p23, size = 2.8, fontface = "italic") +
      # overall KW p in top-left
      annotate("text", x = 0.6, y = b3 + yr * 0.03,
               label = kw_lab, size = 2.8, fontface = "italic", hjust = 0) +
      coord_cartesian(ylim = c(NA, b3 + yr * 0.10))
  }
  
  p
}

########################################################################
# FIGURE 2 DATA
########################################################################

dat_fig2 <- dat_with_pet %>%
  filter(!is.na(avg_c_k2) & is.finite(avg_c_k2)) %>%
  filter(group %in% c('Lean Control', 'Type 2 Diabetes', 'Obese Control')) %>%
  mutate(Cohort = case_when(
    group == "Type 2 Diabetes" ~ "T2D\n(RENAL-HEIRitage)",
    group == "Lean Control"    ~ "Healthy Control\n(CROCODILE)",
    group == "Obese Control"   ~ "Obese Control\n(CROCODILE)"
  ))

cat("Figure 2 data:", nrow(dat_fig2), "rows\n")

# Print IDs for each group
for (grp in unique(dat_fig2$Cohort)) {
  ids <- dat_fig2 %>% filter(Cohort == grp) %>% pull(record_id)
  cat("\n", grp, "IDs (n =", length(ids), "):\n")
  cat(paste(ids, collapse = ", "), "\n")
}

########################################################################
# FIGURE 2 PANELS
########################################################################

fig2_colors <- c(
  "Healthy Control\n(CROCODILE)" = cols$control,
  "Obese Control\n(CROCODILE)"   = cols$obese_ctrl,
  "T2D\n(RENAL-HEIRitage)"       = cols$t2d
)

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

########################################################################
# SAVE — combined + individual panels
########################################################################

save_fig(fig2,  "Figure2_T2D_vs_Controls", width = 12, height = 5)
save_fig(fig2b, "Figure2b_k2",             width = 4,  height = 5)
save_fig(fig2c, "Figure2c_F",              width = 4,  height = 5)
save_fig(fig2d, "Figure2d_k2_F",           width = 4,  height = 5)

cat("Figure 2 saved!\n")



