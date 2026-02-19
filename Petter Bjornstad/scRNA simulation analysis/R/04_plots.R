################################################################################
# 04_plots.R
#
# PURPOSE:
#   Generate publication-quality plots comparing NEBULA vs DESeq2 vs edgeR
#   across all simulation parameters:
#     1. Power curves vs n_subjects_per_arm, effect_size, prop_de
#     2. FDR control plots (nominal vs observed)
#     3. Type I error rate under null (prop_de == 0)
#     4. Compute time comparison
#     5. AUC-ROC heatmap across parameter combinations
#     6. Power vs individual-level variance (indiv_var_label)
#     7. Power vs within-subject correlation (corr_cells)
#
# INPUT (from S3):
#   bucket: scrna
#   prefix: Projects/Paired scRNA simulation analysis/results/benchmark/
#     benchmark_avg.rds
#
# OUTPUT (to S3):
#   bucket: scrna
#   prefix: Projects/Paired scRNA simulation analysis/results/plots/
#     power_curves.pdf
#     fdr_control.pdf
#     t1e_null.pdf
#     compute_time.pdf
#     auc_heatmap.pdf
#     indiv_var_effect.pdf
#     corr_cells_effect.pdf
#
# USAGE:
#   Rscript 04_plots.R
################################################################################

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(optparse)
  library(patchwork)
  library(ggridges)
  library(viridis)
  library(scales)
  library(aws.s3)
})

# ── S3 / Multi-user setup ─────────────────────────────────────────────────────
setup_s3 <- function() {
  user <- Sys.info()[["user"]]

  if (user == "choiyej") {
    keys_path <- "/Users/choiyej/Library/CloudStorage/OneDrive-UW/YC_RK Lab/KPMP/s3/keys.json"
  } else if (user %in% c("rameshsh", "yejichoi", "pylell")) {
    keys_path <- "/gscratch/scrubbed/yejichoi/keys.json"
  } else {
    stop("Unknown user '", user, "'. Add credentials path to setup_s3().")
  }

  keys <- jsonlite::fromJSON(keys_path)
  Sys.setenv(
    AWS_ACCESS_KEY_ID     = keys$access_key,
    AWS_SECRET_ACCESS_KEY = keys$secret_key,
    AWS_S3_ENDPOINT       = "s3.kopah.uw.edu"
  )
  message(sprintf("S3 configured for user '%s'", user))
}

setup_s3()

# S3 paths
S3_BUCKET    <- "scrna"
S3_BASE      <- "Projects/Paired scRNA simulation analysis/results/"
S3_BENCH_PFX <- paste0(S3_BASE, "benchmark/")
S3_PLOTS_PFX <- paste0(S3_BASE, "plots/")

# ── S3 helper functions ───────────────────────────────────────────────────────
s3write_using_region <- function(FUN, ..., object, bucket,
                                  region = NULL, opts = NULL, filename = NULL) {
  ext  <- if (!is.null(filename)) tools::file_ext(filename) else tools::file_ext(object)
  tmp  <- tempfile(fileext = if (nchar(ext) > 0) paste0(".", ext) else "")
  on.exit(unlink(tmp))
  FUN(..., file = tmp)
  args <- list(file = tmp, object = object, bucket = bucket)
  if (!is.null(region)) args$region <- region
  if (!is.null(opts))   args        <- c(args, opts)
  do.call(aws.s3::put_object, args)
}

# Wrapper: ggsave uses 'filename' arg instead of 'file'
ggsave_to_s3 <- function(plot, s3_object, width = 10, height = 8,
                          bucket = S3_BUCKET, region = "") {
  tmp <- tempfile(fileext = ".pdf")
  on.exit(unlink(tmp))
  ggsave(filename = tmp, plot = plot, width = width, height = height)
  aws.s3::put_object(file = tmp, object = s3_object, bucket = bucket,
                     region = region)
  message(sprintf("  Saved: s3://%s/%s", bucket, s3_object))
}

# ── CLI args ──────────────────────────────────────────────────────────────────
option_list <- list(
  make_option("--fdr_nominal", type = "double", default = 0.05)
)
opt <- parse_args(OptionParser(option_list = option_list))

# ── Load benchmark data from S3 ───────────────────────────────────────────────
message("── [06] Loading benchmark_avg from S3 ──")
avg <- s3readRDS(
  object = paste0(S3_BENCH_PFX, "benchmark_avg.rds"),
  bucket = S3_BUCKET,
  region = ""
)

METHOD_COLORS <- c(nebula = "#E64B35", deseq2 = "#4DBBD5", edger = "#00A087")
METHOD_LABELS <- c(nebula = "NEBULA-LN", deseq2 = "DESeq2 PB", edger = "edgeR PB")

base_theme <- theme_bw(base_size = 12) +
  theme(strip.background = element_rect(fill = "grey90"),
        legend.position  = "bottom")

# ── 1. Power curves ───────────────────────────────────────────────────────────
message("  Plotting power curves...")

pw_dat <- avg %>%
  filter(prop_de > 0) %>%
  mutate(effect_size_label = factor(effect_size_label, levels = c("med", "high")))

p_power <- ggplot(pw_dat,
                  aes(x = n_subjects_per_arm, y = mean_power,
                      colour = method, fill = method,
                      group = method)) +
  geom_ribbon(aes(ymin = mean_power - sd_power,
                  ymax = mean_power + sd_power), alpha = 0.15, colour = NA) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  facet_grid(prop_de ~ effect_size_label,
             labeller = label_both) +
  scale_colour_manual(values = METHOD_COLORS, labels = METHOD_LABELS) +
  scale_fill_manual(values   = METHOD_COLORS, labels = METHOD_LABELS) +
  scale_y_continuous(limits = c(0, 1), labels = percent) +
  labs(title  = "Power vs Sample Size",
       x      = "Subjects per arm (each with PRE + POST)",
       y      = "Power (mean ± SD across replicates)",
       colour = NULL, fill = NULL) +
  base_theme

ggsave_to_s3(p_power, paste0(S3_PLOTS_PFX, "power_curves.pdf"),
             width = 10, height = 8)

# ── 2. FDR control ────────────────────────────────────────────────────────────
message("  Plotting FDR control...")

fdr_dat <- avg %>% filter(prop_de > 0)

p_fdr <- ggplot(fdr_dat,
                aes(x = factor(n_subjects_per_arm), y = mean_fdr,
                    colour = method, group = method)) +
  geom_hline(yintercept = opt$fdr_nominal, linetype = "dashed",
             colour = "red", linewidth = 0.7) +
  geom_point(size = 2, position = position_dodge(0.4)) +
  geom_errorbar(aes(ymin = mean_fdr - sd_fdr, ymax = mean_fdr + sd_fdr),
                width = 0.3, position = position_dodge(0.4)) +
  facet_grid(prop_de ~ effect_size_label, labeller = label_both) +
  scale_colour_manual(values = METHOD_COLORS, labels = METHOD_LABELS) +
  scale_y_continuous(labels = percent) +
  labs(title  = "Observed FDR vs Nominal 5%",
       x      = "Subjects per arm",
       y      = "Observed FDR",
       colour = NULL) +
  base_theme

ggsave_to_s3(p_fdr, paste0(S3_PLOTS_PFX, "fdr_control.pdf"),
             width = 10, height = 8)

# ── 3. Type I error under null ────────────────────────────────────────────────
message("  Plotting Type I error...")

t1e_dat <- avg %>% filter(prop_de == 0)

p_t1e <- ggplot(t1e_dat,
                aes(x = factor(n_subjects_per_arm), y = mean_t1e,
                    colour = method, group = method)) +
  geom_hline(yintercept = opt$fdr_nominal, linetype = "dashed",
             colour = "red", linewidth = 0.7) +
  geom_point(size = 2, position = position_dodge(0.4)) +
  geom_errorbar(aes(ymin = mean_t1e - sd_t1e, ymax = mean_t1e + sd_t1e),
                width = 0.3, position = position_dodge(0.4)) +
  facet_grid(indiv_var_label ~ corr_cells,
             labeller = label_both) +
  scale_colour_manual(values = METHOD_COLORS, labels = METHOD_LABELS) +
  scale_y_continuous(labels = percent) +
  labs(title  = "Type I Error Rate (null: prop_de = 0)",
       x      = "Subjects per arm",
       y      = "Type I error rate",
       colour = NULL) +
  base_theme

ggsave_to_s3(p_t1e, paste0(S3_PLOTS_PFX, "t1e_null.pdf"),
             width = 10, height = 8)

# ── 4. Compute time ───────────────────────────────────────────────────────────
message("  Plotting compute time...")

time_dat <- avg %>%
  select(method, n_subjects_per_arm, cells_label, mean_compute_s) %>%
  mutate(method = factor(method))

p_time <- ggplot(time_dat,
                 aes(x = factor(n_subjects_per_arm), y = mean_compute_s,
                     colour = method, group = method)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  facet_wrap(~cells_label, labeller = label_both) +
  scale_colour_manual(values = METHOD_COLORS, labels = METHOD_LABELS) +
  scale_y_log10(labels = label_number(suffix = " s")) +
  labs(title  = "Compute Time per Scenario",
       x      = "Subjects per arm",
       y      = "Wall time (seconds, log scale)",
       colour = NULL) +
  base_theme

ggsave_to_s3(p_time, paste0(S3_PLOTS_PFX, "compute_time.pdf"),
             width = 8, height = 5)

# ── 5. Individual-level variance effect on power ──────────────────────────────
message("  Plotting individual variance effect on power...")

iv_dat <- avg %>%
  filter(prop_de > 0, effect_size_label == "med") %>%
  mutate(indiv_var_label = factor(indiv_var_label, levels = c("no","med","high")))

p_iv <- ggplot(iv_dat,
               aes(x = indiv_var_label, y = mean_power,
                   colour = method, group = method)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.5) +
  geom_errorbar(aes(ymin = mean_power - sd_power,
                    ymax = mean_power + sd_power),
                width = 0.2) +
  facet_grid(n_subjects_per_arm ~ prop_de, labeller = label_both) +
  scale_colour_manual(values = METHOD_COLORS, labels = METHOD_LABELS) +
  scale_y_continuous(limits = c(0,1), labels = percent) +
  labs(title  = "Power vs Individual-Level Variance",
       x      = "Individual variance level",
       y      = "Power",
       colour = NULL) +
  base_theme

ggsave_to_s3(p_iv, paste0(S3_PLOTS_PFX, "indiv_var_effect.pdf"),
             width = 10, height = 8)

# ── 6. Cell-cell correlation effect on power ─────────────────────────────────
message("  Plotting correlation effect on power...")

corr_dat <- avg %>%
  filter(prop_de > 0, effect_size_label == "med")

p_corr <- ggplot(corr_dat,
                 aes(x = factor(corr_cells), y = mean_power,
                     colour = method, group = method)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.5) +
  facet_grid(n_subjects_per_arm ~ prop_de, labeller = label_both) +
  scale_colour_manual(values = METHOD_COLORS, labels = METHOD_LABELS) +
  scale_y_continuous(limits = c(0,1), labels = percent) +
  labs(title  = "Power vs Within-Subject Cell-Cell Correlation",
       x      = "Correlation among cells within subject",
       y      = "Power",
       colour = NULL) +
  base_theme

ggsave_to_s3(p_corr, paste0(S3_PLOTS_PFX, "corr_cells_effect.pdf"),
             width = 10, height = 8)

# ── 7. AUC-ROC heatmap ────────────────────────────────────────────────────────
message("  Plotting AUC-ROC heatmap...")

auc_dat <- avg %>%
  filter(prop_de > 0) %>%
  mutate(param_combo = paste0(effect_size_label, "\n",
                              indiv_var_label, " var\n",
                              "corr=", corr_cells))

p_auc <- ggplot(auc_dat,
                aes(x = param_combo, y = method, fill = mean_auc_roc)) +
  geom_tile(colour = "white") +
  scale_fill_viridis(option = "C", limits = c(0.5, 1),
                     name = "AUC-ROC") +
  facet_wrap(~ n_subjects_per_arm, labeller = label_both) +
  labs(title = "AUC-ROC Heatmap",
       x = "Parameter combination", y = "Method") +
  base_theme +
  theme(axis.text.x = element_text(size = 7, angle = 30, hjust = 1))

ggsave_to_s3(p_auc, paste0(S3_PLOTS_PFX, "auc_heatmap.pdf"),
             width = 14, height = 6)

message(sprintf("── [06] All plots saved to s3://%s/%s ──", S3_BUCKET, S3_PLOTS_PFX))
