# ---
# limma baseline analysis (unadjusted and adjusted) for CKM stages
# Author: Darwin Del Castillo
# Date: `r lubridate::today()`
# ---

#### Plasma Baseline Analysis using limma ####

#######################
# Unadjusted Analysis #
#######################

# Identify plasma metabolite columns
seq_plasma <- is_seq_plasma(names(cleaned_log_baseline_plasma_data))

# create design matrix for unadjusted analysis (ANOVA-style)
design_unadj <- model.matrix(~ ckm_syn_base_factor,
                             data = cleaned_log_baseline_plasma_data)

# prepare response matrix
ymat <- t(cleaned_log_baseline_plasma_data[, seq_plasma])

# fit unadjusted models
fit <- lmFit(ymat, design_unadj)
fit <- eBayes(fit)

# Test for overall group differences (ANOVA F-test)
# This tests if any of the groups differ from each other
results_ckm_anova <- topTable(fit, coef = 2:3, number = nrow(ymat))
results_ckm_anova <- results_ckm_anova[order(results_ckm_anova$P.Value), ]

# Individual contrasts (pairwise comparisons)
# Stage 2+ vs Stage 2
results_ckm_stage2plus <- topTable(fit, coef = 2, number = nrow(ymat))
# Stage 3 vs Stage 2
results_ckm_stage3 <- topTable(fit, coef = 3, number = nrow(ymat))

# Save unadjusted results
write.csv(results_ckm_anova, "results/baseline_unadj_anova.csv",
          row.names = TRUE)
write.csv(results_ckm_stage2plus, "results/baseline_unadj_stage2plus.csv",
          row.names = TRUE)
write.csv(results_ckm_stage3, "results/baseline_unadj_stage3.csv",
          row.names = TRUE)

# Summary of unadjusted results
cat("\n========== UNADJUSTED ANALYSIS SUMMARY ==========\n")
cat("\nNumber of metabolites with FDR < 0.05:\n")
cat("  Overall ANOVA (any group difference):",
    sum(results_ckm_anova$adj.P.Val < 0.05), "\n")
cat("  Stage 2+ vs Stage 2:",
    sum(results_ckm_stage2plus$adj.P.Val < 0.05), "\n")
cat("  Stage 3 vs Stage 2:",
    sum(results_ckm_stage3$adj.P.Val < 0.05), "\n")

cat("\nTop 10 metabolites - Overall ANOVA:\n")
print(head(results_ckm_anova, 10))
cat("\nTop 10 metabolites - Stage 2+ vs Stage 2:\n")
print(head(results_ckm_stage2plus, 10))
cat("\nTop 10 metabolites - Stage 3 vs Stage 2:\n")
print(head(results_ckm_stage3, 10))

#####################
# Adjusted Analysis #
#####################

# keeping observations with nonmissing values covariates
cleaned_baseline_plasma_data_adj <- cleaned_baseline_plasma_data |>
  filter(!is.na(hb_a1c) & !is.na(log_trig) & !is.na(sbp) & !is.na(si_1_ins0)) |>
  mutate(ckm_syn_base_factor = factor(ckm_syn_base,
                                      levels = c("Stage 2",
                                                 "Stage 2+",
                                                 "Stage 3")))

# log transform the adjusted dataset
seq_plasma_temp <- is_seq_plasma(names(cleaned_baseline_plasma_data_adj))
cleaned_log_baseline_plasma_data_adj <- cleaned_baseline_plasma_data_adj |>
  mutate(across(all_of(names(cleaned_baseline_plasma_data_adj)[seq_plasma_temp]),
                ~log(.x)))

# make design matrix (categorical CKM stages + covariates)
design_adj <- model.matrix(~ ckm_syn_base_factor +
                             hb_a1c +
                             log_trig +
                             sbp +
                             si_1_ins0,
                           data = cleaned_baseline_plasma_data_adj)

# create adjusted analysis matrices
seq_plasma_adj <- is_seq_plasma(names(cleaned_log_baseline_plasma_data_adj))
ymat_adj <- t(cleaned_log_baseline_plasma_data_adj[, seq_plasma_adj])

# fit adjusted model
fit_adj <- lmFit(ymat_adj, design_adj)
fit_adj <- eBayes(fit_adj)

# Test for overall group differences - adjusted
results_ckm_adj_anova <- topTable(fit_adj, coef = 2:3, number = nrow(ymat_adj))
results_ckm_adj_anova <- results_ckm_adj_anova[order(results_ckm_adj_anova$P.Value), ]

# Individual contrasts - adjusted
# Stage 2+ vs Stage 2
results_ckm_adj_stage2plus <- topTable(fit_adj, coef = 2, number = nrow(ymat_adj))
# Stage 3 vs Stage 2
results_ckm_adj_stage3 <- topTable(fit_adj, coef = 3, number = nrow(ymat_adj))

# Save adjusted results
write.csv(results_ckm_adj_anova, "results/baseline_adj_anova.csv",
          row.names = TRUE)
write.csv(results_ckm_adj_stage2plus, "results/baseline_adj_stage2plus.csv",
          row.names = TRUE)
write.csv(results_ckm_adj_stage3, "results/baseline_adj_stage3.csv",
          row.names = TRUE)

# Summary of adjusted results
cat("\n========== ADJUSTED ANALYSIS SUMMARY ==========\n")
cat("\nNumber of metabolites with FDR < 0.05:\n")
cat("  Overall ANOVA (any group difference):",
    sum(results_ckm_adj_anova$adj.P.Val < 0.05), "\n")
cat("  Stage 2+ vs Stage 2:",
    sum(results_ckm_adj_stage2plus$adj.P.Val < 0.05), "\n")
cat("  Stage 3 vs Stage 2:",
    sum(results_ckm_adj_stage3$adj.P.Val < 0.05), "\n")

cat("\nTop 10 metabolites - Overall ANOVA (adjusted):\n")
print(head(results_ckm_adj_anova, 10))
cat("\nTop 10 metabolites - Stage 2+ vs Stage 2 (adjusted):\n")
print(head(results_ckm_adj_stage2plus, 10))
cat("\nTop 10 metabolites - Stage 3 vs Stage 2 (adjusted):\n")
print(head(results_ckm_adj_stage3, 10))

# Diagnostic plots for unadjusted analysis
png("results/baseline_unadj_diagnostic_plots.png", width = 1200,
    height = 1000, res = 120)
par(mfrow = c(2, 2), mar = c(5, 4, 4, 2) + 0.5, xpd = FALSE)

# Volcano plot - Stage 2+ vs Stage 2
plot(results_ckm_stage2plus$logFC,
     -log10(results_ckm_stage2plus$P.Value),
     xlab = "Log2 Fold Change", ylab = "-log10(p-value)",
     main = "Volcano: Stage 2+ vs Stage 2 (Unadjusted)",
     pch = 20,
     col = ifelse(results_ckm_stage2plus$adj.P.Val < 0.05, "red", "gray"))
abline(h = -log10(0.05), col = "blue", lty = 2)
top_metabolites <- rownames(results_ckm_stage2plus)[1:5]
top_metabolites <- gsub("\\.in\\.uM.*", "", top_metabolites)
text(results_ckm_stage2plus$logFC[1:5],
     -log10(results_ckm_stage2plus$P.Value[1:5]),
     labels = top_metabolites, cex = 0.6, pos = 4, offset = 0.3)

# P-value histogram
hist(results_ckm_stage2plus$P.Value, breaks = 50,
     main = "P-value Distribution: Stage 2+ vs Stage 2",
     xlab = "P-value", col = "lightblue")

# Volcano plot - Stage 3 vs Stage 2
plot(results_ckm_stage3$logFC, -log10(results_ckm_stage3$P.Value),
     xlab = "Log2 Fold Change", ylab = "-log10(p-value)",
     main = "Volcano: Stage 3 vs Stage 2 (Unadjusted)",
     pch = 20,
     col = ifelse(results_ckm_stage3$adj.P.Val < 0.05, "red", "gray"))
abline(h = -log10(0.05), col = "blue", lty = 2)
top_metabolites <- rownames(results_ckm_stage3)[1:5]
top_metabolites <- gsub("\\.in\\.uM.*", "", top_metabolites)
text(results_ckm_stage3$logFC[1:5],
     -log10(results_ckm_stage3$P.Value[1:5]),
     labels = top_metabolites, cex = 0.6, pos = 4, offset = 0.3)

# QQ plot
qqnorm(fit$t[, 2], main = "QQ Plot: Stage 2+ vs Stage 2")
qqline(fit$t[, 2])

par(mfrow = c(1, 1))
dev.off()

# Diagnostic plots for adjusted analysis
png("results/baseline_adj_diagnostic_plots.png", width = 1200,
    height = 1000, res = 120)
par(mfrow = c(2, 2), mar = c(5, 4, 4, 2) + 0.5, xpd = FALSE)

# Volcano plot - Stage 2+ vs Stage 2 (adjusted)
plot(results_ckm_adj_stage2plus$logFC,
     -log10(results_ckm_adj_stage2plus$P.Value),
     xlab = "Log2 Fold Change", ylab = "-log10(p-value)",
     main = "Volcano: Stage 2+ vs Stage 2 (Adjusted)",
     pch = 20,
     col = ifelse(results_ckm_adj_stage2plus$adj.P.Val < 0.05, "red", "gray"))
abline(h = -log10(0.05), col = "blue", lty = 2)
top_metabolites <- rownames(results_ckm_adj_stage2plus)[1:5]
top_metabolites <- gsub("\\.in\\.uM.*", "", top_metabolites)
text(results_ckm_adj_stage2plus$logFC[1:5],
     -log10(results_ckm_adj_stage2plus$P.Value[1:5]),
     labels = top_metabolites, cex = 0.6, pos = 4, offset = 0.3)

# P-value histogram
hist(results_ckm_adj_stage2plus$P.Value, breaks = 50,
     main = "P-value Distribution: Stage 2+ vs Stage 2 (Adj)",
     xlab = "P-value", col = "lightblue")

# Volcano plot - Stage 3 vs Stage 2 (adjusted)
plot(results_ckm_adj_stage3$logFC,
     -log10(results_ckm_adj_stage3$P.Value),
     xlab = "Log2 Fold Change", ylab = "-log10(p-value)",
     main = "Volcano: Stage 3 vs Stage 2 (Adjusted)",
     pch = 20,
     col = ifelse(results_ckm_adj_stage3$adj.P.Val < 0.05, "red", "gray"))
abline(h = -log10(0.05), col = "blue", lty = 2)
top_metabolites <- rownames(results_ckm_adj_stage3)[1:5]
top_metabolites <- gsub("\\.in\\.uM.*", "", top_metabolites)
text(results_ckm_adj_stage3$logFC[1:5],
     -log10(results_ckm_adj_stage3$P.Value[1:5]),
     labels = top_metabolites, cex = 0.6, pos = 4, offset = 0.3)

# QQ plot
qqnorm(fit_adj$t[, 2], main = "QQ Plot: Stage 2+ vs Stage 2 (Adj)")
qqline(fit_adj$t[, 2])

par(mfrow = c(1, 1))
dev.off()

cat("\n========== PLASMA FILES SAVED ==========\n")
cat("Results files:\n")
cat("  - baseline_unadj_anova.csv\n")
cat("  - baseline_unadj_stage2plus.csv\n")
cat("  - baseline_unadj_stage3.csv\n")
cat("  - baseline_adj_anova.csv\n")
cat("  - baseline_adj_stage2plus.csv\n")
cat("  - baseline_adj_stage3.csv\n")
cat("Diagnostic plots:\n")
cat("  - baseline_unadj_diagnostic_plots.png\n")
cat("  - baseline_adj_diagnostic_plots.png\n")

#### Urine Baseline Analysis using limma ####

#######################
# Unadjusted Analysis #
#######################

# Identify urine metabolite columns
seq_urine <- is_seq_urine(names(cleaned_log_baseline_urine_data))

# create design matrix for unadjusted analysis (ANOVA-style)
design_unadj_urine <- model.matrix(~ ckm_syn_base_factor,
                                   data = cleaned_log_baseline_urine_data)

# prepare response matrix
ymat_urine <- t(cleaned_log_baseline_urine_data[, seq_urine])

# fit unadjusted models
fit_urine <- lmFit(ymat_urine, design_unadj_urine)
fit_urine <- eBayes(fit_urine)

# Test for overall group differences (ANOVA F-test)
results_ckm_anova_urine <- topTable(fit_urine, coef = 2:3,
                                    number = nrow(ymat_urine))
results_ckm_anova_urine <- results_ckm_anova_urine[
  order(results_ckm_anova_urine$P.Value), ]

# Individual contrasts (pairwise comparisons)
# Stage 2+ vs Stage 2
results_ckm_stage2plus_urine <- topTable(fit_urine, coef = 2,
                                         number = nrow(ymat_urine))
# Stage 3 vs Stage 2
results_ckm_stage3_urine <- topTable(fit_urine, coef = 3,
                                     number = nrow(ymat_urine))

# Save unadjusted results
write.csv(results_ckm_anova_urine,
          "results/baseline_urine_unadj_anova.csv", row.names = TRUE)
write.csv(results_ckm_stage2plus_urine,
          "results/baseline_urine_unadj_stage2plus.csv", row.names = TRUE)
write.csv(results_ckm_stage3_urine,
          "results/baseline_urine_unadj_stage3.csv", row.names = TRUE)

# Summary of unadjusted results
cat("\n========== URINE UNADJUSTED ANALYSIS SUMMARY ==========\n")
cat("\nNumber of metabolites with FDR < 0.05:\n")
cat("  Overall ANOVA (any group difference):",
    sum(results_ckm_anova_urine$adj.P.Val < 0.05), "\n")
cat("  Stage 2+ vs Stage 2:",
    sum(results_ckm_stage2plus_urine$adj.P.Val < 0.05), "\n")
cat("  Stage 3 vs Stage 2:",
    sum(results_ckm_stage3_urine$adj.P.Val < 0.05), "\n")

cat("\nTop 10 metabolites - Overall ANOVA:\n")
print(head(results_ckm_anova_urine, 10))
cat("\nTop 10 metabolites - Stage 2+ vs Stage 2:\n")
print(head(results_ckm_stage2plus_urine, 10))
cat("\nTop 10 metabolites - Stage 3 vs Stage 2:\n")
print(head(results_ckm_stage3_urine, 10))

#####################
# Adjusted Analysis #
#####################

# keeping observations with nonmissing values covariates
cleaned_baseline_urine_data_adj <- cleaned_baseline_urine_data |>
  filter(!is.na(hb_a1c) & !is.na(log_trig) &
         !is.na(sbp) & !is.na(si_1_ins0)) |>
  mutate(ckm_syn_base_factor = factor(ckm_syn_base,
                                      levels = c("Stage 2",
                                                 "Stage 2+",
                                                 "Stage 3")))

# log transform the adjusted dataset
seq_urine_temp <- is_seq_urine(names(cleaned_baseline_urine_data_adj))
cleaned_log_baseline_urine_data_adj <- cleaned_baseline_urine_data_adj |>
  mutate(across(all_of(names(cleaned_baseline_urine_data_adj)[seq_urine_temp]),
                ~log(.x)))

# make design matrix (categorical CKM stages + covariates)
design_adj_urine <- model.matrix(~ ckm_syn_base_factor +
                                  hb_a1c +
                                  log_trig +
                                  sbp +
                                  si_1_ins0,
                                data = cleaned_log_baseline_urine_data_adj)

# create adjusted analysis matrices
seq_urine_adj <- is_seq_urine(names(cleaned_log_baseline_urine_data_adj))
ymat_adj_urine <- t(cleaned_log_baseline_urine_data_adj[, seq_urine_adj])

# fit adjusted model
fit_adj_urine <- lmFit(ymat_adj_urine, design_adj_urine)
fit_adj_urine <- eBayes(fit_adj_urine)

# Test for overall group differences - adjusted
results_ckm_adj_anova_urine <- topTable(fit_adj_urine, coef = 2:3,
                                        number = nrow(ymat_adj_urine))
results_ckm_adj_anova_urine <- results_ckm_adj_anova_urine[
  order(results_ckm_adj_anova_urine$P.Value), ]

# Individual contrasts - adjusted
# Stage 2+ vs Stage 2
results_ckm_adj_stage2plus_urine <- topTable(fit_adj_urine, coef = 2,
                                             number = nrow(ymat_adj_urine))
# Stage 3 vs Stage 2
results_ckm_adj_stage3_urine <- topTable(fit_adj_urine, coef = 3,
                                         number = nrow(ymat_adj_urine))

# Save adjusted results
write.csv(results_ckm_adj_anova_urine,
          "results/baseline_urine_adj_anova.csv", row.names = TRUE)
write.csv(results_ckm_adj_stage2plus_urine,
          "results/baseline_urine_adj_stage2plus.csv", row.names = TRUE)
write.csv(results_ckm_adj_stage3_urine,
          "results/baseline_urine_adj_stage3.csv", row.names = TRUE)

# Summary of adjusted results
cat("\n========== URINE ADJUSTED ANALYSIS SUMMARY ==========\n")
cat("\nNumber of metabolites with FDR < 0.05:\n")
cat("  Overall ANOVA (any group difference):",
    sum(results_ckm_adj_anova_urine$adj.P.Val < 0.05), "\n")
cat("  Stage 2+ vs Stage 2:",
    sum(results_ckm_adj_stage2plus_urine$adj.P.Val < 0.05), "\n")
cat("  Stage 3 vs Stage 2:",
    sum(results_ckm_adj_stage3_urine$adj.P.Val < 0.05), "\n")

cat("\nTop 10 metabolites - Overall ANOVA (adjusted):\n")
print(head(results_ckm_adj_anova_urine, 10))
cat("\nTop 10 metabolites - Stage 2+ vs Stage 2 (adjusted):\n")
print(head(results_ckm_adj_stage2plus_urine, 10))
cat("\nTop 10 metabolites - Stage 3 vs Stage 2 (adjusted):\n")
print(head(results_ckm_adj_stage3_urine, 10))

# Diagnostic plots for unadjusted analysis
png("results/baseline_urine_unadj_diagnostic_plots.png", width = 1200,
    height = 1000, res = 120)
par(mfrow = c(2, 2), mar = c(5, 4, 4, 2) + 0.5, xpd = FALSE)

# Volcano plot - Stage 2+ vs Stage 2
plot(results_ckm_stage2plus_urine$logFC,
     -log10(results_ckm_stage2plus_urine$P.Value),
     xlab = "Log2 Fold Change", ylab = "-log10(p-value)",
     main = "Volcano: Stage 2+ vs Stage 2 (Urine, Unadj)",
     pch = 20,
     col = ifelse(results_ckm_stage2plus_urine$adj.P.Val < 0.05,
                  "red", "gray"))
abline(h = -log10(0.05), col = "blue", lty = 2)
top_metabolites <- rownames(results_ckm_stage2plus_urine)[1:5]
top_metabolites <- gsub("\\.in\\..*", "", top_metabolites)
text(results_ckm_stage2plus_urine$logFC[1:5],
     -log10(results_ckm_stage2plus_urine$P.Value[1:5]),
     labels = top_metabolites, cex = 0.6, pos = 4, offset = 0.3)

# P-value histogram
hist(results_ckm_stage2plus_urine$P.Value, breaks = 50,
     main = "P-value Distribution: Stage 2+ vs Stage 2 (Urine)",
     xlab = "P-value", col = "lightblue")

# Volcano plot - Stage 3 vs Stage 2
plot(results_ckm_stage3_urine$logFC,
     -log10(results_ckm_stage3_urine$P.Value),
     xlab = "Log2 Fold Change", ylab = "-log10(p-value)",
     main = "Volcano: Stage 3 vs Stage 2 (Urine, Unadj)",
     pch = 20,
     col = ifelse(results_ckm_stage3_urine$adj.P.Val < 0.05,
                  "red", "gray"))
abline(h = -log10(0.05), col = "blue", lty = 2)
top_metabolites <- rownames(results_ckm_stage3_urine)[1:5]
top_metabolites <- gsub("\\.in\\..*", "", top_metabolites)
text(results_ckm_stage3_urine$logFC[1:5],
     -log10(results_ckm_stage3_urine$P.Value[1:5]),
     labels = top_metabolites, cex = 0.6, pos = 4, offset = 0.3)

# QQ plot
qqnorm(fit_urine$t[, 2], main = "QQ Plot: Stage 2+ vs Stage 2 (Urine)")
qqline(fit_urine$t[, 2])

par(mfrow = c(1, 1))
dev.off()

# Diagnostic plots for adjusted analysis
png("results/baseline_urine_adj_diagnostic_plots.png", width = 1200,
    height = 1000, res = 120)
par(mfrow = c(2, 2), mar = c(5, 4, 4, 2) + 0.5, xpd = FALSE)

# Volcano plot - Stage 2+ vs Stage 2 (adjusted)
plot(results_ckm_adj_stage2plus_urine$logFC,
     -log10(results_ckm_adj_stage2plus_urine$P.Value),
     xlab = "Log2 Fold Change", ylab = "-log10(p-value)",
     main = "Volcano: Stage 2+ vs Stage 2 (Urine, Adj)",
     pch = 20,
     col = ifelse(results_ckm_adj_stage2plus_urine$adj.P.Val < 0.05,
                  "red", "gray"))
abline(h = -log10(0.05), col = "blue", lty = 2)
top_metabolites <- rownames(results_ckm_adj_stage2plus_urine)[1:5]
top_metabolites <- gsub("\\.in\\..*", "", top_metabolites)
text(results_ckm_adj_stage2plus_urine$logFC[1:5],
     -log10(results_ckm_adj_stage2plus_urine$P.Value[1:5]),
     labels = top_metabolites, cex = 0.6, pos = 4, offset = 0.3)

# P-value histogram
hist(results_ckm_adj_stage2plus_urine$P.Value, breaks = 50,
     main = "P-value Distribution: Stage 2+ vs Stage 2 (Urine, Adj)",
     xlab = "P-value", col = "lightblue")

# Volcano plot - Stage 3 vs Stage 2 (adjusted)
plot(results_ckm_adj_stage3_urine$logFC,
     -log10(results_ckm_adj_stage3_urine$P.Value),
     xlab = "Log2 Fold Change", ylab = "-log10(p-value)",
     main = "Volcano: Stage 3 vs Stage 2 (Urine, Adj)",
     pch = 20,
     col = ifelse(results_ckm_adj_stage3_urine$adj.P.Val < 0.05,
                  "red", "gray"))
abline(h = -log10(0.05), col = "blue", lty = 2)
top_metabolites <- rownames(results_ckm_adj_stage3_urine)[1:5]
top_metabolites <- gsub("\\.in\\..*", "", top_metabolites)
text(results_ckm_adj_stage3_urine$logFC[1:5],
     -log10(results_ckm_adj_stage3_urine$P.Value[1:5]),
     labels = top_metabolites, cex = 0.6, pos = 4, offset = 0.3)

# QQ plot
qqnorm(fit_adj_urine$t[, 2],
       main = "QQ Plot: Stage 2+ vs Stage 2 (Urine, Adj)")
qqline(fit_adj_urine$t[, 2])

par(mfrow = c(1, 1))
dev.off()

cat("\n========== URINE FILES SAVED ==========\n")
cat("Results files:\n")
cat("  - baseline_urine_unadj_anova.csv\n")
cat("  - baseline_urine_unadj_stage2plus.csv\n")
cat("  - baseline_urine_unadj_stage3.csv\n")
cat("  - baseline_urine_adj_anova.csv\n")
cat("  - baseline_urine_adj_stage2plus.csv\n")
cat("  - baseline_urine_adj_stage3.csv\n")
cat("Diagnostic plots:\n")
cat("  - baseline_urine_unadj_diagnostic_plots.png\n")
cat("  - baseline_urine_adj_diagnostic_plots.png\n")
