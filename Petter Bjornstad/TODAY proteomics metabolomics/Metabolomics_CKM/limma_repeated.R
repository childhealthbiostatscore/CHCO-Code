# ---
# Mixed Model Analysis using limma for CKM Progression (unadjusted and adjusted)
# Author: Darwin Del Castillo
# Date: `r lubridate::today()`
# ---

pacman::p_load(limma)

###############
### Plasma ####
###############

# Check data structure
cat("Data structure:\n")
str(plasma_long)
cat("\nNumber of unique subjects:", length(unique(plasma_long$releaseid)), "\n")
cat("Number of observations:", nrow(plasma_long), "\n")

# Identify all metabolite columns
metabolite_cols <- str_subset(names(plasma_long), ".in.uM")
cat("\nNumber of metabolites found:", length(metabolite_cols), "\n")
cat("First few metabolite names:", head(metabolite_cols), "\n")

# Extract just the metabolite data
metabolite_data <- plasma_long[, c("releaseid", "visit", metabolite_cols)]

# Transpose so metabolites are rows and samples are columns
expression_matrix <- t(metabolite_data[, -c(1, 2)])  # Exclude releaseid and visit for transposition
colnames(expression_matrix) <- paste(metabolite_data$releaseid, metabolite_data$visit, sep = "_")

cat("\nExpression matrix dimensions: ", dim(expression_matrix), "\n")
cat("Rows (metabolites):", nrow(expression_matrix), "\n")
cat("Columns (samples):", ncol(expression_matrix), "\n")

# Create metadata dataframe
metadata <- plasma_long[, c("releaseid", "visit", "progress_ckm")] |>
  as.data.frame() |>
  distinct(releaseid, visit, .keep_all = TRUE)
rownames(metadata) <- paste(metadata$releaseid, metadata$visit, sep = "_")

# Verify alignment
# if(!all(colnames(expression_matrix) == rownames(metadata))) {
#   stop("Sample order mismatch between expression matrix and metadata!")
# }

cat("\nFirst few samples in expression matrix:\n")
print(expression_matrix[1:min(3, nrow(expression_matrix)), 1:min(6, ncol(expression_matrix))])

# Create design matrix
metadata$visit_factor <- factor(metadata$visit)
metadata$progress_ckm <- factor(metadata$progress_ckm, levels = c(0, 1))

# Design with interaction term
design <- model.matrix(~ progress_ckm * visit_factor, data = metadata)
colnames(design) <- c("Intercept", "ProgressCKM", "Visit2", "ProgressCKM_Visit2")

# Alternative simpler model without interaction:
#design <- model.matrix(~ progress_CKM + visit_factor, data = metadata)
#colnames(design) <- c("Intercept", "ProgressCKM", "Visit2")

cat("\nDesign matrix dimensions:", dim(design), "\n")
cat("Design matrix (first few rows):\n")
print(head(design))

# Set up blocking factor for repeated measures
block <- factor(metadata$releaseid)
cat("\nBlocking structure (subject IDs):\n")
print(table(block))

# Estimate the correlation between repeated measures
cat("\nEstimating intra-subject correlation across all metabolites...\n")
correlation <- duplicateCorrelation(expression_matrix, design, block = block)
cat("Consensus correlation:", round(correlation$consensus.correlation, 4), "\n")

# Fit the mixed model with blocking and correlation structure
cat("\nFitting mixed model for", nrow(expression_matrix), "metabolites...\n")
fit <- lmFit(expression_matrix, design, block = block,
             correlation = correlation$consensus.correlation)

# Apply empirical Bayes moderation
fit <- eBayes(fit)

# Extract results for different contrasts
cat("\n========== RESULTS ==========\n")

# 1. Main effect of CKM progression (averaged across visits)
cat("\n1. Main effect of CKM progression (top 20 proteins):\n")
results_progression <- topTable(fit, coef = "ProgressCKM", adjust.method = "BH", number = Inf)
print(head(results_progression, 20))

# Save full results
write.csv(results_progression, "results/plasma_limma_results_progression.csv", row.names = TRUE)
cat("\nFull results saved to 'plasma_limma_results_progression.csv'\n")

# 2. Visit effect (change from visit 1 to visit 2)
cat("\n2. Visit effect - Visit 2 vs Visit 1 (top 20 proteins):\n")
results_visit <- topTable(fit, coef = "Visit2", adjust.method = "BH", number = Inf)
print(head(results_visit, 20))
write.csv(results_visit, "results/plasma_limma_results_visit.csv", row.names = TRUE)

# 3. Interaction effect
if("ProgressCKM_Visit2" %in% colnames(design)) {
  cat("\n3. Interaction - differential progression effect by visit (top 20):\n")
  results_interaction <- topTable(fit, coef = "ProgressCKM_Visit2", adjust.method = "BH", number = Inf)
  print(head(results_interaction, 20))
  write.csv(results_interaction, "results/plasma_limma_results_interaction.csv", row.names = TRUE)
}

# Create custom contrasts for specific comparisons
contrast.matrix <- makeContrasts(
  Progression_Visit1 = ProgressCKM,  # Effect at visit 1
  Progression_Visit2 = ProgressCKM + ProgressCKM_Visit2,  # Effect at visit 2
  Progression_Average = ProgressCKM + 0.5*ProgressCKM_Visit2,  # Average effect
  Change_in_Effect = ProgressCKM_Visit2,  # How effect changes over time
  levels = design
)

# Fit contrasts
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

cat("\n4. Custom contrasts:\n")
cat("\n4a. CKM progression effect at Visit 1 (top 20):\n")
results_v1 <- topTable(fit2, coef = "Progression_Visit1", adjust.method = "BH", number = Inf)
print(head(results_v1, 20))

cat("\n4b. CKM progression effect at Visit 2 (top 20):\n")
results_v2 <- topTable(fit2, coef = "Progression_Visit2", adjust.method = "BH", number = Inf)
print(head(results_v2, 20))

# Summary of significant findings
cat("\n========== SUMMARY OF SIGNIFICANT FINDINGS ==========\n")
cat("\nNumber of metabolites with FDR < 0.05:\n")
cat("  Main progression effect:", sum(results_progression$adj.P.Val < 0.05), "\n")
cat("  Visit effect:", sum(results_visit$adj.P.Val < 0.05), "\n")
if("ProgressCKM_Visit2" %in% colnames(design)) {
  cat("  Interaction effect:", sum(results_interaction$adj.P.Val < 0.05), "\n")
}

# Diagnostic plots for plasma
png("results/plasma_diagnostic_plots.png", width = 1200, height = 1000,
    res = 120)
par(mfrow = c(2, 2), mar = c(5, 4, 4, 2) + 0.5, xpd = FALSE)

# 1. Volcano plot - Progression effect
plot(results_progression$logFC, -log10(results_progression$P.Value),
     xlab = "Log2 Fold Change", ylab = "-log10(p-value)",
     main = "Volcano Plot: CKM Progression Effect (Plasma)",
     pch = 20, col = ifelse(results_progression$adj.P.Val < 0.05, "red", "gray"))
abline(h = -log10(0.05), col = "blue", lty = 2)
# Add labels with intelligent positioning to avoid cutoff
top_metabolites <- rownames(results_progression)[1:5]
top_metabolites <- gsub("\\.in\\.uM.*", "", top_metabolites)  # Shorten names
text(results_progression$logFC[1:5], -log10(results_progression$P.Value[1:5]),
     labels = top_metabolites, cex = 0.6, pos = 4, offset = 0.3)

# 2. P-value histogram
hist(results_progression$P.Value, breaks = 50,
     main = "P-value Distribution: Progression Effect (Plasma)",
     xlab = "P-value", col = "lightblue")

# 3. MA plot
plotMA(fit, coef = "ProgressCKM", main = "MA Plot: CKM Progression (Plasma)")

# 4. QQ plot
qqnorm(fit$t[, "ProgressCKM"], main = "QQ Plot: Test Statistics (Plasma)")
qqline(fit$t[, "ProgressCKM"])

par(mfrow = c(1, 1))
dev.off()

cat("\nAnalysis complete!\n")
cat("Results files saved:\n")
cat("  - plasma_limma_results_progression.csv\n")
cat("  - plasma_limma_results_visit.csv\n")
cat("  - plasma_limma_results_interaction.csv\n")
cat("  - plasma_diagnostic_plots.png\n")

#############
### Urine ###
#############

# Check data structure
cat("Data structure:\n")
str(urine_long)
cat("\nNumber of unique subjects:", length(unique(urine_long$releaseid)), "\n")
cat("Number of observations:", nrow(urine_long), "\n")

# Identify all metabolite columns (uM/mM.Creatinine only)
metabolite_cols <- str_subset(names(urine_long),
                              "\\.in\\.uM/mM\\.Creatinine")
cat("\nNumber of metabolites found:", length(metabolite_cols), "\n")
cat("First few metabolite names:", head(metabolite_cols), "\n")

# Extract just the metabolite data
metabolite_data <- urine_long[, c("releaseid", "visit", metabolite_cols)]

# Transpose so metabolites are rows and samples are columns
expression_matrix <- t(metabolite_data[, -c(1, 2)])  # Exclude releaseid and visit for transposition
colnames(expression_matrix) <- paste(metabolite_data$releaseid, metabolite_data$visit, sep = "_")

cat("\nExpression matrix dimensions: ", dim(expression_matrix), "\n")
cat("Rows (metabolites):", nrow(expression_matrix), "\n")
cat("Columns (samples):", ncol(expression_matrix), "\n")

# Create metadata dataframe
metadata <- urine_long[, c("releaseid", "visit", "progress_ckm")] |>
  as.data.frame() |>
  distinct(releaseid, visit, .keep_all = TRUE)
rownames(metadata) <- paste(metadata$releaseid, metadata$visit, sep = "_")

cat("\nFirst few samples in expression matrix:\n")
print(expression_matrix[1:min(3, nrow(expression_matrix)), 1:min(6, ncol(expression_matrix))])

# Create design matrix
metadata$visit_factor <- factor(metadata$visit)
metadata$progress_ckm <- factor(metadata$progress_ckm, levels = c(0, 1))

# Design with interaction term
design <- model.matrix(~ progress_ckm * visit_factor, data = metadata)
colnames(design) <- c("Intercept", "ProgressCKM", "Visit2", "ProgressCKM_Visit2")

# Alternative simpler model without interaction:
#design <- model.matrix(~ progress_CKM + visit_factor, data = metadata)
#colnames(design) <- c("Intercept", "ProgressCKM", "Visit2")

cat("\nDesign matrix dimensions:", dim(design), "\n")
cat("Design matrix (first few rows):\n")
print(head(design))

# Set up blocking factor for repeated measures
block <- factor(metadata$releaseid)
cat("\nBlocking structure (subject IDs):\n")
print(table(block))

# Estimate the correlation between repeated measures
cat("\nEstimating intra-subject correlation across all metabolites...\n")
correlation <- duplicateCorrelation(expression_matrix, design, block = block)
cat("Consensus correlation:", round(correlation$consensus.correlation, 4), "\n")

# Fit the mixed model with blocking and correlation structure
cat("\nFitting mixed model for", nrow(expression_matrix), "metabolites...\n")
fit <- lmFit(expression_matrix, design, block = block,
             correlation = correlation$consensus.correlation)

# Apply empirical Bayes moderation
fit <- eBayes(fit)

# Extract results for different contrasts
cat("\n========== RESULTS ==========\n")

# 1. Main effect of CKM progression (averaged across visits)
cat("\n1. Main effect of CKM progression (top 20 proteins):\n")
results_progression <- topTable(fit, coef = "ProgressCKM", adjust.method = "BH", number = Inf)
print(head(results_progression, 20))

# Save full results
write.csv(results_progression, "results/urine_limma_results_progression.csv", row.names = TRUE)
cat("\nFull results saved to 'urine_limma_results_progression.csv'\n")

# 2. Visit effect (change from visit 1 to visit 2)
cat("\n2. Visit effect - Visit 2 vs Visit 1 (top 20 proteins):\n")
results_visit <- topTable(fit, coef = "Visit2", adjust.method = "BH", number = Inf)
print(head(results_visit, 20))
write.csv(results_visit, "results/urine_limma_results_visit.csv", row.names = TRUE)

# 3. Interaction effect (if included in model)
if("ProgressCKM_Visit2" %in% colnames(design)) {
  cat("\n3. Interaction - differential progression effect by visit (top 20):\n")
  results_interaction <- topTable(fit, coef = "ProgressCKM_Visit2", adjust.method = "BH", number = Inf)
  print(head(results_interaction, 20))
  write.csv(results_interaction, "results/urine_limma_results_interaction.csv", row.names = TRUE)
}

# Create custom contrasts for specific comparisons
contrast.matrix <- makeContrasts(
  Progression_Visit1 = ProgressCKM,  # Effect at visit 1
  Progression_Visit2 = ProgressCKM + ProgressCKM_Visit2,  # Effect at visit 2
  Progression_Average = ProgressCKM + 0.5*ProgressCKM_Visit2,  # Average effect
  Change_in_Effect = ProgressCKM_Visit2,  # How effect changes over time
  levels = design
)

# Fit contrasts
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

cat("\n4. Custom contrasts:\n")
cat("\n4a. CKM progression effect at Visit 1 (top 20):\n")
results_v1 <- topTable(fit2, coef = "Progression_Visit1", adjust.method = "BH", number = Inf)
print(head(results_v1, 20))

cat("\n4b. CKM progression effect at Visit 2 (top 20):\n")
results_v2 <- topTable(fit2, coef = "Progression_Visit2", adjust.method = "BH", number = Inf)
print(head(results_v2, 20))

# Summary of significant findings
cat("\n========== SUMMARY OF SIGNIFICANT FINDINGS ==========\n")
cat("\nNumber of metabolites with FDR < 0.05:\n")
cat("  Main progression effect:", sum(results_progression$adj.P.Val < 0.05), "\n")
cat("  Visit effect:", sum(results_visit$adj.P.Val < 0.05), "\n")
if("ProgressCKM_Visit2" %in% colnames(design)) {
  cat("  Interaction effect:", sum(results_interaction$adj.P.Val < 0.05), "\n")
}

# Diagnostic plots for urine
png("results/urine_diagnostic_plots.png", width = 1200, height = 1000,
    res = 120)
par(mfrow = c(2, 2), mar = c(5, 4, 4, 2) + 0.5, xpd = FALSE)

# 1. Volcano plot - Progression effect
plot(results_progression$logFC, -log10(results_progression$P.Value),
     xlab = "Log2 Fold Change", ylab = "-log10(p-value)",
     main = "Volcano Plot: CKM Progression Effect (Urine)",
     pch = 20,
     col = ifelse(results_progression$adj.P.Val < 0.05, "red", "gray"))
abline(h = -log10(0.05), col = "blue", lty = 2)
# Add labels with intelligent positioning to avoid cutoff
top_metabolites <- rownames(results_progression)[1:5]
top_metabolites <- gsub("\\.in\\..*", "", top_metabolites)  # Shorten names
text(results_progression$logFC[1:5], -log10(results_progression$P.Value[1:5]),
     labels = top_metabolites, cex = 0.6, pos = 4, offset = 0.3)

# 2. P-value histogram
hist(results_progression$P.Value, breaks = 50,
     main = "P-value Distribution: Progression Effect (Urine)",
     xlab = "P-value", col = "lightblue")

# 3. MA plot
plotMA(fit, coef = "ProgressCKM", main = "MA Plot: CKM Progression (Urine)")

# 4. QQ plot
qqnorm(fit$t[, "ProgressCKM"], main = "QQ Plot: Test Statistics (Urine)")
qqline(fit$t[, "ProgressCKM"])

par(mfrow = c(1, 1))
dev.off()

cat("\nAnalysis complete!\n")
cat("Results files saved:\n")
cat("  - urine_limma_results_progression.csv\n")
cat("  - urine_limma_results_visit.csv\n")
cat("  - urine_limma_results_interaction.csv\n")
cat("  - urine_diagnostic_plots.png\n")