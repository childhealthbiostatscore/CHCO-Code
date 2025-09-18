# Mixed Model Analysis using limma for CKM Syndrome Progression
# Analyzing multiple proteins with repeated measures (visits 1 and 2)

# Load required libraries
library(limma)
library(tidyverse)

# Read the data (adjust path as needed)
data <- soma_long

# Example data structure with multiple proteins:
# data <- data.frame(
#   RELEASEID = c("65-10210", "65-10210", "65-10368", "65-10368", "65-10550", "65-10550"),
#   visit = c(1, 2, 1, 2, 1, 2),
#   progress_CKM = c(0, 0, 0, 0, 1, 1),
#   seq.10000.28 = c(501.3, 500.6, 483.2, 540.1, 526.0, 495.6),
#   seq.10000.29 = c(234.5, 245.6, 223.4, 256.7, 267.8, 234.5),
#   seq.10000.30 = c(1023.4, 1045.6, 987.6, 1034.5, 1123.4, 1098.7)
# )

# Check data structure
cat("Data structure:\n")
str(data)
cat("\nNumber of unique subjects:", length(unique(data$RELEASEID)), "\n")
cat("Number of observations:", nrow(data), "\n")

# Identify all protein columns
protein_cols <- grep("^seq\\.", names(data), value = TRUE)
cat("\nNumber of proteins found:", length(protein_cols), "\n")
cat("First few protein names:", head(protein_cols), "\n")

# Create unique sample identifiers
data$sample_id <- paste(data$RELEASEID, "v", data$visit, sep = "")

# Method 1: Reshape to wide format for limma (proteins in rows, samples in columns)
# Extract just the protein data
protein_data <- data[, c("sample_id", protein_cols)]

# Transpose so proteins are rows and samples are columns
expression_matrix <- t(protein_data[, -1])  # Exclude sample_id for transposition
colnames(expression_matrix) <- protein_data$sample_id

cat("\nExpression matrix dimensions: ", dim(expression_matrix), "\n")
cat("Rows (proteins):", nrow(expression_matrix), "\n")
cat("Columns (samples):", ncol(expression_matrix), "\n")

# Create metadata dataframe aligned with expression matrix columns
metadata <- data[, c("sample_id", "RELEASEID", "visit", "progress_CKM")]
metadata <- metadata[!duplicated(metadata$sample_id), ]
rownames(metadata) <- metadata$sample_id

# Ensure metadata order matches expression matrix columns
metadata <- metadata[colnames(expression_matrix), ]

# Verify alignment
if(!all(colnames(expression_matrix) == rownames(metadata))) {
  stop("Sample order mismatch between expression matrix and metadata!")
}

cat("\nFirst few samples in expression matrix:\n")
print(expression_matrix[1:min(3, nrow(expression_matrix)), 1:min(6, ncol(expression_matrix))])

# Create design matrix
metadata$visit_factor <- factor(metadata$visit)
metadata$progress_CKM <- factor(metadata$progress_CKM, levels = c(0, 1))

# Design with interaction term
design <- model.matrix(~ progress_CKM * visit_factor, data = metadata)
colnames(design) <- c("Intercept", "ProgressCKM", "Visit2", "ProgressCKM_Visit2")

# Alternative simpler model without interaction:
#design <- model.matrix(~ progress_CKM + visit_factor, data = metadata)
#colnames(design) <- c("Intercept", "ProgressCKM", "Visit2")

cat("\nDesign matrix dimensions:", dim(design), "\n")
cat("Design matrix (first few rows):\n")
print(head(design))

# Set up blocking factor for repeated measures
block <- factor(metadata$RELEASEID)
cat("\nBlocking structure (subject IDs):\n")
print(table(block))

# Estimate the correlation between repeated measures
cat("\nEstimating intra-subject correlation across all proteins...\n")
correlation <- duplicateCorrelation(expression_matrix, design, block = block)
cat("Consensus correlation:", round(correlation$consensus.correlation, 4), "\n")

# Fit the mixed model with blocking and correlation structure
cat("\nFitting mixed model for", nrow(expression_matrix), "proteins...\n")
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
write.csv(results_progression, "limma_results_progression.csv", row.names = TRUE)
cat("\nFull results saved to 'limma_results_progression.csv'\n")

# 2. Visit effect (change from visit 1 to visit 2)
cat("\n2. Visit effect - Visit 2 vs Visit 1 (top 20 proteins):\n")
results_visit <- topTable(fit, coef = "Visit2", adjust.method = "BH", number = Inf)
print(head(results_visit, 20))
write.csv(results_visit, "limma_results_visit.csv", row.names = TRUE)

# 3. Interaction effect (if included in model)
if("ProgressCKM_Visit2" %in% colnames(design)) {
  cat("\n3. Interaction - differential progression effect by visit (top 20):\n")
  results_interaction <- topTable(fit, coef = "ProgressCKM_Visit2", adjust.method = "BH", number = Inf)
  print(head(results_interaction, 20))
  write.csv(results_interaction, "limma_results_interaction.csv", row.names = TRUE)
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
cat("\nNumber of proteins with FDR < 0.05:\n")
cat("  Main progression effect:", sum(results_progression$adj.P.Val < 0.05), "\n")
cat("  Visit effect:", sum(results_visit$adj.P.Val < 0.05), "\n")
if("ProgressCKM_Visit2" %in% colnames(design)) {
  cat("  Interaction effect:", sum(results_interaction$adj.P.Val < 0.05), "\n")
}

# Volcano plot for progression effect
par(mfrow = c(2, 2))

# 1. Volcano plot - Progression effect
plot(results_progression$logFC, -log10(results_progression$P.Value),
     xlab = "Log2 Fold Change", ylab = "-log10(p-value)",
     main = "Volcano Plot: CKM Progression Effect",
     pch = 20, col = ifelse(results_progression$adj.P.Val < 0.05, "red", "gray"))
abline(h = -log10(0.05), col = "blue", lty = 2)
text(results_progression$logFC[1:5], -log10(results_progression$P.Value[1:5]),
     labels = rownames(results_progression)[1:5], cex = 0.7, pos = 2)

# 2. P-value histogram
hist(results_progression$P.Value, breaks = 50,
     main = "P-value Distribution: Progression Effect",
     xlab = "P-value", col = "lightblue")

# 3. MA plot
plotMA(fit, coef = "ProgressCKM", main = "MA Plot: CKM Progression")

# 4. QQ plot
qqnorm(fit$t[, "ProgressCKM"], main = "QQ Plot: Test Statistics")
qqline(fit$t[, "ProgressCKM"])

par(mfrow = c(1, 1))

# Example: Examine specific protein of interest
protein_of_interest <- "seq.10000.28"
if(protein_of_interest %in% rownames(results_progression)) {
  cat("\n========== DETAILED RESULTS FOR", protein_of_interest, "==========\n")
  
  cat("\nProgression effect:\n")
  print(results_progression[protein_of_interest, ])
  
  cat("\nVisit effect:\n")
  print(results_visit[protein_of_interest, ])
  
  if("ProgressCKM_Visit2" %in% colnames(design)) {
    cat("\nInteraction effect:\n")
    print(results_interaction[protein_of_interest, ])
  }
  
  # Plot individual protein
  protein_values <- expression_matrix[protein_of_interest, ]
  
  par(mfrow = c(1, 2))
  
  # Boxplot by group and visit
  boxplot(protein_values ~ metadata$progress_CKM:metadata$visit_factor,
          xlab = "Progress_CKM:Visit", ylab = "Expression",
          main = paste(protein_of_interest, "by Group and Visit"),
          col = c("lightblue", "lightcoral"))
  
  # Individual trajectories for this protein
  plot_data <- data.frame(
    value = protein_values,
    visit = metadata$visit,
    progress = metadata$progress_CKM,
    subject = metadata$RELEASEID
  )
  
  plot(NA, xlim = c(0.8, 2.2), ylim = range(protein_values),
       xlab = "Visit", ylab = "Expression",
       main = paste(protein_of_interest, "Trajectories"))
  
  for(subj in unique(plot_data$subject)) {
    subj_data <- plot_data[plot_data$subject == subj, ]
    subj_data <- subj_data[order(subj_data$visit), ]
    col_use <- ifelse(subj_data$progress[1] == 1, "red", "blue")
    lines(subj_data$visit, subj_data$value, col = col_use, type = "b")
  }
  legend("topright", c("No Progression", "Progression"), 
         col = c("blue", "red"), lty = 1)
  
  par(mfrow = c(1, 1))
}

cat("\nAnalysis complete!\n")
cat("Results files saved:\n")
cat("  - limma_results_progression.csv\n")
cat("  - limma_results_visit.csv\n")
cat("  - limma_results_interaction.csv\n")