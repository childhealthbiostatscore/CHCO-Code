# Install and load required packages
install.packages("lme4")
library(lme4)

# Load your scRNA-seq data (replace 'your_expression_matrix' with your actual data)
# Assuming your data has columns 'Gene', 'Cell', 'Condition', and 'Expression'
# For example:
# sc_data <- data.frame(Gene = rownames(your_expression_matrix),
#                       Cell = colnames(your_expression_matrix),
#                       Condition = c("Condition1", "Condition1", "Condition2", "Condition2", ...),
#                       Expression = as.vector(your_expression_matrix))

# Filter data to include only Condition 1 and Condition 2
condition_data <- subset(sc_data, Condition %in% c("Condition1", "Condition2"))

# Fit the mixed-effects model comparing Condition 1 vs Condition 2
model <- lmer(Expression ~ Gene + Condition + (1|Patient), data = condition_data)

# Summarize the model results
summary(model)

# Assuming 'model' is the fitted mixed-effects model

# Extract the coefficient estimates and their associated p-values
coefficients_table <- summary(model)$coefficients

# Filter significant genes based on the p-value threshold (e.g., p < 0.05)
significant_genes <- coefficients_table[coefficients_table[, "Pr(>|t|)"] < 0.05, ]

# Display the significant genes
print(significant_genes)



# Assuming 'model' is the fitted mixed-effects model

# Extract the coefficient estimates and their associated p-values
coefficients_table <- summary(model)$coefficients

# Extract the p-values from the table
p_values <- coefficients_table[, "Pr(>|t|)"]

# Perform FDR correction
adjusted_p_values <- p.adjust(p_values, method = "fdr")

# Add the adjusted p-values to the coefficients table
coefficients_table$Adj_P_Value <- adjusted_p_values

# Filter significant genes based on the adjusted p-value threshold (e.g., adjusted p < 0.05)
significant_genes <- coefficients_table[coefficients_table$Adj_P_Value < 0.05, ]

# Display the significant genes
print(significant_genes)


