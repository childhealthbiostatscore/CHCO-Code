# Visualize how individual variation affects disease detection
library(ggplot2)

set.seed(123)
n_cells <- 100
n_individuals <- 4

# Low individual variation
low_var_data <- data.frame(
  expression = c(
    rnorm(n_cells*2, mean = 5, sd = 0.5),   # Controls (similar)
    rnorm(n_cells*2, mean = 10, sd = 0.5)   # Patients (similar)
  ),
  group = rep(c("Control", "Patient"), each = n_cells*2),
  individual = rep(1:4, each = n_cells),
  scenario = "Low Individual Variation (0.1)"
)

# High individual variation  
high_var_data <- data.frame(
  expression = c(
    c(rnorm(n_cells, mean = 3, sd = 0.5),   # Control 1 (low expresser)
      rnorm(n_cells, mean = 7, sd = 0.5)),  # Control 2 (high expresser)
    c(rnorm(n_cells, mean = 8, sd = 0.5),   # Patient 1 (low expresser)
      rnorm(n_cells, mean = 12, sd = 0.5))  # Patient 2 (high expresser)
  ),
  group = rep(c("Control", "Patient"), each = n_cells*2),
  individual = rep(1:4, each = n_cells),
  scenario = "High Individual Variation (0.5)"
)

# Combine and plot
combined_data <- rbind(low_var_data, high_var_data)

ggplot(combined_data, aes(x = group, y = expression, fill = group)) +
  geom_violin(alpha = 0.5) +
  geom_boxplot(width = 0.2, outlier.size = 0.5) +
  facet_wrap(~scenario) +
  theme_minimal() +
  labs(title = "Impact of Individual Variation on Disease Signal",
       subtitle = "Same disease effect, different individual variation")