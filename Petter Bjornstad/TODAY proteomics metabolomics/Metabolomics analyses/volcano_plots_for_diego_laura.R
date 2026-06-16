
### more visually appealing volcano plots for Diego/Laura posters ###
full_results <- datasets[["HTN"]]

# ---- P-value columns ----
full_results$p_adjusted   <- full_results$adj.p.value
full_results$p10_adjusted   <- -log10(pmax(full_results$p_adjusted,   1e-10))

# ---- Color assignment ----
full_results$color_adjusted <- case_when(
  full_results$p_adjusted < 0.05 & full_results$estimate > 1 ~ "#FF2020",
  full_results$p_adjusted < 0.05 & full_results$estimate < 1 ~ "#1E90FF",
  full_results$metabolite_type == "urine"                    ~ "#DAA520",
  full_results$metabolite_type == "plasma"                   ~ "#9B59B6",
  TRUE                                                       ~ "gray70"
)

# ---- Significant subsets ----
significant_df_adjusted   <- full_results[full_results$p_adjusted   < 0.05, ]

# ---- X-axis limits ----
x_min <- min(full_results$estimate)
x_max <- max(full_results$estimate)

p_adjusted <- ggplot(full_results, aes(x = estimate, y = p10_adjusted, color = color_adjusted)) +
  geom_point(alpha = 0.7) +
  scale_color_identity() +
  theme_minimal() +
  labs(
    x     = "HR",
    y     = "-log10(P-Value)",
    color = "LogFC Direction") +
  theme(
    panel.grid.major  = element_blank(),
    panel.grid.minor  = element_blank(),
    axis.text         = element_text(size = 13),
    axis.title        = element_text(size = 15),
    strip.text        = element_blank()) +
  scale_x_continuous(limits = c(x_min, x_max), expand = expansion(mult = 0.2)) +
  coord_cartesian(clip = "off") +
  facet_wrap(~ metabolite_type, ncol = 2, scales="free") +
  geom_text_repel(
    data         = significant_df_adjusted,
    aes(label    = term),
    size         = 4,
    color        = "black",
    max.overlaps = Inf
  )
p_adjusted
