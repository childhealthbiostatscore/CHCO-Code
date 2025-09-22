############################################################
#                                                          #
#   Author: Long Yuan                                      #
#   Email : lyuan13@jhmi.edu                               #
#                                                          #
#           __                                             #
#       (___()'`;                                          #
#       /,    /`                                           #
#       \\\"--\\                                            #
#                                                          #
#   Notebook Description:                                  #
#   [Pseudopair with Hungarian Matching for Data with      #
#    Missingness] - R Version                              #
############################################################

# Load required libraries
library(dplyr)
library(readr)
library(VIM)  # for mice imputation
library(mice)
library(prcomp)
library(ggplot2)
library(ggrepel)
library(RColorBrewer)
library(clue)  # for Hungarian algorithm (solve_LSAP)

# Read and prepare data
df_full <- read_csv('pb90_rhrh2improve_clinical_subset.csv')

desired_cols <- c('record_id', 'date', 'improve_id', 'rh_id', 'rh2_id',
                  'sex', 'study', 'group', 'hba1c', 'eGFR_CKD_epi', 'age', 'bmi')

df_subset <- df_full %>%
  select(all_of(desired_cols)) %>%
  mutate(date = as.Date(date))

# --- RH-* co-enrolled with IMPROVE ---
rh_ids <- c('RH-59-T', 'RH-60-T', 'RH-66-T')
improve_ids <- c('IT_07', 'IT_08', 'IT_10')
co_enrolled_1 <- df_subset %>% 
  filter(record_id %in% c(rh_ids, improve_ids))

# --- RH2-* co-enrolled in RH-* ---
rh2_ids <- c('RH2-14-T', 'RH2-19-T')
rh_ids_from_rh2 <- c('RH-23-T', 'RH-67-T')
co_enrolled_2 <- df_subset %>% 
  filter(record_id %in% c(rh2_ids, rh_ids_from_rh2))

# --- IMPROVE patients with complete pairs ---
multi_visit_ids <- df_subset %>%
  count(record_id) %>%
  filter(n > 1) %>%
  pull(record_id)

co_enrolled_3 <- df_subset %>% 
  filter(record_id %in% c('IT_11', 'IT_12'))

# Combine all
complete_pairs_subset <- bind_rows(co_enrolled_1, co_enrolled_2, co_enrolled_3) %>%
  arrange(record_id, date)

print(complete_pairs_subset)

# Encode sex and prepare for imputation
df_subset <- df_subset %>%
  mutate(sex_encoded = ifelse(sex == "Male", 1, 0))

predictors <- df_subset %>%
  select(hba1c, eGFR_CKD_epi, age, bmi, sex_encoded)

# Perform multiple imputation using mice
mice_result <- mice(predictors, m = 5, method = 'pmm', seed = 0, printFlag = FALSE)
imputed_data <- complete(mice_result)

# Update df_subset with imputed values
df_subset$hba1c <- imputed_data$hba1c
df_subset$eGFR_CKD_epi <- imputed_data$eGFR_CKD_epi

# Prepare features for PCA
X <- df_subset %>%
  select(hba1c, eGFR_CKD_epi, age, bmi, sex_encoded) %>%
  scale()

# Perform PCA
pca_result <- prcomp(X, center = FALSE, scale. = FALSE)
df_subset$pca1 <- pca_result$x[, 1]
df_subset$pca2 <- pca_result$x[, 2]

# Create PCA plot with labels
p1 <- ggplot(df_subset, aes(x = pca1, y = pca2, color = sex, shape = study)) +
  geom_point(size = 3) +
  geom_text_repel(aes(label = record_id), size = 3, fontface = "bold") +
  labs(x = "PC-1", y = "PC-2", title = "PCA Plot with Record IDs") +
  theme_minimal() +
  theme(legend.position = "right")

print(p1)
ggsave('mice_pca.png', p1, width = 10, height = 8, dpi = 600)

# Prepare for Hungarian matching
complete_ids <- c('RH-59-T','RH-60-T','RH-66-T','IT_07','IT_08','IT_10',
                  'RH2-14-T','RH2-19-T','RH-23-T','RH-67-T','IT_11','IT_12')

# Exclude complete pairs
df <- df_subset %>%
  filter(!record_id %in% complete_ids)

# Split by median date
median_date <- median(df$date)
baseline <- df %>% filter(date <= median_date)
followup <- df %>% filter(date > median_date)

cat("Baseline:", nrow(baseline), "Followup:", nrow(followup), "\n")

# Build cost matrix for Hungarian algorithm
n_f <- nrow(followup)
n_b <- nrow(baseline)
C <- matrix(0, nrow = n_f, ncol = n_b)
BIG <- 1e6

# Calculate costs for each pair
for (i in 1:n_f) {
  for (j in 1:n_b) {
    post <- followup[i, ]
    base <- baseline[j, ]
    
    # Sex must match
    sex_pen <- ifelse(post$sex == base$sex, 0, BIG)
    
    # BMI percentile difference
    bmi_pen <- post$bmi - base$bmi
    
    # Timepoint: ideally 12-48 mo apart
    months <- as.numeric(post$date - base$date) / 30.44
    dev <- abs(months - 12)
    time_pen <- ifelse(months >= 12 & months <= 48, dev, dev + BIG/10)
    
    # Age: ideally 0–4 y older
    ad <- post$age - base$age
    age_pen <- ifelse(ad >= 0 & ad <= 4, abs(ad - 1), abs(ad - 1) + BIG/1000)
    
    # HbA1c: post should be lower
    hba_pen <- post$hba1c - base$hba1c
    
    # eGFR: if baseline <90, post must be higher
    egfr_pen <- ifelse(base$eGFR_CKD_epi < 90, 
                       max(0, base$eGFR_CKD_epi - post$eGFR_CKD_epi), 
                       0)
    
    C[i, j] <- sex_pen + bmi_pen + time_pen + age_pen + hba_pen + egfr_pen
  }
}

# Solve assignment using Hungarian algorithm
assignment <- solve_LSAP(C)

# Create pairs dataframe
pairs_df <- data.frame(
  baseline_id = baseline$record_id[assignment],
  followup_id = followup$record_id,
  cost = C[cbind(1:n_f, assignment)]
) %>%
  arrange(cost)

print(pairs_df)

# Filter good pairs (cost < 1000)
good_pairs_df <- pairs_df %>% filter(cost < 1000)

# Define complete pairs
complete_pairs_df <- data.frame(
  record_id_pre = c('RH-59-T', 'RH-60-T', 'RH-66-T', 'RH-23-T', 'RH-67-T', 'IT_11', 'IT_12'),
  date_pre = as.Date(c('2019-07-11', '2019-07-26', '2019-10-03', '2018-11-16', 
                       '2019-12-16', '2019-12-12', '2020-01-09')),
  record_id_post = c('IT_07', 'IT_08', 'IT_10', 'RH2-14-T', 'RH2-19-T', 'IT_11', 'IT_12'),
  date_post = as.Date(c('2021-01-28', '2021-02-09', '2020-11-10', '2022-10-31',
                        '2023-03-15', '2021-07-21', '2021-09-21'))
)

# Create pseudopairs visualization
create_arrow_plot <- function() {
  p <- ggplot(df_subset, aes(x = pca1, y = pca2)) +
    geom_point(color = 'lightgray', alpha = 0.7, size = 2) +
    theme_minimal()
  
  # Add arrows for good pairs
  for (i in 1:nrow(good_pairs_df)) {
    pre_data <- df_subset %>% filter(record_id == good_pairs_df$baseline_id[i])
    post_data <- df_subset %>% filter(record_id == good_pairs_df$followup_id[i])
    
    if (nrow(pre_data) == 1 && nrow(post_data) == 1) {
      p <- p + geom_segment(
        x = pre_data$pca1, y = pre_data$pca2,
        xend = post_data$pca1, yend = post_data$pca2,
        arrow = arrow(length = unit(0.2, "cm")),
        color = 'steelblue', alpha = 0.7, size = 1
      )
    }
  }
  
  # Add arrows for complete pairs
  for (i in 1:nrow(complete_pairs_df)) {
    pre_data <- df_subset %>% 
      filter(record_id == complete_pairs_df$record_id_pre[i], 
             date == complete_pairs_df$date_pre[i])
    post_data <- df_subset %>% 
      filter(record_id == complete_pairs_df$record_id_post[i],
             date == complete_pairs_df$date_post[i])
    
    if (nrow(pre_data) == 1 && nrow(post_data) == 1) {
      p <- p + geom_segment(
        x = pre_data$pca1, y = pre_data$pca2,
        xend = post_data$pca1, yend = post_data$pca2,
        arrow = arrow(length = unit(0.2, "cm")),
        color = 'crimson', alpha = 0.8, size = 1.2
      )
    }
  }
  
  # Get all IDs for annotation
  all_ids <- unique(c(good_pairs_df$baseline_id, good_pairs_df$followup_id,
                      complete_pairs_df$record_id_pre, complete_pairs_df$record_id_post))
  annotated <- df_subset %>% filter(record_id %in% all_ids)
  
  # Add points and labels for annotated data
  p <- p + 
    geom_point(data = annotated, aes(color = sex, shape = study), size = 3) +
    geom_text_repel(data = annotated, aes(label = record_id), 
                    size = 3, fontface = "bold") +
    labs(title = "Complete and Pseudopairs", x = "PC1", y = "PC2") +
    theme(legend.position = "right")
  
  return(p)
}

p2 <- create_arrow_plot()
print(p2)
ggsave('pseudopairs.png', p2, width = 10, height = 8, dpi = 600)

# Calculate distances
features <- c('hba1c', 'eGFR_CKD_epi', 'age', 'bmi', 'sex_encoded')
df_scaled <- df_subset
df_scaled[features] <- scale(df_scaled[features])

# Function to compute pairwise distance
compute_distance <- function(id1, id2, date1 = NULL, date2 = NULL, df) {
  if (is.null(date1)) {
    pre <- df %>% filter(record_id == id1)
    post <- df %>% filter(record_id == id2)
  } else {
    pre <- df %>% filter(record_id == id1, date == date1)
    post <- df %>% filter(record_id == id2, date == date2)
  }
  
  if (nrow(pre) != 1 || nrow(post) != 1) return(NA)
  
  pre_vec <- as.numeric(pre[features])
  post_vec <- as.numeric(post[features])
  return(sqrt(sum((post_vec - pre_vec)^2)))
}

# Calculate distances for good pairs
good_pairs_df$distance <- mapply(compute_distance, 
                                 good_pairs_df$baseline_id, 
                                 good_pairs_df$followup_id,
                                 MoreArgs = list(df = df_scaled))

# Calculate distances for complete pairs
complete_pairs_df$distance <- mapply(compute_distance,
                                     complete_pairs_df$record_id_pre,
                                     complete_pairs_df$record_id_post,
                                     complete_pairs_df$date_pre,
                                     complete_pairs_df$date_post,
                                     MoreArgs = list(df = df_scaled))

print("Good pairs with distances:")
print(good_pairs_df)
print("Complete pairs with distances:")
print(complete_pairs_df)

# Create distance comparison plot
df_dist <- bind_rows(
  complete_pairs_df %>% 
    select(record_id_pre, record_id_post, distance) %>%
    mutate(type = 'complete'),
  good_pairs_df %>% 
    select(baseline_id, followup_id, distance) %>%
    rename(record_id_pre = baseline_id, record_id_post = followup_id) %>%
    mutate(type = 'pseudo')
) %>%
  filter(!is.na(distance))

p3 <- ggplot(df_dist, aes(x = type, y = distance, fill = type)) +
  geom_violin(alpha = 0.7) +
  geom_boxplot(width = 0.1, fill = "white", alpha = 0.8) +
  scale_fill_brewer(type = "qual", palette = "Set2") +
  labs(title = "Distribution of Pairwise Euclidean Distances",
       x = "Pair Type", y = "Distance") +
  theme_minimal() +
  theme(legend.position = "none")

print(p3)
ggsave('distance_comparison.png', p3, width = 6, height = 8, dpi = 600)

# Statistical tests
dist_complete <- df_dist %>% filter(type == 'complete') %>% pull(distance)
dist_pseudo <- df_dist %>% filter(type == 'pseudo') %>% pull(distance)

# Mann-Whitney U test
wilcox_result <- wilcox.test(dist_complete, dist_pseudo, alternative = "two.sided")

# Welch's t-test
t_result <- t.test(dist_complete, dist_pseudo, var.equal = FALSE)

stat_results <- data.frame(
  Test = c("Mann-Whitney U", "Welch's t-test"),
  Statistic = c(wilcox_result$statistic, t_result$statistic),
  p_value = c(wilcox_result$p.value, t_result$p.value)
)

print("Statistical test results:")
print(stat_results)

# Compute Standardized Mean Differences (SMD) for Love Plot
compute_smd <- function(df, complete_pairs, pseudo_pairs, features) {
  smd_results <- data.frame(
    feature = features,
    smd_complete = NA,
    smd_pseudo = NA
  )
  
  for (i in 1:length(features)) {
    feat <- features[i]
    
    # Complete pairs differences
    comp_diffs <- c()
    for (j in 1:nrow(complete_pairs)) {
      pre <- df %>% filter(record_id == complete_pairs$record_id_pre[j])
      post <- df %>% filter(record_id == complete_pairs$record_id_post[j])
      if (nrow(pre) == 1 && nrow(post) == 1) {
        comp_diffs <- c(comp_diffs, post[[feat]] - pre[[feat]])
      }
    }
    
    # Pseudo pairs differences
    pseudo_diffs <- c()
    for (j in 1:nrow(pseudo_pairs)) {
      pre <- df %>% filter(record_id == pseudo_pairs$baseline_id[j])
      post <- df %>% filter(record_id == pseudo_pairs$followup_id[j])
      if (nrow(pre) == 1 && nrow(post) == 1) {
        pseudo_diffs <- c(pseudo_diffs, post[[feat]] - pre[[feat]])
      }
    }
    
    # Calculate SMD
    pooled_std <- sd(c(comp_diffs, pseudo_diffs), na.rm = TRUE)
    if (pooled_std > 0) {
      smd_results$smd_complete[i] <- mean(comp_diffs, na.rm = TRUE) / pooled_std
      smd_results$smd_pseudo[i] <- mean(pseudo_diffs, na.rm = TRUE) / pooled_std
    }
  }
  
  return(smd_results)
}

smd_df <- compute_smd(df_subset, complete_pairs_df, good_pairs_df, features)

# Create Love Plot
smd_long <- smd_df %>%
  tidyr::pivot_longer(cols = c(smd_complete, smd_pseudo), 
                      names_to = "type", values_to = "smd") %>%
  mutate(type = gsub("smd_", "", type))

p4 <- ggplot(smd_long, aes(x = feature, y = smd, color = type, group = type)) +
  geom_line(size = 1) +
  geom_point(size = 3) +
  geom_hline(yintercept = c(-0.1, 0.1), linetype = "dashed", color = "gray") +
  geom_hline(yintercept = 0, linetype = "solid", color = "black") +
  labs(title = "Love Plot: Standardized Mean Differences",
       x = "Feature", y = "Standardized Mean Difference",
       color = "Pair Type") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(p4)
ggsave('love_plot.png', p4, width = 8, height = 6, dpi = 600)

cat("Analysis complete! All plots saved.\n")