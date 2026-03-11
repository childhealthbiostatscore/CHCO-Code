## QC script for attempt_clean.R edits
## Variables checked: age, height, sex, creatinine_s, cystatin_c_s

vars <- c("age", "height", "sex", "creatinine_s", "cystatin_c_s")

cat("=== EXPECTED RANGES ===\n")
cat("age:          0–100 years\n")
cat("height:       100–250 cm\n")
cat("creatinine_s: 0.1–15 mg/dL\n")
cat("cystatin_c_s: 0.1–10 mg/L\n")
cat("sex:          Male / Female only\n\n")

# ── 1. Duplicate record_id × visit ──────────────────────────────────────────
cat("=== DUPLICATES (record_id × visit) ===\n")
dups <- merged_data %>%
  count(record_id, visit) %>%
  filter(n > 1)
if (nrow(dups) == 0) {
  cat("No duplicates found.\n\n")
} else {
  cat(nrow(dups), "duplicate record_id × visit combinations:\n")
  print(dups)
  cat("\n")
}

# ── 2. Missing visits (subjects with < expected number of visits) ────────────
cat("=== VISITS PER SUBJECT ===\n")
visit_counts <- merged_data %>%
  count(record_id, name = "n_visits")
print(table(visit_counts$n_visits))
cat("\n")

subjects_missing_baseline <- merged_data %>%
  group_by(record_id) %>%
  summarise(has_baseline = any(visit == "baseline", na.rm = TRUE)) %>%
  filter(!has_baseline)
cat(nrow(subjects_missing_baseline), "subject(s) missing a baseline visit\n\n")

# ── 3. Missingness summary ───────────────────────────────────────────────────
cat("=== MISSINGNESS SUMMARY ===\n")
miss_summary <- merged_data %>%
  summarise(across(all_of(vars), ~ sum(is.na(.)))) %>%
  pivot_longer(everything(), names_to = "variable", values_to = "n_missing") %>%
  mutate(pct_missing = round(n_missing / nrow(merged_data) * 100, 1))
print(miss_summary)
cat("\n")

# ── 4. Missingness by visit ──────────────────────────────────────────────────
cat("=== MISSINGNESS BY VISIT ===\n")
miss_by_visit <- merged_data %>%
  group_by(visit) %>%
  summarise(across(all_of(vars), ~ sum(is.na(.))), n = n()) %>%
  arrange(visit)
print(miss_by_visit)
cat("\n")

# ── 5. Range checks ─────────────────────────────────────────────────────────
cat("=== RANGE CHECKS ===\n")

numeric_vars <- c("age", "height", "creatinine_s", "cystatin_c_s")
range_summary <- merged_data %>%
  summarise(across(all_of(numeric_vars), list(
    min  = ~ min(.x, na.rm = TRUE),
    max  = ~ max(.x, na.rm = TRUE),
    mean = ~ round(mean(.x, na.rm = TRUE), 2),
    sd   = ~ round(sd(.x, na.rm = TRUE), 2)
  ))) %>%
  pivot_longer(everything(),
               names_to  = c("variable", "stat"),
               names_sep = "_(?=[^_]+$)") %>%
  pivot_wider(names_from = stat, values_from = value)
print(range_summary)
cat("\n")

# Flag out-of-range values
cat("--- Out-of-range flags ---\n")
out_of_range <- merged_data %>%
  filter(
    age < 0 | age > 100 |
    height < 100 | height > 250 |
    creatinine_s < 0.1 | creatinine_s > 15 |
    cystatin_c_s < 0.1 | cystatin_c_s > 10
  ) %>%
  select(record_id, visit, any_of(numeric_vars))
if (nrow(out_of_range) == 0) {
  cat("No out-of-range values found.\n\n")
} else {
  cat(nrow(out_of_range), "row(s) with out-of-range values:\n")
  print(out_of_range)
  cat("\n")
}

# ── 6. Sex values ────────────────────────────────────────────────────────────
cat("=== SEX VALUES ===\n")
print(table(merged_data$sex, useNA = "ifany"))
unexpected_sex <- merged_data %>%
  filter(!sex %in% c("Male", "Female") | is.na(sex))
cat(nrow(unexpected_sex), "row(s) with unexpected/missing sex values\n\n")

cat("=== QC COMPLETE ===\n")
