library(redcapAPI)
library(tidyverse)
library(Hmisc)
home_dir <- switch(Sys.info()["sysname"],
  "Darwin" = "/Users/tim/Library/CloudStorage/OneDrive-TheUniversityofColoradoDenver/Vigers/BDC/Janet Snell-Bergeon/CALICO",
  "Linux" = "/home/tim/OneDrive/Vigers/BDC/Janet Snell-Bergeon/CALICO",
  "Windows" = "C:/Users/Tim/OneDrive - The University of Colorado Denver/Vigers/BDC/Janet Snell-Bergeon/CALICO"
)
setwd(home_dir)
# Open REDCap
unlockREDCap(c(rcon = "CALICO"),
  keyring = "API_KEYs",
  envir = 1,
  url = "https://redcap.ucdenver.edu/api/"
)
# Function for converting "UNK", etc. to NA in numeric fields
isMissingSpecial <- function(x, ...) {
  is.na(x) | x == "" | x == "UNK" | x == "OTH" | x == "NA" | x == "PM"
}
# Import history and clinic visits separately for cleanliness
history <- exportReportsTyped(rcon,
  report_id = 131532, warn_zero_coded = F,
  na = list(number = isMissingSpecial, radio = isMissingSpecial)
)
clinic <- exportReportsTyped(rcon,
  report_id = 131533, warn_zero_coded = F,
  na = list(number = isMissingSpecial, radio = isMissingSpecial)
)
# Diagnostic visits are 0 months since diagnosis
clinic$cv_monthssincepcosdx[clinic$cv_visittype == "Diagnostic Visit"] <- 0
# Completed  entries only, drop some columns
history <- history %>%
  filter(history_complete == "Complete") %>%
  select(-redcap_repeat_instrument, -redcap_repeat_instance, -data_entryname)
clinic <- clinic %>%
  filter(clinical_visit_complete == "Complete") %>%
  select(-redcap_repeat_instrument)
# Merge
df <- full_join(history, clinic, by = join_by(record_number)) %>%
  select(record_number, redcap_repeat_instance, everything())
# Weight categories
df <- df %>%
  mutate(
    weight_perc = case_when(
      !is.na(cv_bmi_percentile) & cv_bmi_percentile < 85 ~ "Normal weight",
      !is.na(cv_bmi_percentile) & cv_bmi_percentile >= 85 &
        cv_bmi_percentile < 95 ~ "Overweight",
      !is.na(cv_bmi_percentile) & cv_bmi_percentile >= 95 ~ "Obese",
      !is.na(cv_bmi) & cv_bmi < 25 ~ "Normal weight",
      !is.na(cv_bmi) & cv_bmi >= 25 & cv_bmi < 30 ~ "Overweight",
      !is.na(cv_bmi) & cv_bmi >= 30 ~ "Obese",
      .default = NA
    ),
    weight_raw = cut(cv_bmi,
      breaks = c(-Inf, 25, 30, Inf), right = F,
      labels = c("Normal weight", "Overweight", "Obese")
    )
  )
df$weight_perc <- factor(df$weight_perc,
  levels = c("Normal weight", "Overweight", "Obese")
)
df$weight_raw <- factor(df$weight_raw,
  levels = c("Normal weight", "Overweight", "Obese")
)
# Combined weight category
df$weight_cat <- coalesce(df$weight_perc, df$weight_raw)
# Overweight yes no categories
df <- df %>%
  mutate(
    overweight_perc = case_when(
      !is.na(cv_bmi_percentile) & cv_bmi_percentile >= 85 ~ "Yes",
      !is.na(cv_bmi_percentile) & cv_bmi_percentile < 85 ~ "No",
      is.na(cv_bmi_percentile) & cv_bmi >= 25 ~ "Yes",
      is.na(cv_bmi_percentile) & cv_bmi < 25 ~ "No",
      .default = NA
    ),
    overweight_raw = ifelse(cv_bmi >= 25, "Yes", "No")
  )
df$overweight_perc <- factor(df$overweight_perc, levels = c("No", "Yes"))
df$overweight_raw <- factor(df$overweight_raw, levels = c("No", "Yes"))
# Combined race column
races <- c(
  "Caucasian", "African American", "Asian", "Pacific Islander",
  "American Indian or Alaska Native", "Other", "Unknown/Not recorded",
  "Not applicable", "Other", "Premenarchal"
)
df$combined_race <- apply(
  df[, paste0("race___", c(1:5, 60, "unk"))], 1,
  function(r) {
    w <- which(r == "Checked")
    if (length(w) == 0) {
      return(NA)
    } else if (length(w) > 1) {
      return("More than one")
    } else {
      return(races[w])
    }
  }
)
df$combined_race <- factor(df$combined_race)
# Combined family history - luckily the column numbers for family and non-family
# match up
hist_labels <- sub(
  "Non-parent family history \\(include all history, not only from the diagnostic visit\\)\\. Including siblings, first cousins, aunts, uncles, and grandparents. \\(choice\\=",
  "", as.character(label(df[, grep("pcosdx_famhx___", colnames(df))]))
)
hist_labels <- sub("\\)", "", hist_labels)
hist_labels <- sub(" \\{pcosdx_famhx_family_otherdx\\}", "", hist_labels)
hist_labels <- paste0("Any family history of ", hist_labels, "?")
hist_vars <- sub("pcosdx_famhx___", "", colnames(df[, grep(
  "pcosdx_famhx___", colnames(df)
)]))
any_hist <- data.frame(lapply(hist_vars, function(v) {
  parent <- df[, paste0("pcosdx_famhx_parent___", v)]
  nonparent <- df[, paste0("pcosdx_famhx___", v)]
  any <- (parent == "Checked") + (nonparent == "Checked")
  any <- factor(any, levels = 0:2, labels = c("No", "Yes", "Yes"))
  return(any)
}))
colnames(any_hist) <- paste0("pcosdx_any_famhx___", hist_vars)
df <- cbind(df, any_hist)
# Mental health screening
df$mental_health_screening <-
  rowSums(!is.na(df[, c("cv_phq2", "cv_phq8", "cv_phq9", "cv_cesd20")]))
df$mental_health_screening <- factor(df$mental_health_screening,
  levels = 0:4, labels = c("No", "Yes", "Yes", "Yes", "Yes")
)
# Mental health diagnosis variables
# Convert to numeric
df$cv_phq9[df$cv_phq9 %in%
  c("performed, score not included", "performed, results not listed")] <- NA
df$cv_phq9[df$cv_phq9 == "negative"] <- 0
df$cv_cesd20[df$cv_cesd20 == "AN"] <- NA
df$cv_phq9 <- as.numeric(df$cv_phq9)
df$cv_cesd20 <- as.numeric(df$cv_cesd20)
df <- df %>%
  mutate(
    depression = cv_newdx___16 == "Checked" | pcosdx_pmh___16 == "Checked" |
      cv_phq2 >= 3 | cv_phq8 >= 10 | cv_phq9 >= 10 | cv_cesd20 >= 16,
    anxiety =
      factor(cv_newdx___15 == "Checked" | pcosdx_pmh___15 == "Checked",
        levels = c(F, T), labels = c("No", "Yes")
      ),
    bed =
      factor(cv_newdx___17 == "Checked" | pcosdx_pmh___17 == "Checked",
        levels = c(F, T), labels = c("No", "Yes")
      ),
    red =
      factor(cv_newdx___18 == "Checked" | pcosdx_pmh___18 == "Checked",
        levels = c(F, T), labels = c("No", "Yes")
      ),
    adhd =
      factor(cv_newdx___19 == "Checked" | pcosdx_pmh___19 == "Checked",
        levels = c(F, T), labels = c("No", "Yes")
      )
  )
df$depression[is.na(df$depression)] <- F
df$depression <-
  factor(df$depression, levels = c(F, T), labels = c("No", "Yes"))
# Any mental health counseling
df$pcosdx_any_mentalhealthcounseling <-
  rowSums(df[, c(paste0("pcosdx_mentalhealthcounseling___", 1:3))] == "Checked")
df$pcosdx_any_mentalhealthcounseling <-
  factor(df$pcosdx_any_mentalhealthcounseling,
    levels = 0:3,
    labels = c("No", "Yes", "Yes", "Yes")
  )
# Calculate age for visits where that's not entered
df$cv_age <- (df$cv_monthssincepcosdx / 12) + df$pcosdx_age
# Labels
label(df$depression) <- "Depression?"
label(df$anxiety) <- "Anxiety?"
label(df$bed) <- "Binge Eating Disorder?"
label(df$red) <- "Restrictive Eating Disorder?"
label(df$adhd) <- "ADHD?"
label(df$pcosdx_any_mentalhealthcounseling) <-
  "Any mental health counseling in the 12 months prior to diagnosis?"
# LARC
df$larc <- df$cv_medications___10 == "Checked" |
  df$cv_medications___11 == "Checked"
df$larc <- factor(df$larc, levels = c(F, T), labels = c("No", "Yes"))
# Same as LARC, but for estrogen-containing medications (EC)
df$ec <- df$cv_medications___5 == "Checked" |
  df$cv_medications___6 == "Checked" | df$cv_medications___7 == "Checked"
df$ec <- factor(df$ec, levels = c(F, T), labels = c("No", "Yes"))
# Metformin
df$metformin <- df$cv_medications___1 == "Checked"
df$metformin <- factor(df$metformin, levels = c(F, T), labels = c("No", "Yes"))
# Lifestyle medicine means not on Metformin or EC
df$lifestyle <- df$metformin == "No" & df$ec == "No"
df$lifestyle <- factor(df$lifestyle, levels = c(F, T), labels = c("No", "Yes"))
# Weight loss medication
df$weight_loss <-
  df$cv_medications___25 == "Checked" | df$cv_medications___26 == "Checked" |
    df$cv_medications___27 == "Checked" | df$cv_medications___28 == "Checked" |
    df$cv_medications___29 == "Checked"
df$weight_loss <-
  factor(df$weight_loss, levels = c(F, T), labels = c("No", "Yes"))
# Progesterone containing medication
df$progesterone <-
  df$cv_medications___9 == "Checked" | df$cv_medications___10 == "Checked" |
    df$cv_medications___11 == "Checked" | df$cv_medications___12 == "Checked"
df$progesterone <-
  factor(df$progesterone, levels = c(F, T), labels = c("No", "Yes"))
# Age group at diagnosis
df$age_group <- cut(df$pcosdx_age, c(-Inf, 15, Inf),
  right = F, labels = c("< 15 years", ">= 15 years")
)
# Convert columns to numeric
df$cv_a1c <- suppressWarnings(as.numeric(df$cv_a1c))
df$redcap_repeat_instance <- as.numeric(df$redcap_repeat_instance)
df$pcosdx_age <- as.numeric(df$pcosdx_age)
# Fix/add labels
label(df$insur_type) <- "Insurance Type"
# label(df$Region) <- "Region"
label(df$ethnicity) <- "Ethnicity"
label(df$pcosdx_age) <- "Age at time of PCOS diagnosis"
label(df$age_group) <- "Age Group at Diagnosis"
label(df$combined_race) <- "Race"
label(df$weight_perc) <- "Weight Category (by percentile and raw value)"
label(df$weight_raw) <- "Weight Category (by raw value)"
label(df$weight_cat) <- "Weight Category (by percentile then raw value)"
label(df$overweight_perc) <- "Overweight Status (by percentile and raw value)"
label(df$overweight_raw) <- "Overweight Status (by raw value)"
label(df$cv_a1c) <- "HbA1C"
label(df$larc) <- "On LARC?"
label(df$ec) <- "On Estrogen-Containing Medication?"
label(df$metformin) <- "On Metformin?"
label(df$lifestyle) <- "On Lifestyle Medicine?"
label(df$weight_loss) <- "On Weight Loss Medication?"
label(df$progesterone) <- "On Progesterone-Containing Medication?"
label(df$pcosdx_age) <- "Age at time of PCOS diagnosis"
label(df$mental_health_screening) <-
  "PHQ-2, PHQ-8, PHQ-9 or CED-S Score Available?"
label(df[, grep("___unk", colnames(df))]) <-
  as.list(sub(
    "choice=NA", "choice=Unknown/Not recorded",
    label(df[, grep("___unk", colnames(df))])
  ))
label(df$cv_age) <- "Age at Clinic Visit"
label(df$cv_bmi) <- "BMI (kg/m2)"
# Drop unused levels
l <- label(df)
df <- droplevels(df)
label(df) <- as.list(l)
# Save
save(df, file = "./Data_Clean/analysis_data.RData")
