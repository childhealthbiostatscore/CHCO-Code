## KDIGO CKD staging in long-format data
## One row per person-day; columns assumed: 
##   id = "releaseid", time = "days", eGFR = "ckd_gfr",
##   UACR either in mg/g ("UAlbCreat_mg_g") or as a ratio needing ×1000 ("UAlbCreat")
## Requires: data.table (fast non‑equi joins)

library(data.table)

kdigo_classify <- function(
    df,
    id_col   = "releaseid",
    day_col  = "days",
    egfr_col = "ckd_gfr",
    uacr_mgg_col = "UAlbCreat_mg_g",   # set to NULL if not present
    uacr_raw_col = "UAlbCreat",        # set to NULL if not present
    window_days  = 365L
){
  stopifnot(requireNamespace("data.table", quietly = TRUE))
  DT <- as.data.table(df)
  
  # Harmonize column names into a working table
  DT <- DT[, .(
    id   = get(id_col),
    day  = as.integer(get(day_col)),
    egfr = suppressWarnings(as.numeric(get(egfr_col))),
    uacr_mgg = if (!is.null(uacr_mgg_col) && uacr_mgg_col %in% names(DT))
      suppressWarnings(as.numeric(get(uacr_mgg_col))) else NA_real_,
    uacr_raw = if (!is.null(uacr_raw_col) && uacr_raw_col %in% names(DT))
      suppressWarnings(as.numeric(get(uacr_raw_col))) else NA_real_
  )]
  
  # Create UACR in mg/g (use UAlbCreat_mg_g if present, else convert a small ratio by ×1000)
  DT[, uacr := fifelse(!is.na(uacr_mgg), uacr_mgg,
                       fifelse(!is.na(uacr_raw),
                               # Heuristic: values < 5 are almost surely a ratio (e.g., 0.030 -> 30 mg/g)
                               fifelse(uacr_raw < 5, uacr_raw * 1000, uacr_raw),
                               NA_real_))]
  
  # GFR categories (KDIGO G1–G5)  :contentReference[oaicite:1]{index=1}
  DT[, G_cat := fcase(
    !is.na(egfr) & egfr >= 90,                    "G1",
    !is.na(egfr) & egfr >= 60 & egfr < 90,        "G2",
    !is.na(egfr) & egfr >= 45 & egfr < 60,        "G3a",
    !is.na(egfr) & egfr >= 30 & egfr < 45,        "G3b",
    !is.na(egfr) & egfr >= 15 & egfr < 30,        "G4",
    !is.na(egfr) & egfr < 15,                     "G5",
    default = NA_character_)]
  
  # Albuminuria categories (ACR mg/g): A1 <30; A2 30–300; A3 ≥300  :contentReference[oaicite:2]{index=2}
  DT[, A_cat := fcase(
    !is.na(uacr) & uacr < 30,                     "A1",
    !is.na(uacr) & uacr >= 30 & uacr < 300,       "A2",
    !is.na(uacr) & uacr >= 300,                   "A3",
    default = NA_character_)]
  
  # Lookup table for KDIGO heat‑map risk (low / moderate / high / very_high)  :contentReference[oaicite:3]{index=3}
  risk_lut <- data.table(
    G_cat = rep(c("G1","G2","G3a","G3b","G4","G5"), each = 3),
    A_cat = rep(c("A1","A2","A3"), times = 6),
    risk4 = c(
      # G1
      "low", "moderate", "high",
      # G2
      "low", "moderate", "high",
      # G3a
      "moderate", "high", "very_high",
      # G3b
      "high", "very_high", "very_high",
      # G4
      "very_high", "very_high", "very_high",
      # G5
      "very_high", "very_high", "very_high"
    )
  )
  risk_lut[, risk_rank4 := match(risk4, c("low","moderate","high","very_high"))]
  risk_lut[, risk3 := fifelse(risk4 == "very_high", "very_high",
                              fifelse(risk4 %chin% c("moderate","high"), "mod_high", "low"))]
  risk_lut[, tier3 := match(risk3, c("low","mod_high","very_high")) - 1L]  # 0,1,2
  
  setkey(risk_lut, G_cat, A_cat)
  
  ## ---------- Same‑day combined risk (if both available on the same day) ----------
  same_day <- DT[!is.na(G_cat) & !is.na(A_cat), .(id, day, G_cat, A_cat)]
  same_day <- risk_lut[same_day, on = .(G_cat, A_cat)]
  # If multiple results on a day, take the worst (max rank)
  same_day <- same_day[, .(
    risk_rank4_sd = max(risk_rank4, na.rm = TRUE)
  ), by = .(id, day)]
  same_day[, risk4_sd := c("low","moderate","high","very_high")[risk_rank4_sd]]
  same_day[, risk3_sd := fifelse(risk4_sd == "very_high", "very_high",
                                 fifelse(risk4_sd %chin% c("moderate","high"), "mod_high", "low"))]
  
  ## ---------- Cross‑day pairing within +/- window_days; assign to the earlier day ----------
  eg <- DT[!is.na(G_cat), .(id, G_cat, g_day = day, start = day, end = day)]
  ua <- DT[!is.na(A_cat), .(id, A_cat, a_day = day, start = day - window_days, end = day + window_days)]
  setkey(eg, id, start, end)
  setkey(ua, id, start, end)
  
  pairs <- foverlaps(eg, ua, nomatch = 0L)  # eg overlaps ua windows
  # Compute risk for each pair and assign to the earlier of the two measurements
  pairs <- risk_lut[pairs, on = .(G_cat, A_cat)]
  pairs[, earlier_day := pmin(g_day, a_day)]
  paired_earlier <- pairs[, .(
    risk_rank4_pair = max(risk_rank4, na.rm = TRUE)  # worst risk if multiple pairs map to same earlier_day
  ), by = .(id, day = earlier_day)]
  paired_earlier[, risk4_pair := c("low","moderate","high","very_high")[risk_rank4_pair]]
  paired_earlier[, risk3_pair := fifelse(risk4_pair == "very_high", "very_high",
                                         fifelse(risk4_pair %chin% c("moderate","high"), "mod_high", "low"))]
  
  ## ---------- Merge back to all person-days and select the worst risk on each day ----------
  setkey(DT, id, day)
  OUT <- copy(DT)
  OUT <- same_day[OUT, on = .(id, day)]
  setnames(OUT, old = c("risk_rank4_sd","risk4_sd","risk3_sd"),
           new = c("risk_rank4_sd","risk4_sd","risk3_sd"))
  OUT <- paired_earlier[OUT, on = .(id, day)]
  setnames(OUT, old = c("risk_rank4_pair","risk4_pair","risk3_pair"),
           new = c("risk_rank4_pair","risk4_pair","risk3_pair"))
  
  # Final combined risk for the day = worst of same-day and paired-within-365 (if any)
  OUT[, risk_rank4 := fifelse(!is.na(risk_rank4_sd) & !is.na(risk_rank4_pair),
                              pmax(risk_rank4_sd, risk_rank4_pair),
                              fifelse(!is.na(risk_rank4_sd), risk_rank4_sd,
                                      fifelse(!is.na(risk_rank4_pair), risk_rank4_pair, NA_integer_)))]
  OUT[, kdigo_risk4 := ifelse(is.na(risk_rank4), NA_character_,
                              c("low","moderate","high","very_high")[risk_rank4])]
  OUT[, kdigo_risk3 := fifelse(is.na(kdigo_risk4), NA_character_,
                               fifelse(kdigo_risk4 == "very_high", "very_high",
                                       fifelse(kdigo_risk4 %chin% c("moderate","high"), "mod_high", "low")))]
  OUT[, kdigo_tier := match(kdigo_risk3, c("low","mod_high","very_high")) - 1L]  # 0/1/2 with NA for missing
  
  # Source flag (for auditing)
  OUT[, kdigo_source := fifelse(!is.na(risk_rank4_sd) & !is.na(risk_rank4_pair), "same_day+paired",
                                fifelse(!is.na(risk_rank4_sd), "same_day",
                                        fifelse(!is.na(risk_rank4_pair), "paired_365", NA_character_)))]
  
  ## ---------- Progression days (3-level scheme: low | mod_high | very_high) ----------
  # "Low risk is baseline risk": treat no prior observation as prior=low (tier 0)
  # Progression events: low->mod/high ; low->very_high ; mod/high->very_high
  setorder(OUT, id, day)
  OUT[, last_obs_tier := shift(nafill(kdigo_tier, type = "locf"), type = "lag"), by = id]
  OUT[, last_obs_tier := fifelse(is.na(last_obs_tier), 0L, last_obs_tier)]  # baseline low if none prior
  OUT[, progression_type := fcase(
    !is.na(kdigo_tier) & kdigo_tier == 1L & last_obs_tier == 0L, "low_to_mod_high",
    !is.na(kdigo_tier) & kdigo_tier == 2L & last_obs_tier == 0L, "low_to_very_high",
    !is.na(kdigo_tier) & kdigo_tier == 2L & last_obs_tier == 1L, "mod_high_to_very_high",
    default = NA_character_)]
  OUT[, progression_event := !is.na(progression_type)]
  OUT[, progression_event := as.integer(progression_event)]
  
  # Per-person count of days with a progression event (constant within id)
  prog_counts <- OUT[progression_event == 1L, .N, by = id]
  setnames(prog_counts, "N", "n_progression_days")
  OUT <- prog_counts[OUT, on = "id"]
  OUT[is.na(n_progression_days), n_progression_days := 0L]
  
  # Return augmented long table
  OUT[]
}

## ------------------ Example usage ------------------
## df <- read.csv("your_long_table.csv")
## res <- kdigo_classify(df,
##                       id_col="releaseid", day_col="days",
##                       egfr_col="ckd_gfr",
##                       uacr_mgg_col="UAlbCreat_mg_g",
##                       uacr_raw_col="UAlbCreat")
## head(res[order(id, day)])
## res[progression_event == 1L, .(id, day, progression_type, kdigo_risk3)][order(id, day)]
