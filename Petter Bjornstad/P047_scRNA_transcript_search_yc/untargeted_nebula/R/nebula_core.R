# =============================================================================
# nebula_core.R  —  derive metadata, subset/relevel, untargeted NEBULA, output
# Sourced after setup.R + config.R.
# =============================================================================

# KPMP general (low-res) mapping, identical to transcript_search.qmd
kpmp_general_map <- function(ct) {
  dplyr::case_when(
    ct %in% c("PT-S1/S2", "PT-S3", "aPT")                              ~ "PT",
    ct %in% c("C-TAL-1", "C-TAL-2", "aTAL", "dTAL")                    ~ "TAL",
    ct %in% c("CCD-PC", "CNT-PC", "dCCD-PC", "M-PC", "tPC-IC")          ~ "PC",
    ct %in% c("IC-A", "IC-B", "aIC")                                   ~ "IC",
    ct %in% c("EC-AVR", "EC-GC", "EC-PTC", "EC-AEA", "EC-LYM", "EC/VSMC") ~ "EC",
    ct %in% c("MAC", "MON", "cDC", "pDC")                              ~ "Immune (M)",
    ct %in% c("CD4+ T", "CD8+ T", "B", "NK")                           ~ "Immune (L)",
    ct %in% c("VSMC/P", "FIB")                                         ~ "VSMC/P, FIB",
    ct %in% c("POD")                                                   ~ "POD",
    TRUE                                                               ~ "Other"
  )
}

# ---------------------------------------------------------------------------
# derive_meta(): add all columns the contrasts need, on the FULL object,
# BEFORE splitting. Runs once per dataset in the split step.
# ---------------------------------------------------------------------------
derive_meta <- function(obj, dataset, resolutions, nebula_cfg) {
  md <- obj@meta.data

  # ensure both cell-type resolution columns exist
  spec_col <- resolutions$specific
  gen_col  <- resolutions$general
  if (!spec_col %in% names(md))
    stop("Specific cell-type column '", spec_col, "' missing from ", dataset)
  if (!gen_col %in% names(md)) {
    message("  deriving ", gen_col, " from ", spec_col)
    md[[gen_col]] <- kpmp_general_map(md[[spec_col]])
  }

  if (dataset == "pb90") {
    # DKD groups (eGFR<90 OR uACR>=cutoff)
    md <- md %>%
      mutate(
        dkd_group_30 = dplyr::case_when(
          eGFR_CKD_epi <  90 | acr_u >= 30  ~ "DKD",
          eGFR_CKD_epi >= 90 & acr_u <  30  ~ "non_DKD",
          TRUE ~ NA_character_),
        dkd_group_30 = factor(dkd_group_30, levels = c("non_DKD", "DKD")),
        dkd_group_100 = dplyr::case_when(
          eGFR_CKD_epi <  90 | acr_u >= 100 ~ "DKD",
          eGFR_CKD_epi >= 90 & acr_u <  100 ~ "non_DKD",
          TRUE ~ NA_character_),
        dkd_group_100 = factor(dkd_group_100, levels = c("non_DKD", "DKD"))
      )
  }

  # binary pathology indicators (0 = none, >=1 = present), wherever the
  # severity column exists. Derived for any dataset that has them.
  for (sev in c("arteriosclerosis_sev", "arteriolohyalinosis_sev")) {
    if (sev %in% names(md)) {
      yn <- ifelse(md[[sev]] > 0, "Yes", "No")
      md[[sub("_sev$", "_yn", sev)]] <- factor(yn, levels = c("No", "Yes"))
    }
  }

  # obesity classifications, computed wherever the required inputs exist.
  # Pediatric (age <= 19) uses BMI percentile (bmip); adults use BMI.
  if (all(c("age", "bmi", "bmip") %in% names(md))) {
    md$bmi_obesity <- factor(dplyr::case_when(
      md$age >  19 & md$bmi <  25                   ~ "Normal",
      md$age >  19 & md$bmi >= 25 & md$bmi  < 30    ~ "Overweight",
      md$age >  19 & md$bmi >= 30                   ~ "Obese",
      md$age <= 19 & md$bmip <  85                  ~ "Normal",
      md$age <= 19 & md$bmip >= 85 & md$bmip < 95   ~ "Overweight",
      md$age <= 19 & md$bmip >= 95                  ~ "Obese",
      TRUE ~ NA_character_), levels = c("Normal", "Overweight", "Obese"))
  }
  # body-fat % by DXA (sex-specific cutoffs). Requires sex coded "Male"/"Female".
  if (all(c("dexa_body_fat", "sex") %in% names(md))) {
    md$dxa_obesity <- factor(dplyr::case_when(
      md$dexa_body_fat <  25 & md$sex == "Male"                            ~ "Normal",
      md$dexa_body_fat <  32 & md$sex == "Female"                          ~ "Normal",
      md$dexa_body_fat >= 25 & md$dexa_body_fat < 30 & md$sex == "Male"     ~ "Overweight",
      md$dexa_body_fat >= 32 & md$dexa_body_fat < 35 & md$sex == "Female"   ~ "Overweight",
      md$dexa_body_fat >= 30 & md$sex == "Male"                            ~ "Obese",
      md$dexa_body_fat >= 35 & md$sex == "Female"                          ~ "Obese",
      TRUE ~ NA_character_), levels = c("Normal", "Overweight", "Obese"))
  }
  # waist-to-height ratio classifications
  if (all(c("waistcm", "height") %in% names(md))) {
    md$whtr <- md$waistcm / md$height
    md$whtr_obesity <- factor(dplyr::case_when(
      md$whtr <  0.5                    ~ "Normal",
      md$whtr >= 0.5 & md$whtr < 0.6    ~ "Overweight",
      md$whtr >= 0.6                    ~ "Obese",
      TRUE ~ NA_character_), levels = c("Normal", "Overweight", "Obese"))
    md$whtr_obese_binary <- factor(dplyr::case_when(
      md$whtr_obesity == "Obese"                      ~ "Obese",
      md$whtr_obesity %in% c("Normal", "Overweight")  ~ "Non_Obese",
      TRUE ~ NA_character_), levels = c("Non_Obese", "Obese"))
    md$whtr_ow_obese_binary <- factor(dplyr::case_when(
      md$whtr_obesity == "Normal"                     ~ "Normal",
      md$whtr_obesity %in% c("Overweight", "Obese")   ~ "Overweight_Obese",
      TRUE ~ NA_character_), levels = c("Normal", "Overweight_Obese"))
  }

  obj@meta.data <- md

  if (!nebula_cfg$offset_col %in% names(obj@meta.data))
    stop("Offset column '", nebula_cfg$offset_col, "' missing from ", dataset,
         " — NEBULA needs per-cell size factors (pooled_offset).")
  obj
}

# ---------------------------------------------------------------------------
# apply_subset(): cohort filter -> baseline filter -> transform -> drop NA.
# Returns the subsetted object (or NULL if nothing left).
# ---------------------------------------------------------------------------
apply_subset <- function(obj, contrast, subset_spec) {
  # cohort filter on `group`
  if (!is.null(subset_spec$cohort)) {
    if (!"group" %in% names(obj@meta.data))
      stop("cohort filter requested but no 'group' column present")
    keep <- obj@meta.data$group %in% subset_spec$cohort
    if (!any(keep)) return(NULL)
    obj <- subset(obj, cells = colnames(obj)[keep])
  }
  # baseline (PRE) only — honored at either contrast or subset level
  if (isTRUE(subset_spec$baseline_only) || isTRUE(contrast$baseline_only)) {
    keep <- obj@meta.data$visit %in% "PRE"
    if (!any(keep)) return(NULL)
    obj <- subset(obj, cells = colnames(obj)[keep])
  }
  # transform -> new modelled column, drop non-finite
  if (!is.null(contrast$transform)) {
    fn <- match.fun(contrast$transform)
    obj@meta.data[[contrast$var]] <- fn(obj@meta.data[[contrast$source_var]])
    fin <- is.finite(obj@meta.data[[contrast$var]])
    if (!any(fin)) return(NULL)
    obj <- subset(obj, cells = colnames(obj)[fin])
  }
  obj
}

# ---------------------------------------------------------------------------
# relevel_predictors(): set reference levels so NEBULA estimates are vs the
# intended reference. Returns the object plus a human-readable reference label.
# ---------------------------------------------------------------------------
relevel_predictors <- function(obj, contrast) {
  md <- obj@meta.data
  ref_label <- NA_character_

  if (identical(contrast$type, "categorical")) {
    v <- as.character(md[[contrast$var]])
    f <- factor(v)
    if (!is.null(contrast$reference) && contrast$reference %in% levels(f)) {
      f <- relevel(f, ref = contrast$reference)
      ref_label <- contrast$reference
    } else {
      ref_label <- levels(f)[1]
      if (!is.null(contrast$reference))
        warning("reference '", contrast$reference, "' not found for ",
                contrast$var, "; using '", ref_label, "'")
    }
    md[[contrast$var]] <- f

  } else if (identical(contrast$type, "did")) {
    parts <- trimws(strsplit(contrast$var, "[*:]")[[1]])  # treatment, visit
    tvar  <- parts[1]; vvar <- parts[2]
    tf <- factor(as.character(md[[tvar]]))
    vf <- factor(as.character(md[[vvar]]))
    if (!is.null(contrast$treatment_ref) && contrast$treatment_ref %in% levels(tf))
      tf <- relevel(tf, ref = contrast$treatment_ref)
    if (!is.null(contrast$visit_ref) && contrast$visit_ref %in% levels(vf))
      vf <- relevel(vf, ref = contrast$visit_ref)
    md[[tvar]] <- tf; md[[vvar]] <- vf
    ref_label <- paste0(tvar, "=", levels(tf)[1], "; ", vvar, "=", levels(vf)[1])
  }
  # continuous -> ref_label stays NA

  obj@meta.data <- md
  list(obj = obj, reference = ref_label)
}

# ---------------------------------------------------------------------------
# select_genes(): untargeted gene set = genes expressed in enough cells.
# ---------------------------------------------------------------------------
select_genes <- function(counts, min_cells, min_frac) {
  n        <- ncol(counts)
  n_expr   <- Matrix::rowSums(counts > 0)
  thresh   <- max(min_cells, ceiling(min_frac * n))
  rownames(counts)[n_expr >= thresh]
}

# ---------------------------------------------------------------------------
# run_nebula(): fit NEBULA on the given genes. Returns the raw nebula result.
# ---------------------------------------------------------------------------
run_nebula <- function(obj, form, record_id, features, nebula_cfg, ncore) {
  md        <- obj@meta.data
  vars_used <- all.vars(form)
  keep      <- stats::complete.cases(md[, vars_used, drop = FALSE])
  if (sum(keep) < 2 || dplyr::n_distinct(md[[record_id]][keep]) < 2) return(NULL)

  md_v     <- droplevels(md[keep, , drop = FALSE])
  counts   <- GetAssayData(obj, layer = "counts")
  counts_v <- round(counts[features, keep, drop = FALSE])
  offset_v <- md[[nebula_cfg$offset_col]][keep]
  pred     <- model.matrix(form, data = md_v)

  dg <- group_cell(count = counts_v, id = md_v[[record_id]],
                   pred = pred, offset = offset_v)
  if (is.null(dg)) dg <- list(count = counts_v, id = md_v[[record_id]],
                              pred = pred, offset = offset_v)

  nebula(count = dg$count, id = dg$id, pred = dg$pred, offset = dg$offset,
         ncore = ncore, reml = nebula_cfg$reml, model = nebula_cfg$model)
}

# ---------------------------------------------------------------------------
# count_stats(): n_cells / n_people / median cells-per-person, plus per-group
# cell + people counts for categorical / DiD predictors.
# ---------------------------------------------------------------------------
count_stats <- function(obj, var, person_col, visit_col = NULL) {
  md     <- obj@meta.data
  people <- md[[person_col]]

  if (!is.null(visit_col) && visit_col %in% names(md)) {
    per_person <- as.integer(table(interaction(people, md[[visit_col]], drop = TRUE)))
  } else {
    per_person <- as.integer(table(people))
  }

  stats <- data.frame(
    n_cells                 = ncol(obj),
    n_people                = dplyr::n_distinct(people),
    median_cells_per_person = stats::median(per_person)
  )

  terms      <- trimws(strsplit(var, "[*:+]")[[1]])
  terms      <- terms[terms %in% names(md)]
  group_vars <- terms[!vapply(terms, function(t) is.numeric(md[[t]]), logical(1))]
  if (length(group_vars) > 0) {
    grp       <- do.call(interaction, c(md[group_vars], list(drop = TRUE, sep = "_")))
    cells_by  <- table(grp)
    people_by <- tapply(people, grp, dplyr::n_distinct)
    for (g in names(cells_by)) {
      stats[[paste0("cells_",  g)]] <- as.integer(cells_by[[g]])
      stats[[paste0("people_", g)]] <- as.integer(people_by[[g]])
    }
  }
  stats
}

# ---------------------------------------------------------------------------
# run_one(): the full single-run driver used by the worker.
#   meta describes the run (contrast/subset/resolution/celltype) for the output.
#   Returns a data.frame (one row per gene) or NULL if the run was skipped.
# ---------------------------------------------------------------------------
run_one <- function(obj, contrast, subset_key, subset_spec, resolution_name,
                    ct_col, celltype, dataset_cfg, nebula_cfg, ncore,
                    dataset_name = contrast$dataset) {

  obj <- apply_subset(obj, contrast, subset_spec)
  if (is.null(obj) || ncol(obj) == 0) { message("  skip: no cells after subset"); return(NULL) }

  # skip gracefully if the predictor column(s) are absent in this dataset
  need <- all.vars(as.formula(paste("~", contrast$var)))
  miss <- setdiff(need, names(obj@meta.data))
  if (length(miss) > 0) {
    message("  skip: column(s) not in metadata: ", paste(miss, collapse = ", "))
    return(NULL)
  }

  rl   <- relevel_predictors(obj, contrast)
  obj  <- rl$obj
  ref  <- rl$reference

  # run-level size gates
  n_people <- dplyr::n_distinct(obj@meta.data[[dataset_cfg$person_col]])
  if (ncol(obj) < nebula_cfg$min_total_cells ||
      n_people  < nebula_cfg$min_total_people) {
    message("  skip: too small (cells=", ncol(obj), ", people=", n_people, ")")
    return(NULL)
  }

  # untargeted gene selection
  counts   <- GetAssayData(obj, layer = "counts")
  features <- select_genes(counts, nebula_cfg$min_cells, nebula_cfg$min_frac)
  if (length(features) == 0) { message("  skip: no genes pass expression filter"); return(NULL) }

  form <- as.formula(paste("~", contrast$var))
  res  <- tryCatch(
    run_nebula(obj, form, dataset_cfg$record_id, features, nebula_cfg, ncore),
    error = function(e) { message("  nebula error: ", conditionMessage(e)); NULL })
  if (is.null(res)) return(NULL)

  summ <- res$summary
  # NEBULA may internally drop low-expression genes (its `cpc` filter), so the
  # summary can have fewer rows than `features`. Trust NEBULA's own gene labels;
  # only fall back to `features` when row counts actually match.
  if (is.null(summ$gene)) {
    summ$gene <- if (nrow(summ) == length(features)) features else seq_len(nrow(summ))
  }
  # convergence + overdispersion are returned per retained gene, in summary order
  summ$convergence <- res$convergence
  if (!is.null(res$overdispersion)) {
    od <- as.data.frame(res$overdispersion)
    names(od) <- paste0("overdisp_", names(od))
    summ <- cbind(summ, od)
  }

  stats <- count_stats(obj, contrast$var, dataset_cfg$person_col, dataset_cfg$visit_col)

  desc <- data.frame(
    contrast      = contrast$name,
    var           = contrast$var,
    type          = contrast$type,
    dataset       = dataset_name,
    subset        = subset_key,
    subset_label  = subset_spec$label,
    resolution    = resolution_name,
    celltype_col  = ct_col,
    celltype      = celltype,
    reference_group = ref,
    model_formula = paste(deparse(form), collapse = ""),
    n_genes_tested = length(features),
    stringsAsFactors = FALSE
  )

  out <- cbind(
    desc[rep(1, nrow(summ)), , drop = FALSE],
    summ,
    stats[rep(1, nrow(summ)), , drop = FALSE],
    row.names = NULL
  )
  message("  ok: ", contrast$name, " | ", subset_key, " | ", resolution_name,
          " | ", celltype, " | genes=", length(features),
          " | cells=", stats$n_cells, " | people=", stats$n_people)
  out
}
