# =============================================================================
# config.R  —  Untargeted (whole-transcriptome) NEBULA pipeline configuration
# =============================================================================
#   - global  : S3 locations, local work dir, NEBULA parameters, SLURM settings
#   - datasets: where to load each Seurat object and its id/visit columns
#   - resolutions: the two cell-type columns (general / specific)
#   - contrasts: every model to fit. Each contrast -> ONE final CSV that holds
#                all subsets + all cell types + both resolutions.
#
# A "run" (one SLURM array task / one NEBULA fit) is the cross product:
#     contrast  x  subset  x  resolution  x  cell type
# NEBULA is fit on ALL (sufficiently expressed) genes at once, so the gene set
# is untargeted — there is no gene list to maintain.
#
# Contrast field reference
#   dataset     : "pb90" | "attempt"  (must match a name in `datasets`)
#   type        : "categorical" | "continuous" | "did"
#   var         : model predictor. For "did" use "treatment*visit".
#   reference   : (categorical) the reference LEVEL of `var`. Recorded in output.
#   treatment_ref/visit_ref : (did) reference levels for the two factors.
#   source_var  : (continuous w/ transform) raw column the transform reads.
#   transform   : (continuous) name of a function applied to source_var, e.g.
#                 "log". Leave NULL to model `var` directly.
#   baseline_only : (attempt continuous) TRUE -> keep visit == "PRE" only.
#   subsets     : named list of cohort filters. Each subset:
#                   cohort : character vector of `group` levels to keep,
#                            or NULL for everyone.
#                   label  : human-readable description written to the output.
# =============================================================================

CONFIG <- list(

  # ---- global settings ------------------------------------------------------
  global = list(
    # local scratch dir (on /gscratch) for the job manifest + logs only.
    # The big data (splits) and the results live on S3 (see below).
    work_dir = "/mmfs1/gscratch/togo/yejichoi/project_logs/untargeted_nebula_logs",

    # S3 bucket + key prefixes (everything stored on S3 per project choice)
    s3_bucket      = "togo.projects",
    splits_prefix  = "P047_scRNA_transcript_search_yc/clean_data/splits/",
    partials_prefix= "P047_scRNA_transcript_search_yc/results/temp_csv/untargeted_nebula/partials/",
    results_prefix = "P047_scRNA_transcript_search_yc/results/temp_csv/untargeted_nebula/",

    # NEBULA parameters
    nebula = list(
      model     = "NBLMM",   # negative-binomial lognormal mixed model
      reml      = 1,
      ncore     = 4,         # overridden by SLURM_CPUS_PER_TASK when present
      offset_col= "pooled_offset",
      # untargeted gene filtering (applied per subset x cell type, before fit):
      # keep a gene if it is expressed (count > 0) in at least min_cells cells
      # AND in at least min_frac of the cells in that run.
      min_cells = 5,
      min_frac  = 0.01,
      # run-level gates: skip a (subset x cell type) fit if too small
      min_total_cells  = 20,
      min_total_people = 5
    ),

    # SLURM defaults (used by the .slurm scripts via envsubst-style comments)
    slurm = list(
      account   = "togo",
      partition = "cpu-g2",
      cpus      = 4,
      mem       = "120G",
      time_split   = "08:00:00",
      time_nebula  = "12:00:00",
      time_aggregate = "02:00:00"
    )
  ),

  # ---- cell-type resolutions (low res / high res) ---------------------------
  resolutions = list(
    general  = "KPMP_celltype_general",  # low resolution
    specific = "KPMP_celltype"           # high resolution
  ),

  # ---- datasets -------------------------------------------------------------
  datasets = list(
    pb90 = list(
      s3_bucket  = "core.data",
      s3_object  = "scRNA/Seurat/pb90_processed_full.rds",
      person_col = "record_id",
      record_id  = "record_id",
      visit_col  = NULL
    ),
    attempt = list(
      s3_bucket  = "core.data",
      s3_object  = "scRNA/Seurat/attempt_processed_full.rds",
      person_col = "subject_id",
      record_id  = "subject_id",
      visit_col  = "visit"
    )
  ),

  # ---- contrasts (each -> one CSV) ------------------------------------------
  contrasts = list(

    ## --- PB90 -------------------------------------------------------------

    # Disease group, Lean Control as reference (all cohorts)
    group = list(
      dataset   = "pb90",
      type      = "categorical",
      var       = "group",
      reference = "Lean Control",
      subsets   = list(
        all = list(cohort = NULL, label = "All cohorts")
      )
    ),

    # SGLT2i ever (yes/no), T2D only
    sglt2i_ever = list(
      dataset   = "pb90",
      type      = "categorical",
      var       = "sglt2i_ever",
      reference = "No",
      subsets   = list(
        t2d = list(cohort = "Type 2 Diabetes", label = "T2D only")
      )
    ),

    # SGLT2i at timepoint (yes/no), T2D only
    sglt2i_timepoint = list(
      dataset   = "pb90",
      type      = "categorical",
      var       = "sglt2i_timepoint",
      reference = "No",
      subsets   = list(
        t2d = list(cohort = "Type 2 Diabetes", label = "T2D only")
      )
    ),

    # DKD (eGFR<90 or uACR>=30), T2D only; non_DKD reference
    dkd30 = list(
      dataset   = "pb90",
      type      = "categorical",
      var       = "dkd_group_30",
      reference = "non_DKD",
      subsets   = list(
        t2d = list(cohort = "Type 2 Diabetes", label = "T2D only")
      )
    ),

    # DKD (eGFR<90 or uACR>=100), T2D only; non_DKD reference
    dkd100 = list(
      dataset   = "pb90",
      type      = "categorical",
      var       = "dkd_group_100",
      reference = "non_DKD",
      subsets   = list(
        t2d = list(cohort = "Type 2 Diabetes", label = "T2D only")
      )
    ),

    # Arteriosclerosis present (any grade) vs none; derived from
    # arteriosclerosis_sev (0 -> "No", >=1 -> "Yes"). PB90, all cohorts.
    arteriosclerosis_yn = list(
      dataset   = "pb90",
      type      = "categorical",
      var       = "arteriosclerosis_yn",
      reference = "No",
      subsets   = list(
        all = list(cohort = NULL, label = "All cohorts")
      )
    ),

    # Arteriolar hyalinosis present (any grade) vs none; derived from
    # arteriolohyalinosis_sev (0 -> "No", >=1 -> "Yes"). PB90, all cohorts.
    arteriolohyalinosis_yn = list(
      dataset   = "pb90",
      type      = "categorical",
      var       = "arteriolohyalinosis_yn",
      reference = "No",
      subsets   = list(
        all = list(cohort = NULL, label = "All cohorts")
      )
    ),

    # log(uACR) association — everyone AND T2D only (both kept in one CSV)
    log_uacr = list(
      dataset    = "pb90",
      type       = "continuous",
      var        = "log_acr_u",
      source_var = "acr_u",
      transform  = "log",
      subsets    = list(
        all = list(cohort = NULL,              label = "All cohorts"),
        t2d = list(cohort = "Type 2 Diabetes", label = "T2D only")
      )
    ),

    # eGFR association — everyone AND T2D only (both kept in one CSV)
    egfr = list(
      dataset = "pb90",
      type    = "continuous",
      var     = "eGFR_CKD_epi",
      subsets = list(
        all = list(cohort = NULL,              label = "All cohorts"),
        t2d = list(cohort = "Type 2 Diabetes", label = "T2D only")
      )
    ),

    ## --- ATTEMPT (T1D RCT; PRE/POST; Dapagliflozin vs Placebo) -------------

    # Difference-in-differences: treatment x visit interaction.
    # Reference cell = Placebo & PRE. The interaction term is the DiD estimate.
    attempt_did = list(
      dataset       = "attempt",
      type          = "did",
      var           = "treatment*visit",
      treatment_ref = "Placebo",
      visit_ref     = "PRE",
      subsets       = list(
        all = list(cohort = NULL, label = "All ATTEMPT")
      )
    ),

    # eGFR association at baseline (PRE) only
    attempt_egfr_bl = list(
      dataset       = "attempt",
      type          = "continuous",
      var           = "eGFR_CKD_epi",
      baseline_only = TRUE,
      subsets       = list(
        pre = list(cohort = NULL, label = "PRE (baseline) only")
      )
    ),

    # log(uACR) association at baseline (PRE) only
    attempt_log_uacr_bl = list(
      dataset       = "attempt",
      type          = "continuous",
      var           = "log_uacr",
      source_var    = "emu_urine_acr_mean",
      transform     = "log",
      baseline_only = TRUE,
      subsets       = list(
        pre = list(cohort = NULL, label = "PRE (baseline) only")
      )
    )
  )
)

# =============================================================================
# Variable-driven continuous associations (CMR/imaging + labs + clinical)
# -----------------------------------------------------------------------------
# Each variable below becomes ONE contrast (-> one CSV) with three subsets, all
# kept in the same file:
#   - PB90, all cohorts
#   - PB90, T2D only
#   - ATTEMPT, baseline (PRE) visit only
# Subsets carry their own `dataset` (and `baseline_only`), so a single contrast
# spans both datasets. A run is skipped automatically (logged, no output) when
# the variable is absent from a given dataset's metadata or too few cells/people
# remain — so listing a variable that only exists in one cohort is harmless.
# `%||%` is defined in setup.R, which is always sourced before this file.
# =============================================================================

.assoc_subsets <- list(
  pb90_all    = list(dataset = "pb90",    cohort = NULL,              label = "PB90 all cohorts"),
  pb90_t2d    = list(dataset = "pb90",    cohort = "Type 2 Diabetes", label = "PB90 T2D only"),
  attempt_pre = list(dataset = "attempt", baseline_only = TRUE,       label = "ATTEMPT baseline (PRE)")
)

.assoc_vars <- c(
  # ---- CMR / vascular imaging (im_labels; metadata-present columns only) ---
  "lv_hr", "rvedv", "lvedv", "rvesv", "lvesv",
  "rv_stroke_volume", "lv_stroke_volume",
  "rv_cardiac_output", "lv_cardiac_output",
  "lves_mass", "lved_mass",
  "lv_myo_thickness_total", "lv_myo_thickness_dias",
  "rvef", "lvef", "gls", "gcs", "grs",
  "af_high5pctwss", "af_medianwss", "af_meanwss", "af_pwv",
  "tt1", "tt2", "mv_flow_e", "mv_flow_a", "mv_flow_ea",

  # ---- glycemic & insulin resistance --------------------------------------
  "hba1c", "fasting_glu", "homa_ir", "m_i", "fasting_insulin", "fasting_cpep",

  # ---- adiposity / lipids / blood pressure --------------------------------
  "bmi", "triglycerides", "hdl", "ldl", "cholesterol",
  "waist_hip_ratio", "vat", "leptin_base", "sbp", "dbp", "map",

  # ---- renal hemodynamics (clearance-based) -------------------------------
  "gfr_bsa_plasma", "erpf_bsa_plasma", "ff", "rvr", "glomerular_pressure",

  # ---- kidney biopsy morphometry ------------------------------------------
  "gbm_thick_artmean", "mes_index", "pod_nuc_density", "ifta", "glom_volume_con",

  # ---- additional renal labs ----------------------------------------------
  "creatinine_s", "cystatin_c_s", "bun_s", "sua", "copeptin",
  "hematocrit", "hemoglobin",

  # ---- kidney functional MRI (cortex/kidney perfusion + relaxometry) -------
  "avg_c_k1", "avg_c_k2", "avg_c_r2", "avg_c_t1", "avg_c_f",  # cortex
  "adc_left", "adc_right",                                    # diffusion (ADC)
  "avg_k_fsoc", "avg_k_r2", "avg_k_t1",                       # whole kidney
  "avg_pcascl"                                                # ASL perfusion
)

local({
  assoc <- setNames(
    lapply(.assoc_vars, function(v) list(
      dataset = "pb90",          # default; each subset overrides with its own
      type    = "continuous",
      var     = v,
      subsets = .assoc_subsets
    )),
    paste0("assoc_", .assoc_vars)
  )
  CONFIG$contrasts <<- c(CONFIG$contrasts, assoc)
})

# =============================================================================
# Obesity classifications — categorical group differences
# -----------------------------------------------------------------------------
# Derived in derive_meta() (see nebula_core.R) from age/bmi/bmip (bmi_obesity),
# dexa_body_fat/sex (dxa_obesity), and waistcm/height (whtr_* measures). Each is
# a categorical contrast run on PB90 (all cohorts) and ATTEMPT baseline (PRE).
# A run is skipped where the required inputs are missing in that dataset.
# `reference` is the level all others are compared against.
# =============================================================================
.obesity_subsets <- list(
  pb90_all    = list(dataset = "pb90",    cohort = NULL,        label = "PB90 all cohorts"),
  attempt_pre = list(dataset = "attempt", baseline_only = TRUE, label = "ATTEMPT baseline (PRE)")
)

.obesity_refs <- list(
  bmi_obesity          = "Normal",      # Normal / Overweight / Obese
  dxa_obesity          = "Normal",      # Normal / Overweight / Obese
  whtr_obesity         = "Normal",      # Normal / Overweight / Obese (waist-to-height)
  whtr_obese_binary    = "Non_Obese",   # Non_Obese / Obese
  whtr_ow_obese_binary = "Normal"       # Normal / Overweight_Obese
)

local({
  obes <- setNames(
    lapply(names(.obesity_refs), function(v) list(
      dataset   = "pb90",
      type      = "categorical",
      var       = v,
      reference = .obesity_refs[[v]],
      subsets   = .obesity_subsets
    )),
    names(.obesity_refs)
  )
  CONFIG$contrasts <<- c(CONFIG$contrasts, obes)
})
