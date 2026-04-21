# T1D Adiposity Analysis

This project investigates the relationship between adiposity and kidney single-cell RNA sequencing (scRNA-seq) gene expression in youth with Type 1 Diabetes (T1D) compared to lean healthy controls (HC). Data are drawn from the CROCODILE, ATTEMPT, CASPER, RENAL-HEIR, and PANDA studies. Obesity is classified by both BMI percentile and DXA body fat percentage. Differential gene expression is modeled using NEBULA, and results are visualized through volcano plots, GSEA pathway analyses, and interactive Quarto reports.

All data are stored and retrieved from an S3 bucket (`t1d.adiposity`) hosted at `s3.kopah.uw.edu`. AWS credentials are read from a local `keys.json` file, and user-specific paths are resolved automatically at the top of each script.

## Pipeline overview

The scripts are numbered in the order they should be run. Steps 01–03 handle data preparation and exploratory visualization. Steps 04–05 are standalone analysis notebooks. Step 06 compiles NEBULA run metadata after HPC jobs complete. Step 07 generates publication figures from the compiled results, and step 08 produces a slide-based results report.

## Scripts

### 01_clean_dataset.R

Reads the harmonized SomaScan clinical dataset from the shared drive, filters to T1D and lean HC participants across the five studies, restricts to baseline visits, and creates derived obesity categories (Normal / Overweight / Obese) based on both BMI percentile and DXA body fat cutoffs. Saves two cleaned RDS files to S3: one with SomaScan aptamer columns (`t1d_hc_clinical_data_soma.csv`) and one without (`t1d_hc_clinical_data.csv`).

**When to run:** Once at the start of the project, or whenever the upstream harmonized dataset is updated.

### 02_meta_merge_w_seurat.R

Loads the cleaned clinical data and the integrated Seurat scRNA-seq object (`pb90_attempt_integrated_processed.rds` from the `scrna` bucket). Merges clinical metadata onto the Seurat object's cell-level metadata by `record_id` and `kit_id`, creates the `KPMP_celltype_general` grouping variable, and saves the merged Seurat object to S3 (`t1d_hc_scrna_w_clinical.rds`).

**When to run:** After `01_clean_dataset.R`, or whenever the Seurat object or clinical data is updated.

### 03_scrna_visualizations.R

Generates exploratory scRNA-seq figures from the merged Seurat object, including UMAP plots colored by KPMP cell type and by data source (ATTEMPT vs. PB90), and stacked bar charts of cell counts stratified by disease group, DXA obesity, and BMI obesity. All figures are uploaded to S3 under `results/figures/`.

**When to run:** After `02_meta_merge_w_seurat.R`. Re-run to refresh exploratory figures after any data updates.

### 04_nonbx_analysis.qmd

A Quarto notebook that performs descriptive clinical analyses on the full cohort (not limited to participants with biopsies). Produces summary tables of demographics and clinical measures stratified by BMI and DXA obesity groups, tabulates data availability across key outcomes (mGFR, RPF, Gomez hemodynamics, mpMRI, CGM, blood pressure, PWV), and generates preliminary boxplots comparing outcomes across obesity strata.

**When to run:** Interactively as needed for clinical characterization and data availability assessment.

### 05_proteomics_analysis.qmd

A Quarto notebook stub for SomaScan proteomics analysis. Loads the clinical dataset with aptamer columns from S3 and is intended for differential protein expression analyses.

**When to run:** Interactively, once the proteomics analysis plan is finalized.

### 06_compile_metadata.R

Scans all NEBULA result subdirectories on S3, collects per-run metadata CSVs (one per analysis-by-cell-type combination), and compiles them into a single summary file (`results/nebula/csv/run_metadata_compiled.csv`). Reports which analysis types have completed and which are still pending.

**When to run:** After all (or a batch of) NEBULA SLURM array jobs have finished. Can be re-run incrementally as more jobs complete.

### 07_generate_figures.R

Reads the compiled NEBULA metadata and per-cell-type result CSVs from S3, then generates volcano plots (at both nominal p < 0.05 and FDR < 0.05 thresholds) and GSEA Hallmark pathway lollipop plots for every completed analysis-by-cell-type combination. Figures are saved locally under `Results/Figures/`.

**When to run:** After `06_compile_metadata.R`. Run locally (not on HPC) since it reads from S3 and writes figures to a local directory.

### 08_t1d_adiposity_scRNA_results_vis.qmd

A Quarto RevealJS slide deck that renders an interactive results report for the NEBULA analyses. For each analysis type, it displays volcano plots and pathway enrichment figures organized by cell type group (PT, TAL, immune, endothelial, etc.), pulling pre-generated images from the local `Results/Figures/` directory.

**When to run:** After `07_generate_figures.R`. Render with `quarto render` to produce the HTML slide deck for presentation or sharing.

## nebula/ — HPC job infrastructure

This subfolder contains everything needed to run the NEBULA differential expression analyses as SLURM array jobs on the Hyak cluster.

### nebula/scripts/save_celltype_subsets.R

Pre-splits the full merged Seurat object into per-cell-type RDS files on S3 (`data_clean/subset/`), so each array job loads only the cells it needs rather than the entire object. Subsets are created for both individual KPMP cell types (~44 types) and grouped general cell types (~15 groups).

**When to run:** Once after `02_meta_merge_w_seurat.R`, before submitting any SLURM array jobs.

### nebula/scripts/generate_job_config.R

Generates the tab-delimited config file (`config/job_config.txt`) that maps each SLURM array index to a specific analysis type, cell type, and cell type variable. Covers categorical comparisons (unadjusted and age-adjusted), continuous adiposity measures (T1D-only and all-subjects, unadjusted and adjusted), and interaction models.

**When to run:** Once before the first SLURM submission. Re-run if analysis types or cell type lists change.

### nebula/scripts/generate_job_config_adj_age_sex.R

Same as above but generates a config file specifically for age + sex adjusted models (both categorical and continuous).

**When to run:** Once before submitting the age + sex adjusted batch of SLURM jobs.

### nebula/scripts/generate_job_config_normal_vs_ow_obese.R

Generates a focused config file for Normal vs. Overweight+Obese comparisons only (BMI and DXA, across adjustment levels).

**When to run:** Once before submitting this specific subset of SLURM jobs.

### nebula/scripts/run_nebula_single.R

The core analysis script executed by each SLURM array task. Takes three command-line arguments (analysis type, cell type, cell type variable), loads the appropriate cell type subset from S3, builds the NEBULA model (negative binomial mixed model), and saves per-gene results and run metadata back to S3 under `results/nebula/<analysis_type>/`.

**When to run:** Called automatically by the SLURM scheduler; not run manually.

### nebula/nebula_array.slurm (and checkpoint variants)

SLURM batch submission scripts that launch the array jobs. Each task reads its parameters from the config file and calls `run_nebula_single.R` inside an Apptainer container (`YC_scRNA.sif`). The checkpoint variants (`_ckpt`, `_ckpt2`, `_ckpt3`, `_ckpt_age_sex`) are resumption scripts for restarting failed or incomplete batches.

**When to run:** Submit with `sbatch nebula_array.slurm` after generating the config file and saving cell type subsets.

## End-to-end workflow

1. `01_clean_dataset.R` — clean and upload clinical data
2. `02_meta_merge_w_seurat.R` — merge clinical metadata with Seurat object
3. `03_scrna_visualizations.R` — generate exploratory UMAP and cell count figures
4. `nebula/scripts/save_celltype_subsets.R` — pre-split Seurat object by cell type
5. `nebula/scripts/generate_job_config.R` (and variants) — create SLURM config files
6. `sbatch nebula/nebula_array.slurm` — submit NEBULA jobs to Hyak
7. `06_compile_metadata.R` — compile run metadata after jobs finish
8. `07_generate_figures.R` — generate volcano and pathway plots
9. `08_t1d_adiposity_scRNA_results_vis.qmd` — render interactive results report
10. `04_nonbx_analysis.qmd` / `05_proteomics_analysis.qmd` — run as needed for clinical and proteomics analyses
