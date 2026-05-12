# =============================================================================
# ATTEMPT MS - Supplemental GSEA tables (Fig S5-S9)
# -----------------------------------------------------------------------------
# Reads the per-cell-type GSEA RDS files:
#   <root_path>/ATTEMPT/Results/GSEA/full_<celltype>_gsea.RDS
# and writes one xlsx with 15 sheets covering Fig S5-S9
# (S5=PT, S6=TAL, S7=EC, S8=IC, S9=POD; panels A=reactome, B=hallmark, C=go).
#
# Output:
#   <root_path>/ATTEMPT/Results/GSEA/supplemental_gsea_tables.xlsx
# =============================================================================

suppressPackageStartupMessages({
  library(openxlsx)
  library(dplyr)
})

# -----------------------------------------------------------------------------
# Resolve root_path (same logic as ATTEMPT_supplemental_stats.R)
# -----------------------------------------------------------------------------
user <- Sys.info()[["user"]]
if (user == "choiyej") {
  root_path <- "/Users/choiyej/Library/CloudStorage/OneDrive-UW/Bjornstad/Biostatistics Core Shared Drive"
} else if (user == "pylell") {
  root_path <- "/Users/pylell/Library/CloudStorage/OneDrive-SharedLibraries-UW/Bjornstad/Biostatistics Core Shared Drive"
} else {
  stop("Unknown user: please specify root_path for this user.")
}

gsea_dir <- file.path(root_path, "ATTEMPT", "Results", "GSEA")
out_xlsx <- file.path(gsea_dir, "supplemental_gsea_tables.xlsx")

# -----------------------------------------------------------------------------
# Sheet plan: figure -> cell type, panel -> reference
# -----------------------------------------------------------------------------
celltype_map <- list(
  S5 = "PT",
  S6 = "TAL",
  S7 = "EC",
  S8 = "IC",
  S9 = "POD"
)
panel_map <- list(
  A = "reactome",
  B = "hallmark",
  C = "go"
)

# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------

# Coerce $reactome/$hallmark/$go element into a writeable data frame.
# Some pipelines return a data frame directly; others wrap results in a list
# (e.g., gsea_list$reactome$result). Try both, fall back to as.data.frame.
to_table <- function(obj) {
  if (is.null(obj)) return(NULL)
  if (is.data.frame(obj)) return(obj)
  if (is.list(obj)) {
    # Common nested layouts
    if (!is.null(obj$result) && is.data.frame(obj$result)) return(obj$result)
    if (!is.null(obj$table)  && is.data.frame(obj$table))  return(obj$table)
    # If it's an fgsea result list-of-data.frames, bind
    dfs <- Filter(is.data.frame, obj)
    if (length(dfs) == 1) return(dfs[[1]])
    if (length(dfs) > 1)  return(bind_rows(dfs, .id = "source"))
  }
  tryCatch(as.data.frame(obj), error = function(e) NULL)
}

# Build a 31-char-safe sheet name. openxlsx errors out beyond that limit.
make_sheet_name <- function(fig, panel, cell, ref) {
  nm <- paste0("Fig", fig, panel, "_", cell, "_", ref)
  if (nchar(nm) > 31) nm <- substr(nm, 1, 31)
  nm
}

# -----------------------------------------------------------------------------
# Build workbook
# -----------------------------------------------------------------------------
wb <- createWorkbook()

# README mapping sheet
readme_rows <- list()
header_style <- createStyle(textDecoration = "bold", fontColour = "#FFFFFF",
                            fgFill = "#4472C4", halign = "center",
                            border = "TopBottomLeftRight")

for (fig in names(celltype_map)) {
  cell      <- celltype_map[[fig]]
  rds_path  <- file.path(gsea_dir, paste0("full_", tolower(cell), "_gsea.RDS"))

  if (!file.exists(rds_path)) {
    message("Missing file, skipping: ", rds_path)
    next
  }
  gsea_list <- tryCatch(readRDS(rds_path), error = function(e) {
    message("readRDS failed for ", rds_path, ": ", conditionMessage(e))
    NULL
  })
  if (is.null(gsea_list)) next

  for (panel in names(panel_map)) {
    ref <- panel_map[[panel]]
    raw <- gsea_list[[ref]]
    df  <- to_table(raw)
    if (is.null(df) || nrow(df) == 0) {
      message("Empty/missing $", ref, " in ", basename(rds_path),
              " - sheet will be skipped.")
      next
    }
    sheet_name <- make_sheet_name(fig, panel, cell, ref)

    addWorksheet(wb, sheet_name)
    writeData(wb, sheet_name, df, headerStyle = header_style)
    freezePane(wb, sheet_name, firstRow = TRUE)
    setColWidths(wb, sheet_name, cols = seq_len(ncol(df)), widths = "auto")

    readme_rows[[length(readme_rows) + 1]] <- data.frame(
      Sheet      = sheet_name,
      Figure     = paste0("Fig ", fig, panel),
      Cell_type  = cell,
      Reference  = ref,
      Rows       = nrow(df),
      Source_RDS = basename(rds_path),
      stringsAsFactors = FALSE
    )
  }
}

# README first
readme_df <- if (length(readme_rows)) bind_rows(readme_rows) else
  data.frame(Sheet = character(), Figure = character(), Cell_type = character(),
             Reference = character(), Rows = integer(),
             Source_RDS = character(), stringsAsFactors = FALSE)

# Insert README at position 1
addWorksheet(wb, "README")
writeData(wb, "README", readme_df, headerStyle = header_style)
freezePane(wb, "README", firstRow = TRUE)
setColWidths(wb, "README", cols = 1:6, widths = c(28, 12, 12, 12, 8, 32))
# Move README to front
worksheetOrder(wb) <- c(which(names(wb) == "README"),
                        setdiff(seq_along(names(wb)), which(names(wb) == "README")))

# -----------------------------------------------------------------------------
# Save
# -----------------------------------------------------------------------------
if (!dir.exists(gsea_dir)) {
  stop("GSEA directory does not exist: ", gsea_dir)
}
saveWorkbook(wb, file = out_xlsx, overwrite = TRUE)
message("Wrote: ", out_xlsx)
message("Sheets: ", paste(names(wb), collapse = ", "))
