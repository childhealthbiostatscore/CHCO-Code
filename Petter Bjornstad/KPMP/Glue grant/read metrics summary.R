# Preparing files for KPMP glue grant

library(jsonlite)
library(aws.s3)

user <- Sys.info()[["user"]]

if (user == "choiyej") {
  root_path <- "/Users/choiyej/Library/CloudStorage/OneDrive-SharedLibraries-UW/Laura Pyle - Bjornstad/Biostatistics Core Shared Drive/"
  git_path <- "/Users/choiyej/GitHub/CHCO-Code/Petter Bjornstad/"
} else if (user == "pylell") {
  root_path <- "/Users/pylell/Library/CloudStorage/OneDrive-SharedLibraries-UW/Bjornstad/Biostatistics Core Shared Drive/"
  git_path <- "/Users/pylell/Documents/GitHub/CHCO-Code/Petter Bjornstad/"
} else {
  stop("Unknown user: please specify root path for this user.")
}

keys <- fromJSON("/Users/choiyej/Library/CloudStorage/OneDrive-TheUniversityofColoradoDenver/Bjornstad Pyle Lab/keys.json")

Sys.setenv(
  "AWS_ACCESS_KEY_ID" = keys$MY_ACCESS_KEY,
  "AWS_SECRET_ACCESS_KEY" = keys$MY_SECRET_KEY,
  "AWS_DEFAULT_REGION" = "",
  "AWS_REGION" = "",
  "AWS_S3_ENDPOINT" = "s3.kopah.uw.edu"
)

# Read folder names in kopah
res <- get_bucket_df("scrna", region = "")
keys <- res$Key[grepl("^2025-Sept-18/", res$Key)]

# strip prefix
subpaths <- sub("^2025-Sept-18/", "", keys)

# take only first folder part
subfolders <- unique(sub("/.*", "", subpaths))

# remove empty ones (these were files directly in the folder)
subfolders <- subfolders[subfolders != ""]

s3read_using_region <- function(FUN, ..., object, bucket, region = NULL, opts = NULL, filename = NULL) {
  if (missing(bucket)) {
    bucket <- get_bucketname(object)
  }
  object <- get_objectkey(object)
  
  tmp <- if (is.character(filename)) {
    file.path(tempdir(TRUE), filename)
  } else {
    tempfile(fileext = paste0(".", tools::file_ext(object)))
  }
  
  on.exit(unlink(tmp))
  
  # Add region to opts if provided
  if (!is.null(region)) {
    if (is.null(opts)) {
      opts <- list(region = region)
    } else {
      opts$region <- region
    }
  }
  
  if (is.null(opts)) {
    r <- save_object(bucket = bucket, object = object, file = tmp)
  } else {
    r <- do.call("save_object", c(list(bucket = bucket, 
                                       object = object, 
                                       file = tmp), opts))
  }
  
  return(FUN(tmp, ...))
}

# Read metrics_summary.csv file from kopah folders
metrics_summary_compiled <- data.frame()

for (key in subfolders) {
  
  metrics_summary_df <- s3read_using_region(FUN = read.csv, 
                                     object = paste0("2025-Sept-18/", key, "/", key, ".metrics_summary.csv"), 
                                     bucket = "scrna",
                                     region = "")
  metrics_summary_df <- metrics_summary_df %>%
    mutate(Reads.Mapped.to.Genome = as.numeric(sub("%", "", Reads.Mapped.to.Genome)),
           unassigned.reads = ((100 - Reads.Mapped.to.Genome)/100) * as.numeric(gsub(",", "", Number.of.Reads)),
           uniquely.mapped.reads = (Reads.Mapped.to.Genome/100) * as.numeric(gsub(",", "", Number.of.Reads))
           )
  
  metrics_summary_df$run_name <- key
  
  
  metrics_summary_compiled <- rbind(metrics_summary_compiled, metrics_summary_df)
  
}

write.csv(metrics_summary_compiled, file.path(root_path, "KPMP/Glue Grant/metrics_summary_compiled.csv"), row.names = F)

# # FASTQ file names
# fastq_tbl <- tibble(s3_key = keys) %>%
#   filter(str_ends(s3_key, "\\.fastq\\.gz$")) %>%
#   mutate(
#     run_name  = str_match(s3_key, "^2025-Sept-18/([^/]+)/")[,2],
#     # file name only
#     filename  = basename(s3_key),
#     # R1 vs R2 from suffix
#     read      = str_match(filename, "_R([12])_001\\.fastq\\.gz$")[,2],
#     # common stem without the _R1/_R2 tail (good for pairing)
#     stem      = str_remove(filename, "_R[12]_001\\.fastq\\.gz$"),
#     # the “key” for sample pairing (the barcode bit between first '_' and _S#_)
#     sample_key = str_match(filename, "^[^_]+_(.+?)_S\\d+_R[12]_001\\.fastq\\.gz$")[,2]
#   ) %>%
#   dplyr::select(run_name, sample_key, stem, read, s3_key, filename)
# 
# # wide format with separate R1 / R2 columns
# fastq_pairs <- fastq_tbl %>%
#   pivot_wider(
#     names_from  = read,
#     values_from = c(s3_key, filename),
#     names_prefix = "R"
#   ) %>%
#   # Optional: order nicely
#   arrange(run_name, sample_key, stem)
# 
# fastq_pairs
