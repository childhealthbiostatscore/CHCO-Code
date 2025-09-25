suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(here))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(pbapply))

################################################################################
# Functions
source(here("R/identify_signals_from_summary_stats.R"))

get_top_snp = function(ii, infiles)
{
  sumstats_file = infiles[ii,"sumstats_file"]
  indata        = get_gwas_signals(sumstats_file, p_threshold = 1e-6) %>% select(-phenotype)
  
  return(suppressWarnings(cbind(infiles[ii,], indata)))
}

get_all_sumstat_files = function(progress, outfile)
{
  infiles = as.data.frame(rbindlist(lapply(progress, function(x)
  {
    ids            = x[["summary_stats"]]
    path           = x[["folder_path"  ]]
    analysis       = x[["folder"       ]]
    sumstats_files = paste(path, "summary_statistics", paste(ids, "summary_stats.rds", sep = "_"), sep = "/")

    return(data.frame(analysis = analysis, phenotype = ids, sumstats_file = sumstats_files))
  })), stringsAsFactors = FALSE)

  analyzed_phenotypes_n = 0
  
  if(file.exists(sumstats_report_file) == TRUE)
  {
    outdata = fread(sumstats_report_file, sep = "\t", header = TRUE, data.table = FALSE)
    
    infiles = infiles %>% filter(!sumstats_file %in% outdata$sumstats_file)
    
    analyzed_phenotypes_n = nrow(unique(outdata %>% select(analysis, phenotype)))
  }
  
  cat(paste0("Analyzed phenotypes: "  , analyzed_phenotypes_n, "\n"))
  cat(paste0("Phenotypes to analyze: ", nrow(infiles)        , "\n"))
  
  outdata_updated = as.data.frame(rbindlist(pblapply(1:nrow(infiles), function(ii){get_top_snp(ii, infiles)})), stringsAsFactors = FALSE)

  if(file.exists(sumstats_report_file) == TRUE)
  {
    outdata_updated = rbind(outdata, outdata_updated)
  }
  
  cat(paste0("Total signals: ", nrow(outdata_updated), "\n"))
  
  return(outdata_updated)
}


