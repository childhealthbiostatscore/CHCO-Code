source("R/functions.R")

################################################################################
# Get all input data

sumstats_file                      = here("data/top_snps_description.txt")
id_conversion_file                 = here("data/id_conversion.txt")
phenotype_metabolomics_plasma_file = here("data/phenotypes_metabolomics_plasma.txt")
phenotype_metabolomics_urine_file  = here("data/phenotypes_metabolomics_urine.txt")
phenotype_proteomics_file          = here("data/phenotypes_proteomics.txt")
phenotype_outcomes_file            = here("data/phenotypes_outcomes.txt")
phenotype_ckm_file                 = here("data/phenotypes_ckm.txt")
covariates_file                    = here("data/to_qtl", "covariates.txt")
outcomes_file                      = here("data/to_qtl", "phenotypes.txt")

sumstats                      =  fread(sumstats_file                      , sep = "\t", header = TRUE, data.table = FALSE)
covariates                    =  fread(covariates_file                    , sep = "\t", header = TRUE, data.table = FALSE)
id_conversion                 =  fread(id_conversion_file                 , sep = "\t", header = TRUE, data.table = FALSE)
phenotype_metabolomics_plasma =  fread(phenotype_metabolomics_plasma_file , sep = "\t", header = TRUE, data.table = FALSE)
phenotype_metabolomics_urine  =  fread(phenotype_metabolomics_urine_file  , sep = "\t", header = TRUE, data.table = FALSE)
phenotype_proteomics          =  fread(phenotype_proteomics_file          , sep = "\t", header = TRUE, data.table = FALSE)
phenotype_outcomes            =  fread(outcomes_file                      , sep = "\t", header = TRUE, data.table = FALSE)

################################################################################
# Read model data

sumstats = sumstats %>%
  mutate(file_model = here("output/models", paste(analysis, phenotype, "rds", sep = ".")),
         complete   = file.exists(file_model)
         )

table(sumstats$complete)


to_model = unique(sumstats %>% 
                    filter(complete == TRUE) %>%
                    select(phenotype, TargetFullName, Target, UniProt, EntrezGeneID, EntrezGeneSymbol, analysis, file_model)
                    )

models   = as.data.frame(rbindlist(lapply(1:nrow(to_model), function(ii)
{
  infile    = to_model[ii, "file_model"]
  analysis  = to_model[ii, "analysis"  ]
  phenotype = to_model[ii, "phenotype" ]
  indata    = readRDS(infile)
  
  return(cbind(to_model[ii,], indata$results) %>% select(-a))
})), stringsAsFactors = FALSE)

best_models = models %>% 
  filter(term != "(Intercept)") %>%
  #select(analysis, phenotype, ID, outcome, model, formula, estimate, std_error, p_value, auc, aic, bic) %>%
  group_by(analysis, phenotype, outcome) %>%
  slice(which.min(bic)) %>%
  ungroup()

best_models$lr_qvalue = p.adjust(best_models$lr_pvalue, method = "bonferroni")

table(best_models$model)

fwrite(models     , here("output/models.txt"     ), sep = "\t", col.names = TRUE, row.names = FALSE)
fwrite(best_models %>% mutate(index = 1:nrow(best_models)), here("output/best_models.txt"), sep = "\t", col.names = TRUE, row.names = FALSE)

################################################################################
# Some numbers

# number of genome-wide significant associations
print(nrow(sumstats %>% filter(P <= 5e-8)))
print(nrow(unique(sumstats %>% filter(P <= 5e-8) %>% select(phenotype, analysis))))

print(table(unique(sumstats %>% filter(P <= 5e-8) %>% select(phenotype, analysis))[,"analysis"]))
print(table(unique(sumstats %>% select(phenotype, analysis))[,"analysis"]))

print(table(unique(sumstats %>% filter(P <= 5e-8) %>% select(phenotype, analysis))[,"analysis"]))

# Tested phenotypes
print(table(sort(table((unique(sumstats %>% filter(P <= 5e-8) %>% select(phenotype, ID)))[,"phenotype"]), decreasing = TRUE)))

# What variants are most represented?
print(table(sort(table((unique(sumstats %>% filter(P <= 5e-8) %>% select(phenotype, ID)))[,"ID"]), decreasing = TRUE)))

head(sort(table((unique(sumstats %>% filter(P <= 5e-8) %>% select(phenotype, ID)))[,"ID"]), decreasing = TRUE), n = 10)

# Significant associations
message(paste0("Significant associations (p < 5e-8): ", nrow(sumstats %>% filter(P < 5e-8))))
message(paste0("Significant associations (p < 1e-6): ", nrow(sumstats %>% filter(P < 1e-6))))
message(paste0("Significant associations (p < 1e-4): ", nrow(sumstats %>% filter(P < 1e-4))))

# Models
print(table((best_models %>% filter(p_value   <= 0.05))[,"model"]))
print(table((best_models %>% filter(lr_qvalue <= 0.05))[,"model"]))
print(table((best_models %>% filter(p_value <= 0.05 & lr_qvalue <= 0.05))[,"model"]))

################################################################################
# How many tested variants?
system("wc -l /mmfs1/gscratch/togo/matteo/projects/TODAY_qtl/data/to_qtl/genotypes.pvar", intern = TRUE)
################################################################################
# Manhattan plot with all the genome-wide significant variants

manhattan_plot <- function(df, sig_line = 5e-8, 
                           point_size = 1, title = "Manhattan Plot") {
  # Check columns
  if(!all(c("CHROM","POS","P") %in% names(df))) stop("Data must have CHROM, POS, P columns")
  
  # Ensure chromosome is factor
  df$CHROM <- as.factor(df$CHROM)
  
  # Sort by chromosome and position
  df <- df[order(df$CHROM, df$POS), ]
  
  # Compute cumulative positions
  df$cumPOS <- NA
  offset <- 0
  for(chr in levels(df$CHROM)) {
    df$cumPOS[df$CHROM == chr] <- df$POS[df$CHROM == chr] + offset
    offset <- max(df$cumPOS[df$CHROM == chr])
  }
  
  # Assign alternating colors
  n_chr <- length(levels(df$CHROM))

  # Plot
  layout(matrix(1))
  plot(df$cumPOS, -log10(df$P), pch = 20, col = df$color, cex = point_size,
       xlab = "Chromosome", ylab = "-log10(P)", xaxt = "n", main = title)
  
  # Add chromosome labels at the median position
  axis(1, at = tapply(df$cumPOS, df$CHROM, median), labels = levels(df$CHROM))
  
  # Add genome-wide significance line
  abline(h = -log10(sig_line), col = "red", lty = 2)
}

toplot = merge(sumstats %>% filter(P <= 5e-8), analysis2color)

manhattan_plot(toplot)

################################################################################
# Plot trait analysis:
# boxplots showing important data about models
# Manhattan plot genome-wide
# Manhattan plot focused on lead region

dir.create(here("output/plots"), showWarnings = FALSE)

source(here("R/plot_qtls.R"))

# plot CKM
to_ckm = best_models %>% filter(
  grepl("CKM", outcome) & 
    (lr_qvalue <= 0.05 | p_value <= 1e-3)
  ) %>% arrange(lr_pvalue)

invisible(lapply(1:nrow(to_ckm), function(ii)
{
  try(plot_analysis(
    as.character(to_ckm[ii, "phenotype"]), 
    as.character(to_ckm[ii, "ID"       ]), 
    as.character(to_ckm[ii, "outcome"  ]), 
    sumstats, 
    TRUE))
}))

fwrite(to_ckm %>%
         select(outcome, phenotype, Target, TargetFullName, analysis,
                ID, model, p_value, auc, aic, bic, lr_pvalue, lr_qvalue
                ) %>%
         arrange(desc(auc)),
       here("output/models_filtered.CKM.txt"), sep = "\t",
       col.names = TRUE, row.names = FALSE
       )


# Others
a = plot_analysis("seq.10442.1" , "19:55376727:C:T"   , "TGDLP"  , sumstats, TRUE)
a = plot_analysis("seq.10855.55", "22:50289633:T:C"   , "outcome", sumstats, TRUE)
a = plot_analysis("seq.19590.46", "10:79970723:TACA:T", "TGDLP"  , sumstats, TRUE)
a = plot_analysis("seq.18336.31", "19:44908684:T:C"   , "ANYDLP0", sumstats, TRUE)
a = plot_analysis("seq.3309.2"  , "1:161509955:A:G"   , "outcome", sumstats, TRUE)

a = plot_analysis("", "", "", sumstats, TRUE)




################################################################################
# Plot associations

plot_results = function(ii, best_models)
{
  id        = best_models[ii, "ID"       ]
  outcome   = best_models[ii, "outcome"  ]
  analysis  = best_models[ii, "analysis" ]
  phenotype = best_models[ii, "phenotype"]
  infile    = here("output/models", paste(analysis, phenotype, "rds", sep = "."))
  indata    = readRDS(infile)
  data      = indata[["data"]]
  data$id   = data[,id]
  data$y    = data[,outcome]
  return(indata)
}

ii = 333
a = plot_results(ii, best_models)
  
