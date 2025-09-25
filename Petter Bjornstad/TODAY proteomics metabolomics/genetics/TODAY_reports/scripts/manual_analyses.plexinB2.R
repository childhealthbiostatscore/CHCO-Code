suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(here))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(dplyr))

suppressPackageStartupMessages(library(reticulate    ))

my_lib_path = "/mmfs1/gscratch/togo/matteo/software/R"
.libPaths(my_lib_path)

conda_env_path = "/mmfs1/gscratch/togo/matteo/software/conda/r_deps"
# tell the compiler where to find headers/libraries
Sys.setenv(
  #PYTHONHOME = "/gscratch/scrubbed/mdanto/envs/r_deps/bin/python", # OLD, maybe useless
  CPATH           = paste(paste(conda_env_path, "include", sep = "/"), Sys.getenv("CPATH"          ), sep = ":"),
  LIBRARY_PATH    = paste(paste(conda_env_path, "lib"    , sep = "/"), Sys.getenv("LIBRARY_PATH"   ), sep = ":"),
  LD_LIBRARY_PATH = paste(paste(conda_env_path, "lib"    , sep = "/"), Sys.getenv("LD_LIBRARY_PATH"), sep = ":"),
  PATH            = paste(paste(conda_env_path, "bin"    , sep = "/"), Sys.getenv("PATH"           ), sep = ":")
)

################################################################################
# Get important information

top_snps = fread(here("data/top_snps_description.txt"), sep = "\t", header = TRUE, data.table = FALSE)
indata   = as.data.frame(top_snps %>% filter(TargetFullName == "Plexin-B2") %>%
                           slice(which.min(P))
)

sumstats   = readRDS(indata$sumstats_file)

vcf_file   = paste("/mmfs1/gscratch/togo/matteo/projects/TODAY_genetics/data/genotyping/imputation", paste0("chr", indata$CHROM, "_imputed.vcf.gz"), sep = "/")
my_command = paste("bcftools", "query",
                   "-f", "'%ID[\\t%GT]\\n'",
                   "-r", paste0("chr", indata$CHROM, ":", indata$POS - 1, "-", indata$POS + 1),
                   "-H",
                   vcf_file)

id_conversion = fread("/mmfs1/gscratch/togo/matteo/projects/TODAY_qtl/data/to_qtl/id_conversion.txt", sep = "\t", header = TRUE, data.table = FALSE)

gts           = fread(cmd = my_command, sep = "\t", header = TRUE, data.table = FALSE) %>% select(-1)
colnames(gts) = cleaned <- sub(".*\\](.*?)\\:.*", "\\1", colnames(gts))
gts           = data.frame(genotype_id = colnames(gts), genotype = as.character(as.matrix(gts[1, ]))) %>%
  mutate(genotype = case_when(
    genotype == "0|0" ~ 0,
    genotype %in% c("0|1", "1|0") ~ 0.5,
    genotype == "1|1" ~ 1
  ))
phenotypes    = fread("/mmfs1/gscratch/togo/matteo/projects/TODAY_qtl/data/to_qtl/phenotypes_proteomics.txt", sep = "\t", header = TRUE, data.table = FALSE, select = c("sample_id", indata$phenotype)) %>%
  set_names(c("phenotype_id", "phenotype"))

################################################################################
# Manhattan plot

manhattan_plot <- function(df, genomewideline = 5e-8, max_nonsig = 1e5) {
  # order by chromosome and position
  df <- df[order(df$CHROM, df$POS), ]
  
  # separate significant and non-significant
  sig <- df$P < genomewideline
  nonsig <- !sig
  
  # subsample non-significant variants if needed
  if (sum(nonsig) > max_nonsig) {
    set.seed(1)  # reproducibility
    keep_nonsig <- sample(which(nonsig), max_nonsig)
    df <- df[c(which(sig), keep_nonsig), ]
  }
  
  # recompute order
  df <- df[order(df$CHROM, df$POS), ]
  
  # build cumulative positions
  chr_lengths <- tapply(df$POS, df$CHROM, max)
  chr_offsets <- c(0, cumsum(as.numeric(chr_lengths))[-length(chr_lengths)])
  names(chr_offsets) <- names(chr_lengths)
  df$cum_pos <- df$POS + chr_offsets[as.character(df$CHROM)]
  
  # transform P to -log10
  df$logp <- -log10(df$P)
  
  # colors: alternate chromosomes
  colvec <- rep(c("grey20", "grey60"), length.out = length(unique(df$CHROM)))
  df$col <- colvec[as.numeric(as.factor(df$CHROM))]
  
  # base plot
  plot(df$cum_pos, df$logp, pch = 20, col = df$col,
       xlab = "Chromosome", ylab = "-log10(P)", xaxt = "n")
  
  # chromosome labels in the middle
  axis(1, at = tapply(df$cum_pos, df$CHROM, mean), labels = names(chr_lengths))
  
  # significance threshold line
  abline(h = -log10(genomewideline), col = "blue", lty = 2)
  abline(h = 0, col = "black", lty = 1)
  
  # highlight significant variants in red
  sig <- df$P < genomewideline
  points(df$cum_pos[sig], df$logp[sig], col = "red", pch = 20)
  
  # label the most significant variant
  top <- df[which.min(df$P), ]
  text(top$cum_pos, top$logp, labels = top$ID, pos = 2, col = "#000000")
}

manhattan_plot(sumstats)

manhattan_plot(sumstats %>% filter(CHROM == indata$CHROM & POS >= indata$START - 1e4 & POS <= indata$END + 1e4))

################################################################################
# Boxplot genotypes

race_data           = fread("/mmfs1/gscratch/togo/matteo/projects/TODAY_qtl/data/input/PAT.csv", sep = ",", header = TRUE, data.table = FALSE, select = c("releaseid", "race"))
colnames(race_data) = c("phenotype_id", "race_n")
race2name           = data.frame(race_n = 1:4, race = c("White", "Black", "Hispanic", "Other"), color = c("#00bb00", "#bb0000", "#bbbb00", "#00bbbb"))
race_data           = merge(race_data, race2name)

to_boxplot = merge(id_conversion, gts)
to_boxplot = merge(to_boxplot   , phenotypes)
to_boxplot = merge(to_boxplot   , race_data)

boxplot(phenotype ~ genotype, data = to_boxplot, outline = FALSE, 
        axes = FALSE, ylab = indata$TargetFullName)

axis(2)
axis(1, at = 1:3, 
     labels = c(paste(indata$REF, indata$REF, sep = "/"),
                paste(indata$REF, indata$ALT, sep = "/"),
                paste(indata$ALT, indata$ALT, sep = "/")
                ))


table(to_boxplot$genotype)
nrow(to_boxplot)


boxplot(phenotype ~ race_n + genotype, data = to_boxplot, outline = FALSE, 
        axes = FALSE, ylab = indata$TargetFullName, col = race2name$color)

legend("topright", legend = race2name$race, col = race2name$color, pch = 16, pt.cex = 2.5)

axis(2)
axis(1, at = (1:3) * 4, 
     labels = c(paste(indata$REF, indata$REF, sep = "/"),
                paste(indata$REF, indata$ALT, sep = "/"),
                paste(indata$ALT, indata$ALT, sep = "/")
     ))

################################################################################
# Add outcome to the analysis

primout = fread("/mmfs1/gscratch/togo/matteo/projects/TODAY_qtl/data/input/PRIMOUT.csv", sep = ",", header = TRUE, data.table = FALSE) %>%
  rename(phenotype_id = "releaseid")

to_model = merge(to_boxplot, primout)

# A model where outcome is a function of the genotype
summary(glm(outcome ~ genotype + tx, data = to_model, family = "binomial"))

# A model where outcome is a function of the phenotype
summary(glm(outcome ~ phenotype + tx, data = to_model, family = "binomial"))

# A model where outcome is a function of both
summary(glm(outcome ~ phenotype + genotype + tx, data = to_model, family = "binomial"))

boxplot(phenotype ~ outcome + genotype, data = to_model, outline = FALSE, 
        axes = FALSE, ylab = indata$TargetFullName, col = c("#bbbbbb", "#ff0000"))

legend("topright", legend = paste("outcome", 0:1, sep = " = "), col = c("#bbbbbb", "#ff0000"), pch = 16, pt.cex = 2.5)

axis(2)
axis(1, at = (1:3) * 2, 
     labels = c(paste(indata$REF, indata$REF, sep = "/"),
                paste(indata$REF, indata$ALT, sep = "/"),
                paste(indata$ALT, indata$ALT, sep = "/")
     ))

table(to_model %>% select(genotype, outcome))

to_forest = as.data.frame(rbindlist(lapply((0:2)/2, function(gt)
{
  x = to_model %>% filter(genotype == gt)
  
  message(paste(gt, nrow(x)))
  test = t.test((x %>% filter(outcome == 1))[,"phenotype"],
                (x %>% filter(outcome == 0))[,"phenotype"],
                )
  
  print(str(test))
  
  return(data.frame(genotype = gt,
                    mean1    = test$estimate[[1]],
                    mean2    = test$estimate[[2]],
                    pval     = test$p.value,
                    ci1      = test$conf.int[[1]],
                    ci2      = test$conf.int[[2]]
  ) %>% mutate(delta = mean1 - mean2))
  
})), stringsAsFactors = FALSE)

plot(1,1, type = "n", 
     xlim = range(c(to_forest$ci1, to_forest$ci2)), ylim = c(0.5, 3.5), 
     xlab = "Delta Plexin-B2 (outcome == 1 - outcome == 0)", ylab = "Genotype", axes = FALSE)

segments(x0 = to_forest$ci1, x1 = to_forest$ci2, y0 = 1:3, lwd = 3, col = "#bbbbbb")
points  (x  = to_forest$delta, y = 1:3, pch = 21, cex = 3, bg = "#ffff00")
abline  (v  = 0, lty = "dashed", lwd = 2)

axis(1)
axis(2, at = (1:3), 
     labels = c(paste(indata$REF, indata$REF, sep = "/"),
                paste(indata$REF, indata$ALT, sep = "/"),
                paste(indata$ALT, indata$ALT, sep = "/")
     ))

# Cox PH

suppressPackageStartupMessages(library(survival))
suppressPackageStartupMessages(library(survminer))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(forestplot))
suppressPackageStartupMessages(library(patchwork))

data = to_model

phenotype_quartiles  <- quantile(data$phenotype, probs = c(0.25, 0.75), na.rm = TRUE)

data$phenotype_group <- cut(data$phenotype, 
                            breaks = c(-Inf, phenotype_quartiles[1], phenotype_quartiles[2], Inf),
                            labels = c("Bottom_quartile", "Middle_50pct", "Top_quartile"),
                            include.lowest = TRUE)

# Create survival object
data$days_to_outcome = rowSums(data %>% select(daystocensor, daystopo_s), na.rm = TRUE)
data$phenotype_group <- relevel(data$phenotype_group, ref = "Middle_50pct")

surv_object <- Surv(time = data$days_to_outcome, event = data$outcome)

cox_model <- coxph(surv_object ~ phenotype, data = data)
summary(cox_model)

cox_model <- coxph(surv_object ~ phenotype + genotype, data = data)
summary(cox_model)

# Fit Cox proportional hazards model
cox_model <- coxph(surv_object ~ genotype + phenotype_group, data = data)

# View results
summary(cox_model)


# 1. Kaplan-Meier curves by genotype
fit_genotype <- survfit(surv_object ~ genotype, data = data)

p1 <- ggsurvplot(fit_genotype, 
                 data = data,
                 title = "Survival by Genotype",
                 xlab = "Days to Outcome",
                 ylab = "Survival Probability",
                 pval = TRUE,
                 conf.int = TRUE,
                 risk.table = TRUE,
                 ggtheme = theme_minimal())

print(p1)

# 2. Kaplan-Meier curves by phenotype group
fit_phenotype <- survfit(surv_object ~ phenotype_group, data = data)

p2 <- ggsurvplot(fit_phenotype, 
                 data = data,
                 title = "Survival by Phenotype Group",
                 xlab = "Days to Outcome",
                 ylab = "Survival Probability",
                 pval = TRUE,
                 conf.int = TRUE,
                 risk.table = TRUE,
                 ggtheme = theme_minimal())

print(p2)

# 4. Phenotype for reference allele
toplot = data %>% filter(genotype == 0)
toplot$combined_group <- interaction(toplot$genotype, toplot$phenotype_group, sep = "_")
surv_object_small <- Surv(time = toplot$days_to_outcome, event = toplot$outcome)
fit_combined <- survfit(surv_object_small ~ combined_group, data = toplot)

p4 <- ggsurvplot(fit_combined, 
                 data = toplot,
                 title = "Survival in individuals with homozygous reference genotype",
                 xlab = "Days to Outcome",
                 ylab = "Survival Probability",
                 pval = TRUE,
                 conf.int = TRUE,
                 risk.table = TRUE,
                 ggtheme = theme_minimal())

print(p4)

# 4. Phenotype for heterozyous allele
toplot = data %>% filter(genotype == 0.5)
toplot$combined_group <- interaction(toplot$genotype, toplot$phenotype_group, sep = "_")
surv_object_small <- Surv(time = toplot$days_to_outcome, event = toplot$outcome)
fit_combined <- survfit(surv_object_small ~ combined_group, data = toplot)

p5 <- ggsurvplot(fit_combined, 
                 data = toplot,
                 title = "Survival in individuals with heterozygous genotype",
                 xlab = "Days to Outcome",
                 ylab = "Survival Probability",
                 pval = TRUE,
                 conf.int = TRUE,
                 risk.table = TRUE,
                 ggtheme = theme_minimal())

print(p5)



# 4. Phenotype for homozygous alternative allele
toplot = data %>% filter(genotype == 1)
toplot$combined_group <- interaction(toplot$genotype, toplot$phenotype_group, sep = "_")
surv_object_small <- Surv(time = toplot$days_to_outcome, event = toplot$outcome)
fit_combined <- survfit(surv_object_small ~ combined_group, data = toplot)

p6 <- ggsurvplot(fit_combined, 
                 data = toplot,
                 title = "Survival in individuals with homozygous alternative genotype",
                 xlab = "Days to Outcome",
                 ylab = "Survival Probability",
                 pval = TRUE,
                 conf.int = TRUE,
                 risk.table = TRUE,
                 ggtheme = theme_minimal())

print(p6)


combined_plot <- p4$plot + p5$plot + p6$plot + 
  plot_annotation(title = "Cox Regression Analysis Results",
                  tag_levels = 'A')  # Adds A, B, C labels

print(combined_plot)
