
load_genotypes = function()
{
  # cbftools query to get all the phenotypes
}

clean_genotypes = function(gts)
{
  rownames(gts) = gts[,1]
  ids           = gts[,1]
  gts[,1] = NULL
  
  colnames(gts) = gsub(".*\\]([^:]+):.*", "\\1", colnames(gts))

  gts = gts %>%
    mutate(across(everything(), ~ recode(., 
      "0|0" = 0,
      "0|1" = 1,
      "1|0" = 1,
      "1|1" = 2 
    )))
  
  if(nrow(gts) == 1)
  {
    id  = rownames(gts)[[1]]
    gts = data.frame(genotype_id = colnames(gts),
                     genotype    = as.character(gts[1,])
                     )
    
    colnames(gts)[[2]] = id
  }else
  {
    gts = as.data.frame((t(gts)))
    gts$genotype_id = rownames(gts)
  }
  
  
  return(gts[,c("genotype_id", ids)])
}

get_genotypes = function(sumstats, tmp_folder, genotypes_prefix = here("data/genotypes"))
{
  ids_file   = paste(tmp_folder, "ids.txt", sep = "/")
  out_prefix = paste(tmp_folder, "genotypes", sep = "/")
  vcf        = paste(out_prefix, "vcf.gz", sep = ".")
  
  writeLines((sumstats %>% arrange(CHROM, POS))[,"ID"], ids_file)
  
  # Extract variants
  command = paste("plink2",
                  "--silent" ,
                  "--pfile"  , genotypes_prefix,
                  "--extract", ids_file,
                  "--export" , "vcf", "bgz",
                  "--out"    , out_prefix,
                  "")
  
  system(command)
  
  #bcftools = "/mmfs1/gscratch/togo/matteo/software/conda/r_deps/bin/bcftools"
  bcftools = "bcftools"
  
  # index VCF
  command = paste(bcftools, "index", "-t", vcf)
  
  #system(command)
  
  # Extract variants DF
  command = paste(bcftools, "query", "-H", vcf,
                  "-f", "%ID[,%GT]\n",
                  "")
  
  indata = clean_genotypes(fread(cmd = command, sep = ",", header =TRUE, data.table = FALSE))
  
  return(indata)
}

run_models = function(jj, gt, outcome, totest, covariates_totest)
{
  totest = totest %>%
    rename(gt = all_of(gt),
           y  = all_of(outcome)
           ) %>%
    select(all_of(c("y", "gt", "value", covariates_totest))) %>%
    filter(is.na(y) == FALSE & is.na(gt) == FALSE & is.na(value) == FALSE) %>%
    mutate(y = ifelse(y > 0, yes = 1, no = 0))
  
  null_model = glm(y ~ 1, family = binomial, data = totest)

  # Define models to test
  formulas <- list(
    m1 = as.formula(paste("y", paste(c("gt"                     , covariates_totest), collapse = " + "), sep = " ~ ")),
    m2 = as.formula(paste("y", paste(c(      "value"            , covariates_totest), collapse = " + "), sep = " ~ ")), #y ~ value,
    m3 = as.formula(paste("y", paste(c("gt", "value"            , covariates_totest), collapse = " + "), sep = " ~ ")), #y ~ gt + value,
    m4 = as.formula(paste("y", paste(c("gt", "value", "gt:value", covariates_totest), collapse = " + "), sep = " ~ ")) #y ~ gt + value + gt:value
  )
  
  if(jj > 1){formulas = formulas[c("m1", "m3", "m4")]}
  
  df       = totest
  df$y     = factor(df$y, levels = c(0,1), labels = c("No","Yes"))
  df$gt    = as.numeric(df$gt)

  train_control <- trainControl(
    method = "cv", number = 3,
    classProbs = TRUE,
    summaryFunction = twoClassSummary
  )
  
  # Loop through models
  results <- lapply(names(formulas), function(name) {
    
    f <- formulas[[name]]

    # --- 1. Predictive performance (AUC via caret) ---
    cv_fit <- train(
      f, data = df,
      method = "glm", family = "binomial",
      trControl = train_control,
      metric = "ROC"
    )
    auc <- max(cv_fit$results$ROC)
    
    # --- 2. Coefficients from glm ---
    glm_fit    <- glm(f, data = df, family = binomial)
    coef_table <- tidy(glm_fit)  # broom::tidy gives beta, SE, p-value
    
    list(
      model   = name,
      #formula = deparse(f),
      auc     = auc,
      coef    = as.data.frame(coef_table) %>% filter(term %in% c("gt", "value", "gt:value")),
      aic     = AIC(glm_fit),
      bic     = BIC(glm_fit),
      lr_test = anova(null_model, glm_fit, test = "Chisq")
    )
  })
  
  summary_table <- as.data.frame(rbindlist(lapply(results, function(r) {
    data.frame(
      ID        = gt,
      outcome   = outcome,
      model     = r$model,
      #formula   = r$formula,
      term      = r$coef$term,
      estimate  = r$coef$estimate,
      std_error = r$coef$std.error,
      p_value   = r$coef$p.value,
      auc       = r$auc,
      aic       = r$aic,
      bic       = r$bic,
      lr_pvalue = r$lr_test$`Pr(>Chi)`[[2]]
    )
  })), stringsAsFactors = FALSE)
  
  return(summary_table)
}

find_variants = function(ii, phenotypes, sumstats_full, id_conversion, covariates_to_qtl, phenotype_data, phenotpye_outcomes_to_qtl, tmp_folder)
{
  my_phenotype = phenotypes[ii, "phenotype"]
  my_analysis  = phenotypes[ii, "analysis" ]
  sumstats     = sumstats_full %>% 
    filter(analysis == my_analysis & phenotype == my_phenotype)

  outfile = here("output/models", paste(my_analysis, my_phenotype, "rds", sep = "."))
  
  if(file.exists(outfile) == FALSE)
  {
    message(paste(my_analysis, my_phenotype))
    
    gts = get_genotypes(sumstats, tmp_folder, genotypes_prefix = here("data/genotypes"))

    phenotype_values = phenotype_data[[my_analysis]] %>% 
      select(all_of(c("sample_id", my_phenotype))) %>%
      setNames(c("phenotype_id", "value"))
    
    # combine all data to a single data.frame
    totest = merge(gts, id_conversion)
    totest = merge(totest, covariates_to_qtl, by.x = "genotype_id" , by.y = "IID"      )
    totest = merge(totest, phenotype_values )
    totest = merge(totest, phenotpye_outcomes_to_qtl, by.x = "genotype_id", by.y = "IID" )
    totest = totest %>%
      mutate(age_norm   = (age   - mean(age  , na.rm = TRUE)) / sd(age  , na.rm = TRUE),
             value_norm = (value - mean(value, na.rm = TRUE)) / sd(value, na.rm = TRUE),
             sex        = sex   - 1,
             race1      = race1 - 1,
             race2      = race2 - 1,
             race3      = race3 - 1
      )
    
    vars_totest       = sumstats$ID
    outcomes_totest   = setdiff(colnames(phenotpye_outcomes_to_qtl), c("IID"))
    covariates_totest = c("age_norm", setdiff(colnames(covariates_to_qtl), c("IID", "age")))
    var2outcome       = as.data.frame(expand_grid(gt = vars_totest, y = outcomes_totest))
    results           = as.data.frame(rbindlist(lapply(1:nrow(var2outcome), function(jj)
    {
      gt = var2outcome[jj, "gt"]
      y  = var2outcome[jj, "y" ]
      
      #message(paste(my_phenotype, gt, y))
      
      run_models(jj, gt, y, totest, covariates_totest)
    })), stringsAsFactors = FALSE)
    
    best_model = results %>% select(ID, outcome, model, auc, aic, bic, lr_pvalue) %>%
      group_by(ID, outcome) %>%
      slice(which.min(bic)) %>%
      ungroup()
    
    out = list(data       = totest,
               results    = as.data.frame(results),
               best_model = as.data.frame(best_model),
               variants   = vars_totest
               )
    
    saveRDS(out, outfile)
  }
}
