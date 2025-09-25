source(here("R/build_models_for_predictions.R"))


load_model_data = function(ii, models, my_y)
{
  infile       = models[ii, "file_model"]
  my_phenotype = models[ii, "phenotype" ]
  my_id        = models[ii, "ID"        ]
  my_model     = models[ii, "model"     ]
  indata       = readRDS((infile))

  if(my_model == "m1"){my_vars = c("genotype_id", my_id  )}
  if(my_model == "m2"){my_vars = c("genotype_id", "value_norm")}
  if(my_model == "m3"){my_vars = c("genotype_id", "value_norm", my_id)}
  if(my_model == "m4"){my_vars = c("genotype_id", "value_norm", my_id)}
  
  if(ii == 1)
  {
    indata[["data"]]       = indata[["data"]]
    indata[["data"]][,"y"] = indata[["data"]][,my_y]
    my_vars                = c(my_vars, c("y", "sex", "age_norm", paste0("PC", 1:10)))
  }

  to_model = indata[["data"]][, my_vars]
  colnames(to_model) = sub("value_norm", my_phenotype, colnames(to_model))

  return(to_model)
}

run_predictions = function(my_y, best_models_filtered)
{
  models  = best_models_filtered %>% filter(outcome == my_y)
  my_data = as.data.frame(do.call("cbind", lapply(1:nrow(models), function(ii)
  {
    load_model_data(ii, models, my_y)
  })), stringsAsFactors = FALSE)

  my_data        = my_data[ , !duplicated(names(my_data)) ]   
  disease_status = my_data$y
  genotype_ids   = my_data$genotype_id
  my_data        = my_data %>% select(-genotype_id, -y)
  
  # Run multiple models
  results <- suppressWarnings(build_classification_models(
    data            = my_data,
    outcome_vector  = disease_status,
    sample_ids      = genotype_ids,
    models          = c("glm", "rf", "gbm", "xgbTree", "svmRadial", "nnet",  "glmnet", "ranger", "C5.0", "knn", "LogitBoost"),
    #models          = c("glm"),
    n_cores         = 32,      
    parallel_models = TRUE,
    max_memory_gb   = 1000,   
    cv_folds        = 10,
    cv_repeats      = 2
    ))
  
  saveRDS(results, here("output/predictions", paste(my_y, "rds", sep = ".")))
  
  return(results$results %>% mutate(y = my_y))
}