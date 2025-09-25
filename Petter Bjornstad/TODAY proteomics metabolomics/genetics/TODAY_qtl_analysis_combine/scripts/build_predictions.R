################################################################################
# Take all the models, extract all the one significant for one outcome;
# Get all the elements associated with it, and build a prediction model

source("R/functions.R")

best_models = fread(here("output/best_models.txt"), sep = "\t", header = TRUE, data.table = FALSE)

best_models_filtered = best_models %>% filter(p_value <= 0.05 & lr_pvalue <= 0.05)
outcomes_totest      = sort(unique(best_models_filtered$outcome))

################################################################################
# Run predictions

dir.create(here("output/predictions"), showWarnings = FALSE)

source(here("R/run_predictions.R"))

predictions = as.data.frame(rbindlist(lapply(outcomes_totest, function(my_y)
{
  message(my_y)
  
  run_predictions(my_y, best_models_filtered)
})), stringsAsFactors = FALSE)

fwrite(predictions, here("output", "predictions.txt"), sep = "\t", col.names = TRUE, row.names = FALSE)
