library(caret)
library(pROC)
library(PRROC)
library(randomForest)
library(gbm)
library(e1071)
library(nnet)
library(xgboost)

# Function to build and evaluate multiple classification models
build_classification_models <- function(data, outcome_vector, sample_ids = NULL,
                                        models = c("glm", "rf", "gbm", "xgbTree", "svmRadial", "nnet"),
                                        train_prop = 0.8,
                                        cv_folds = 10,
                                        cv_repeats = 3,
                                        seed = 123) {
  
  # Set seed for reproducibility
  set.seed(seed)
  
  # Validate inputs
  if (nrow(data) != length(outcome_vector)) {
    stop("Number of rows in data must match length of outcome_vector")
  }
  
  if (!is.null(sample_ids) && length(sample_ids) != length(outcome_vector)) {
    stop("Length of sample_ids must match length of outcome_vector")
  }
  
  if (!is.factor(outcome_vector)) {
    outcome_vector <- as.factor(outcome_vector)
  }
  
  # Ensure binary outcome
  if (length(levels(outcome_vector)) != 2) {
    stop("Outcome must be binary (2 levels)")
  }
  
  # Fix factor levels to be valid R variable names for caret
  original_levels <- levels(outcome_vector)
  new_levels <- make.names(original_levels)
  
  # Only change if necessary
  if (!identical(original_levels, new_levels)) {
    cat("Converting factor levels from", paste(original_levels, collapse = ", "), 
        "to", paste(new_levels, collapse = ", "), "for caret compatibility\n")
    levels(outcome_vector) <- new_levels
  }
  
  # Ensure all data columns are numeric
  cat("Converting all input data to numeric...\n")
  
  # Convert all columns to numeric, handling factors appropriately
  for (col_name in names(data)) {
    if (is.factor(data[[col_name]])) {
      # Convert factor to numeric (assumes binary factors become 0/1)
      data[[col_name]] <- as.numeric(as.character(data[[col_name]]))
    } else if (!is.numeric(data[[col_name]])) {
      # Convert other non-numeric types to numeric
      data[[col_name]] <- as.numeric(data[[col_name]])
    }
  }
  
  # Check for any NAs introduced by conversion and warn user
  na_cols <- colnames(data)[colSums(is.na(data)) > 0]
  if (length(na_cols) > 0) {
    cat("Warning: NA values introduced in columns:", paste(na_cols, collapse = ", "), "\n")
    cat("Consider checking your data conversion.\n")
  }
  
  # Create training/testing split
  train_index <- createDataPartition(outcome_vector, p = train_prop, list = FALSE)
  
  # Training data
  train_data <- data[train_index, ]
  train_outcome <- outcome_vector[train_index]
  train_ids <- if (!is.null(sample_ids)) sample_ids[train_index] else NULL
  
  # Testing data
  test_data <- data[-train_index, ]
  test_outcome <- outcome_vector[-train_index]
  test_ids <- if (!is.null(sample_ids)) sample_ids[-train_index] else NULL
  
  # Print split information
  cat("Training set size:", length(train_outcome), "samples\n")
  cat("Test set size:", length(test_outcome), "samples\n")
  cat("Training proportion:", round(length(train_outcome) / length(outcome_vector), 3), "\n\n")
  
  # Set up cross-validation
  ctrl <- trainControl(
    method = "repeatedcv",
    number = cv_folds,
    repeats = cv_repeats,
    classProbs = TRUE,
    summaryFunction = twoClassSummary,
    savePredictions = "final",
    allowParallel = TRUE
  )
  
  # Initialize results list
  results_list <- list()
  model_list <- list()
  
  # Train models
  for (model_name in models) {
    cat("Training", model_name, "...\n")
    
    tryCatch({
      # Train model
      if (model_name == "nnet") {
        # Neural networks need special handling
        model <- train(
          x = train_data,
          y = train_outcome,
          method = model_name,
          trControl = ctrl,
          metric = "ROC",
          trace = FALSE,
          linout = FALSE
        )
      } else if (model_name == "xgbTree") {
        # XGBoost needs special handling for factor outcomes
        model <- train(
          x = train_data,
          y = train_outcome,
          method = model_name,
          trControl = ctrl,
          metric = "ROC",
          verbosity = 0  # Suppress XGBoost output
        )
      } else {
        model <- train(
          x = train_data,
          y = train_outcome,
          method = model_name,
          trControl = ctrl,
          metric = "ROC"
        )
      }
      
      # Get cross-validation results from training
      cv_results <- model$results
      best_tune_idx <- which.max(cv_results$ROC)
      cv_roc <- cv_results$ROC[best_tune_idx]
      cv_sens <- cv_results$Sens[best_tune_idx]
      cv_spec <- cv_results$Spec[best_tune_idx]
      
      # Make predictions on test set
      test_pred_probs <- predict(model, test_data, type = "prob")
      test_pred_class <- predict(model, test_data, type = "raw")
      
      # Make predictions on training set for completeness
      train_pred_probs <- predict(model, train_data, type = "prob")
      train_pred_class <- predict(model, train_data, type = "raw")
      
      # Store predictions with sample IDs
      train_predictions <- data.frame(
        sample_id = if (!is.null(train_ids)) train_ids else 1:length(train_outcome),
        model = model_name,
        actual = train_outcome,
        prob_positive = train_pred_probs[, 2],
        predicted_class = train_pred_class,
        dataset = "training",
        stringsAsFactors = FALSE
      )
      
      test_predictions <- data.frame(
        sample_id = if (!is.null(test_ids)) test_ids else 1:length(test_outcome),
        model = model_name,
        actual = test_outcome,
        prob_positive = test_pred_probs[, 2],
        predicted_class = test_pred_class,
        dataset = "test",
        stringsAsFactors = FALSE
      )
      
      # Calculate optimal threshold using Youden's index on test set
      roc_obj <- roc(test_outcome, test_pred_probs[, 2], quiet = TRUE)
      optimal_threshold <- coords(roc_obj, "best", ret = "threshold", best.method = "youden")
      
      # Create confusion matrix for test set
      cm_test <- confusionMatrix(test_pred_class, test_outcome, positive = levels(test_outcome)[2])
      
      # Extract confusion matrix values for test set
      test_TP <- cm_test$table[2, 2]
      test_TN <- cm_test$table[1, 1]
      test_FP <- cm_test$table[1, 2]
      test_FN <- cm_test$table[2, 1]
      
      # Calculate test set ROC AUC
      test_roc_auc <- auc(roc_obj)
      
      # Calculate test set AURPC
      test_outcome_numeric <- as.numeric(test_outcome == levels(test_outcome)[2])
      pr_obj <- pr.curve(scores.class0 = test_pred_probs[, 2], 
                         weights.class0 = test_outcome_numeric, 
                         curve = FALSE)
      test_aurpc <- pr_obj$auc.integral
      
      # Calculate additional test metrics
      test_sensitivity <- cm_test$byClass["Sensitivity"]
      test_specificity <- cm_test$byClass["Specificity"]
      test_precision <- cm_test$byClass["Precision"]
      test_recall <- test_sensitivity
      test_f1 <- cm_test$byClass["F1"]
      test_ppv <- cm_test$byClass["Pos Pred Value"]
      test_npv <- cm_test$byClass["Neg Pred Value"]
      
      # Store training results (CV results)
      train_results <- data.frame(
        Model = model_name,
        Dataset = "Training_CV",
        TP = NA,  # Not directly available from CV
        TN = NA,
        FP = NA,
        FN = NA,
        ROC_AUC = as.numeric(cv_roc),
        Sensitivity = as.numeric(cv_sens),
        Specificity = as.numeric(cv_spec),
        Precision = NA,  # Not directly available from CV
        Recall = as.numeric(cv_sens),
        AURPC = NA,  # Not directly available from CV
        F1 = NA,  # Not directly available from CV
        PPV = NA,
        NPV = NA,
        Threshold = as.numeric(optimal_threshold),
        stringsAsFactors = FALSE
      )
      
      # Store test results
      test_results <- data.frame(
        Model = model_name,
        Dataset = "Test",
        TP = test_TP,
        TN = test_TN,
        FP = test_FP,
        FN = test_FN,
        ROC_AUC = as.numeric(test_roc_auc),
        Sensitivity = as.numeric(test_sensitivity),
        Specificity = as.numeric(test_specificity),
        Precision = as.numeric(test_precision),
        Recall = as.numeric(test_recall),
        AURPC = as.numeric(test_aurpc),
        F1 = as.numeric(test_f1),
        PPV = as.numeric(test_ppv),
        NPV = as.numeric(test_npv),
        Threshold = as.numeric(optimal_threshold),
        stringsAsFactors = FALSE
      )
      
      # Combine training and test results
      results_list[[paste0(model_name, "_train")]] <- train_results
      results_list[[paste0(model_name, "_test")]] <- test_results
      
      # Store predictions
      all_predictions <- rbind(train_predictions, test_predictions)
      if (!exists("predictions_list")) predictions_list <- list()
      predictions_list[[model_name]] <- all_predictions
      
      # Store model for later use if needed
      model_list[[model_name]] <- model
      
    }, error = function(e) {
      cat("Error training", model_name, ":", e$message, "\n")
    })
  }
  
  # Combine results and predictions
  if (length(results_list) > 0) {
    results_df <- do.call(rbind, results_list)
    rownames(results_df) <- NULL
    
    # Sort by model name and dataset for better organization
    results_df <- results_df[order(results_df$Model, results_df$Dataset), ]
    
    # Combine all predictions
    if (exists("predictions_list")) {
      all_predictions_combined <- do.call(rbind, predictions_list)
      rownames(all_predictions_combined) <- NULL
    } else {
      all_predictions_combined <- NULL
    }
    
    return(list(
      results = results_df,
      predictions = all_predictions_combined,
      models = model_list,
      train_data = train_data,
      test_data = test_data,
      train_outcome = train_outcome,
      test_outcome = test_outcome,
      train_ids = train_ids,
      test_ids = test_ids,
      split_info = list(
        total_samples = length(outcome_vector),
        train_samples = length(train_outcome),
        test_samples = length(test_outcome),
        train_proportion = round(length(train_outcome) / length(outcome_vector), 3)
      )
    ))
  } else {
    stop("No models were successfully trained")
  }
}

# Example usage function for a single model
build_single_model <- function(data, outcome_vector, sample_ids = NULL, model_type = "rf", 
                               train_prop = 0.8, cv_folds = 10, cv_repeats = 3, seed = 123) {
  
  result <- build_classification_models(
    data = data,
    outcome_vector = outcome_vector,
    sample_ids = sample_ids,
    models = model_type,
    train_prop = train_prop,
    cv_folds = cv_folds,
    cv_repeats = cv_repeats,
    seed = seed
  )
  
  return(result)
}

# Helper function to display results nicely
display_results <- function(results_object) {
  cat("Model Performance Summary:\n")
  cat("========================\n\n")
  
  results_df <- results_object$results
  
  # Group by model for cleaner display
  models <- unique(results_df$Model)
  
  for (model_name in models) {
    cat("Model:", model_name, "\n")
    model_results <- results_df[results_df$Model == model_name, ]
    
    for (i in 1:nrow(model_results)) {
      dataset <- model_results$Dataset[i]
      cat("  ", dataset, ":\n")
      cat("    ROC AUC:", round(model_results$ROC_AUC[i], 3), "\n")
      if (!is.na(model_results$F1[i])) {
        cat("    F1 Score:", round(model_results$F1[i], 3), "\n")
      }
      cat("    Sensitivity:", round(model_results$Sensitivity[i], 3), "\n")
      cat("    Specificity:", round(model_results$Specificity[i], 3), "\n")
      if (!is.na(model_results$Threshold[i])) {
        cat("    Threshold:", round(model_results$Threshold[i], 3), "\n")
      }
    }
    cat("---\n")
  }
  
  return(results_df)
}

# Helper function to extract predictions for a specific model
get_model_predictions <- function(results_object, model_name) {
  if (is.null(results_object$predictions)) {
    stop("No predictions found in results object")
  }
  
  model_preds <- results_object$predictions[
    grepl(model_name, results_object$predictions$sample_id) | 
      attr(results_object$predictions, "model") == model_name, 
  ]
  
  return(model_preds)
}

# Example usage:
# Assuming you have a data.frame called 'my_data', outcome vector 'disease_status', 
# and sample IDs 'sample_ids'
# 
# # Run multiple models with sample IDs
# results <- build_classification_models(
#   data = my_data,
#   outcome_vector = disease_status,
#   sample_ids = sample_ids,
#   models = c("glm", "rf", "gbm", "xgbTree", "svmRadial")
# )
# 
# # View results (now includes both training CV and test results)
# print(results$results)
# 
# # Access all predictions with probabilities and thresholds
# print(head(results$predictions))
# 
# # Filter predictions by dataset
# train_preds <- results$predictions[results$predictions$dataset == "training", ]
# test_preds <- results$predictions[results$predictions$dataset == "test", ]
# 
# # Access training and test sample IDs
# print(results$train_ids)
# print(results$test_ids)
# 
# # View split information
# print(results$split_info)
# 
# # Run single model
# rf_results <- build_single_model(
#   data = my_data,
#   outcome_vector = disease_status,
#   sample_ids = sample_ids,
#   model_type = "rf"
# )
# 
# # Display formatted results
# display_results(rf_results)