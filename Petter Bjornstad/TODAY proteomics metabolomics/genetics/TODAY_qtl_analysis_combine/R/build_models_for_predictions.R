library(caret)
library(pROC)
library(PRROC)
library(randomForest)
library(gbm)
library(e1071)
library(nnet)
library(xgboost)
library(doParallel)
library(foreach)
library(glmnet)
library(ranger)
library(C50)
library(naivebayes)
library(class)

# Function to build and evaluate multiple classification models
build_classification_models <- function(data, outcome_vector, sample_ids = NULL,
                                        models = c("glm", "rf", "gbm", "xgbTree", "svmRadial", "nnet", 
                                                   "glmnet", "ranger", "C5.0", "nb", "knn", "LogitBoost"),
                                        train_prop = 0.8,
                                        cv_folds = 10,
                                        cv_repeats = 3,
                                        n_cores = NULL,
                                        max_memory_gb = NULL,
                                        parallel_models = TRUE,
                                        seed = 123) {
  
  # Set seed for reproducibility
  set.seed(seed)
  
  # Setup parallel processing
  if (is.null(n_cores)) {
    n_cores <- min(detectCores() - 1, length(models))  # Reserve 1 core, limit to number of models
  }
  
  if (n_cores > 1) {
    cat("Setting up parallel processing with", n_cores, "cores...\n")
    
    # Register parallel backend
    cl <- makeCluster(n_cores)
    registerDoParallel(cl)
    
    # Memory management for cluster workers
    if (!is.null(max_memory_gb)) {
      memory_per_core <- floor(max_memory_gb * 1024 / n_cores)  # MB per core
      cat("Allocating approximately", memory_per_core, "MB per core\n")
    }
    
    # Ensure cluster cleanup
    on.exit({
      cat("Cleaning up parallel cluster...\n")
      try(stopCluster(cl), silent = TRUE)
      registerDoSEQ()  # Reset to sequential
    }, add = TRUE)
    
  } else {
    cat("Using sequential processing\n")
  }
  
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
  
  # Set up cross-validation with parallel processing
  ctrl <- trainControl(
    method = "repeatedcv",
    number = cv_folds,
    repeats = cv_repeats,
    classProbs = TRUE,
    summaryFunction = twoClassSummary,
    savePredictions = "final",
    allowParallel = (n_cores > 1),  # Enable parallel processing for CV
    verboseIter = FALSE,
    returnData = FALSE  # Save memory by not returning training data
  )
  
  # Initialize results list
  results_list <- list()
  model_list <- list()
  
  # Train models - parallel or sequential based on settings
  if (parallel_models && n_cores > 1 && length(models) > 1) {
    
    cat("Training", length(models), "models in parallel...\n")
    
    # Export necessary objects to cluster
    clusterExport(cl, c("train_data", "train_outcome", "ctrl", "test_data", "test_outcome", 
                        "train_ids", "test_ids"), envir = environment())
    
    # Load required libraries on each worker
    clusterEvalQ(cl, {
      library(caret)
      library(pROC)
      library(PRROC)
      library(randomForest)
      library(gbm)
      library(e1071)
      library(nnet)
      library(xgboost)
      library(glmnet)
      library(ranger)
      library(C50)
      library(naivebayes)
      library(class)
    })
    
    # Export the helper function to cluster
    clusterExport(cl, "train_single_model", envir = environment())
    
    # Train models in parallel with better error handling
    model_results <- foreach(model_name = models, 
                             .packages = c("caret", "pROC", "PRROC"),
                             .errorhandling = "pass") %dopar% {
                               
                               # Capture detailed errors
                               tryCatch({
                                 train_single_model(model_name, train_data, train_outcome, test_data, test_outcome, 
                                                    train_ids, test_ids, ctrl)
                               }, error = function(e) {
                                 list(error = paste("Model", model_name, "failed:", e$message),
                                      model_name = model_name)
                               })
                             }
    
    # Process results with detailed error reporting
    for (i in seq_along(models)) {
      model_name <- models[i]
      result <- model_results[[i]]
      
      if (inherits(result, "try-error")) {
        cat("Error training", model_name, ":", as.character(result), "\n")
      } else if (!is.null(result$error)) {
        cat(result$error, "\n")
      } else if (is.null(result)) {
        cat("Error training", model_name, ": returned NULL result\n")
      } else {
        results_list[[paste0(model_name, "_train")]] <- result$train_results
        results_list[[paste0(model_name, "_test")]] <- result$test_results
        model_list[[model_name]] <- result$model
        
        if (!exists("predictions_list")) predictions_list <- list()
        predictions_list[[model_name]] <- result$predictions
        cat("Successfully trained", model_name, "\n")
      }
    }
    
    # If parallel training failed, fall back to sequential
    if (length(model_list) == 0) {
      cat("Parallel training failed for all models. Falling back to sequential training...\n")
      
      for (model_name in models) {
        cat("Training", model_name, "sequentially...\n")
        
        result <- train_single_model(model_name, train_data, train_outcome, test_data, test_outcome, 
                                     train_ids, test_ids, ctrl)
        
        if (!is.null(result)) {
          results_list[[paste0(model_name, "_train")]] <- result$train_results
          results_list[[paste0(model_name, "_test")]] <- result$test_results
          model_list[[model_name]] <- result$model
          
          if (!exists("predictions_list")) predictions_list <- list()
          predictions_list[[model_name]] <- result$predictions
          cat("Successfully trained", model_name, "\n")
        }
      }
    }
    
  } else {
    # Sequential training (original approach)
    cat("Training models sequentially...\n")
    
    for (model_name in models) {
      cat("Training", model_name, "...\n")
      
      result <- train_single_model(model_name, train_data, train_outcome, test_data, test_outcome, 
                                   train_ids, test_ids, ctrl)
      
      if (!is.null(result)) {
        results_list[[paste0(model_name, "_train")]] <- result$train_results
        results_list[[paste0(model_name, "_test")]] <- result$test_results
        model_list[[model_name]] <- result$model
        
        if (!exists("predictions_list")) predictions_list <- list()
        predictions_list[[model_name]] <- result$predictions
        cat("Successfully trained", model_name, "\n")
      }
    }
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
                               train_prop = 0.8, cv_folds = 10, cv_repeats = 3, 
                               n_cores = NULL, max_memory_gb = NULL, seed = 123) {
  
  result <- build_classification_models(
    data = data,
    outcome_vector = outcome_vector,
    sample_ids = sample_ids,
    models = model_type,
    train_prop = train_prop,
    cv_folds = cv_folds,
    cv_repeats = cv_repeats,
    n_cores = n_cores,
    max_memory_gb = max_memory_gb,
    parallel_models = FALSE,  # Single model doesn't need model-level parallelization
    seed = seed
  )
  
  return(result)
}

# Helper function to train a single model (used for both parallel and sequential)
train_single_model <- function(model_name, train_data, train_outcome, test_data, test_outcome, 
                               train_ids, test_ids, ctrl) {
  
  tryCatch({
    # Train model with specific configurations
    if (model_name == "nnet") {
      model <- train(
        x = train_data, y = train_outcome, method = model_name,
        trControl = ctrl, metric = "ROC", trace = FALSE, linout = FALSE
      )
    } else if (model_name == "xgbTree") {
      model <- train(
        x = train_data, y = train_outcome, method = model_name,
        trControl = ctrl, metric = "ROC", verbosity = 0
      )
    } else if (model_name == "rf") {
      # Optimize Random Forest for high-performance computing
      model <- train(
        x = train_data, y = train_outcome, method = model_name,
        trControl = ctrl, metric = "ROC",
        ntree = 1000,  # More trees for better performance on clusters
        importance = TRUE
      )
    } else {
      model <- train(
        x = train_data, y = train_outcome, method = model_name,
        trControl = ctrl, metric = "ROC"
      )
    }
    
    # Get CV results
    cv_results <- model$results
    best_tune_idx <- which.max(cv_results$ROC)
    cv_roc <- cv_results$ROC[best_tune_idx]
    cv_sens <- cv_results$Sens[best_tune_idx]
    cv_spec <- cv_results$Spec[best_tune_idx]
    
    # Make predictions
    test_pred_probs <- predict(model, test_data, type = "prob")
    test_pred_class <- predict(model, test_data, type = "raw")
    train_pred_probs <- predict(model, train_data, type = "prob")
    train_pred_class <- predict(model, train_data, type = "raw")
    
    # Store predictions
    train_predictions <- data.frame(
      sample_id = if (!is.null(train_ids)) train_ids else 1:length(train_outcome),
      model = model_name, actual = train_outcome, prob_positive = train_pred_probs[, 2],
      predicted_class = train_pred_class, dataset = "training", stringsAsFactors = FALSE
    )
    
    test_predictions <- data.frame(
      sample_id = if (!is.null(test_ids)) test_ids else 1:length(test_outcome),
      model = model_name, actual = test_outcome, prob_positive = test_pred_probs[, 2],
      predicted_class = test_pred_class, dataset = "test", stringsAsFactors = FALSE
    )
    
    # Calculate metrics
    roc_obj <- roc(test_outcome, test_pred_probs[, 2], quiet = TRUE)
    optimal_threshold <- coords(roc_obj, "best", ret = "threshold", best.method = "youden")
    cm_test <- confusionMatrix(test_pred_class, test_outcome, positive = levels(test_outcome)[2])
    
    test_TP <- cm_test$table[2, 2]; test_TN <- cm_test$table[1, 1]
    test_FP <- cm_test$table[1, 2]; test_FN <- cm_test$table[2, 1]
    test_roc_auc <- auc(roc_obj)
    
    test_outcome_numeric <- as.numeric(test_outcome == levels(test_outcome)[2])
    pr_obj <- pr.curve(scores.class0 = test_pred_probs[, 2], 
                       weights.class0 = test_outcome_numeric, curve = FALSE)
    test_aurpc <- pr_obj$auc.integral
    
    # Create results dataframes
    train_results <- data.frame(
      Model = model_name, Dataset = "Training_CV", TP = NA, TN = NA, FP = NA, FN = NA,
      ROC_AUC = as.numeric(cv_roc), Sensitivity = as.numeric(cv_sens), 
      Specificity = as.numeric(cv_spec), Precision = NA, Recall = as.numeric(cv_sens),
      AURPC = NA, F1 = NA, PPV = NA, NPV = NA, Threshold = as.numeric(optimal_threshold),
      stringsAsFactors = FALSE
    )
    
    test_results <- data.frame(
      Model = model_name, Dataset = "Test", TP = test_TP, TN = test_TN, FP = test_FP, FN = test_FN,
      ROC_AUC = as.numeric(test_roc_auc), Sensitivity = as.numeric(cm_test$byClass["Sensitivity"]),
      Specificity = as.numeric(cm_test$byClass["Specificity"]), Precision = as.numeric(cm_test$byClass["Precision"]),
      Recall = as.numeric(cm_test$byClass["Sensitivity"]), AURPC = as.numeric(test_aurpc),
      F1 = as.numeric(cm_test$byClass["F1"]), PPV = as.numeric(cm_test$byClass["Pos Pred Value"]),
      NPV = as.numeric(cm_test$byClass["Neg Pred Value"]), Threshold = as.numeric(optimal_threshold),
      stringsAsFactors = FALSE
    )
    
    return(list(
      model = model,
      train_results = train_results,
      test_results = test_results,
      predictions = rbind(train_predictions, test_predictions)
    ))
    
  }, error = function(e) {
    cat("Error training", model_name, ":", e$message, "\n")
    return(NULL)
  })
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
# # High-performance cluster usage examples:
# 
# # Option 1: Sequential processing (most reliable)
# results <- build_classification_models(
#   data = my_data,
#   outcome_vector = disease_status,
#   sample_ids = sample_ids,
#   models = c("glm", "rf", "gbm", "xgbTree"),
#   n_cores = 32,              # Still uses parallel CV within each model
#   parallel_models = FALSE    # Train models sequentially
# )
# 
# # Option 2: Try parallel first, fallback to sequential
# results <- build_classification_models(
#   data = my_data,
#   outcome_vector = disease_status,
#   sample_ids = sample_ids,
#   models = c("glm", "rf", "xgbTree"),
#   n_cores = 8,
#   max_memory_gb = 32,
#   parallel_models = TRUE     # Will fallback to sequential if parallel fails
# )
# 
# # Option 3: Conservative cluster settings
# results <- build_classification_models(
#   data = my_data,
#   outcome_vector = disease_status,
#   sample_ids = sample_ids,
#   models = c("rf", "xgbTree"),
#   n_cores = 4,               # Fewer cores
#   cv_folds = 3,              # Fewer CV folds
#   cv_repeats = 1,            # Single repeat
#   parallel_models = FALSE
# )
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