# This is a function for performing cross validation (CV) to select an optimal model 
# using the ElasticNet. By default uses leave one out (LOO) CV, but k-fold CV
# can also be used by setting cv_method = "kfold" and folds = k. See the trainControl 
# function in caret for additional details. There are two options for the output:
# out = "min.error" produces the model with the lowest CV error, and 
# out = "1se.error" produces all "acceptable" models (CV error within 
# 1 standard error of the minimum). 

# This version was created on 7/22/25 to handle ordinal predictors as well.

easy_elasticnet <- function(data, outcome, predictors,
                            n_alphas = 10, n_lambdas = 100, max_coef = NULL,
                            model_type = "gaussian", time = NULL,
                            cv_method = "loo", folds = NULL, out = "1se.error",
                            cores = 4, seed = 3654,
                            penalty.factor = rep(1, ncol(predictors))) {
  
  require(ensr)
  require(glmnet)
  require(dplyr)
  require(ordinalNet)
  
  df <- data
  set.seed(seed)
  
  # Fix names
  colnames(df) <- make.names(colnames(df), unique = TRUE, allow_ = FALSE)
  preds <- make.names(predictors, unique = TRUE, allow_ = FALSE)
  outcome <- make.names(outcome, unique = TRUE, allow_ = FALSE)
  X <- predictors
  
  # --- Outcome handling ---
  if (model_type == "ordinal") {
    Y <- as.ordered(df[[outcome]])
    idx <- intersect(which(complete.cases(Y)), which(complete.cases(X)))
    X <- data.matrix(X[idx, ])
    Y <- Y[idx]
    
    # Fit ordinalNet model
    fit <- ordinalNet(
      x = X,
      y = Y,
      alpha = 0.5,                      # elastic net
      family = "cumulative",
      link = "logit",
      nLambda = n_lambdas,
      penaltyFactors = penalty.factor
    )
    
    # Extract the model with minimal BIC (for simplicity)
    best_idx <- which.min(fit$loglikPen)
    coefs <- coef(fit, matrix = TRUE)[[best_idx]]
    
    # Get selected features
    selected <- rownames(coefs)[apply(coefs, 1, function(x) any(x != 0))]
    selected <- setdiff(selected, c("(Intercept)", "(Thresholds)"))
    
    return(selected)
    
  } else if (model_type %in% c("cox", "binomial", "gaussian")) {
    
    # Outcome matrix
    if (model_type == "cox") {
      Y <- cbind(time = df[[time]], status = df[[outcome]])
      colnames(Y) <- c("time", "status")
      idx <- intersect(which(complete.cases(Y)), which(complete.cases(X)))
      X <- data.matrix(X[idx, ])
      Y <- data.matrix(Y[idx, ])
      
    } else {
      Y <- df[[outcome]]
      idx <- intersect(which(complete.cases(Y)), which(complete.cases(X)))
      X <- data.matrix(X[idx, ])
      Y <- as.numeric(Y[idx])
    }
    
    # Cross-validation setup
    if (cv_method == "loo") {
      folds <- nrow(X)
    } else if (cv_method != "kfold") {
      stop("Please select either LOO or k-fold CV. If k-fold, specify folds.")
    }
    
    if (!is.null(cores)) {
      require(doParallel)
      registerDoParallel(cores)
      p <- TRUE
    } else {
      p <- FALSE
    }
    
    list2env(list(X = X, Y = Y, n_alphas = n_alphas, n_lambdas = n_lambdas,
                  model_type = model_type, folds = folds, p = p,
                  max_coef = max_coef), .GlobalEnv)
    
    # Run ensr
    if (!is.null(max_coef)) {
      e <- ensr(X, Y, alphas = seq(0, 1, length = n_alphas),
                nlambda = n_lambdas, family = model_type,
                nfolds = folds, parallel = p, pmax = max_coef)
    } else {
      e <- ensr(X, Y, alphas = seq(0, 1, length = n_alphas),
                nlambda = n_lambdas, family = model_type,
                nfolds = folds, parallel = p)
    }
    
    # Select optimal model
    res <- summary(e)
    min_err <- min(res$cvm, na.rm = TRUE)
    se_err <- sd(res$cvm, na.rm = TRUE) / sqrt(sum(!is.na(res$cvm)))
    
    if (out == "min.error") {
      good_mods <- which.min(res$cvm)
    } else if (out == "1se.error") {
      good_mods <- which(res$cvm <= (min_err + se_err))
    }
    
    params <- data.frame(res[good_mods, ])
    params <- params[which(params$nzero == min(params$nzero)), ]
    params <- params[which.min(params$cvm), ]
    
    # Refit with selected alpha/lambda
    a <- params$alpha
    l <- params$lambda
    mod <- glmnet(y = Y, x = X, alpha = a, lambda = l,
                  family = model_type, penalty.factor = penalty.factor)
    
    selected <- as.matrix(coef(mod))
    selected <- rownames(selected)[selected[, 1] != 0]
    selected <- selected[selected != "(Intercept)"]
    
    # Clean up
    rm(X, Y, n_alphas, n_lambdas, model_type, folds, p, envir = .GlobalEnv)
    return(selected)
  } else {
    stop("model_type must be one of: 'gaussian', 'binomial', 'cox', 'ordinal'")
  }
}