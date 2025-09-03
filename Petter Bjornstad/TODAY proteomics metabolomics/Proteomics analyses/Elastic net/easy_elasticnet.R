# This is a function for performing cross validation (CV) to select an optimal model 
# using the ElasticNet. By default uses leave one out (LOO) CV, but k-fold CV
# can also be used by setting cv_method = "kfold" and folds = k. See the trainControl 
# function in caret for additional details. There are two options for the output:
# out = "min.error" produces the model with the lowest CV error, and 
# out = "1se.error" produces all "acceptable" models (CV error within 
# 1 standard error of the minimum). 
easy_elasticnet = function(data,outcome,predictors,
                           n_alphas = 10,n_lambdas = 100,max_coef = NULL,
                           model_type = "gaussian",time = NULL,
                           cv_method = "loo",folds = NULL,out = "1se.error",
                           cores = 4,seed = 3654,penalty.factor=rep(1, ncol(predictors))){
  require(ensr)
  df = data
  # Random seed
  set.seed(seed)
  # Fix names if necessary
  #colnames(df) = make.names(colnames(df),unique = T,allow_ = F)
  #preds = make.names(predictors,unique = T,allow_ = F)
  #outcome = make.names(outcome,unique = T,allow_ = F)
  # Predictor matrix
  #X = data.frame(df[,preds])
  #X = data.frame(df[,predictors])
  X = predictors
  # Outcome matrix depending on model type
  if(model_type == "cox"){
    # Outcome matrix
    Y = cbind(time = df[,time], status = df[,outcome])
    colnames(Y) <- c("time","status")
    # Complete cases
    idx = intersect(which(complete.cases(Y)),which(complete.cases(X)))
    X = data.matrix(X[idx,])
    Y = data.matrix(Y[idx,])
    # Remove variables without any variance
    #near_zero = caret::nearZeroVar(X)
    #if(length(near_zero)>0){
    #  X = X[,-near_zero]
    #}
  } else if (model_type == "binomial" | model_type == "gaussian"){
    # Outcome matrix
    Y = df[,outcome]
    # Complete cases
    idx = intersect(which(complete.cases(Y)),which(complete.cases(X)))
    X = data.matrix(X[idx,])
    Y = as.numeric(Y[idx])
    # Remove variables without any variance
    #near_zero = caret::nearZeroVar(X)
    #if(length(near_zero)>0){
    #  X = X[,-near_zero]
    #}
  }
  # CV parameters
  if(cv_method == "loo"){
    folds = nrow(X)
  } else if (cv_method != "kfold"){
    stop("Please select either LOO or k-fold CV. If you are selecting k-fold CV, please specify the number of folds.")
  }
  # Parallel - recommended
  if(!is.null(cores)){
    require(doParallel)
    p = TRUE
    registerDoParallel(cores)
  } else {p = FALSE}
  # Add variables to global environment. Something about ensr doesn't work with 
  # function in function variable scoping.
  list2env(list(X=X,Y=Y,n_alphas=n_alphas,n_lambdas=n_lambdas,
                model_type=model_type,folds=folds,p=p,max_coef=max_coef),.GlobalEnv)
  # Grid search with glmnet - super slow
  #Y <- Surv(Y[time], Y[outcome])
  if(!is.null(max_coef)){
    e = ensr(X,Y,alphas = seq(0, 1, length = n_alphas),nlambda = n_lambdas,
             family = model_type,nfolds = folds,parallel = p,pmax = max_coef)
  } else {
    e = ensr(X,Y,alphas = seq(0, 1, length = n_alphas),nlambda = n_lambdas,
             family = model_type,nfolds = folds,parallel = p)
  }
  # Get alpha and lambdas
  res = summary(e)
  min_err = min(res$cvm,na.rm = T)
  se_err = sd(res$cvm,na.rm = T)/sqrt(sum(!is.na(res$cvm)))
  if (out == "min.error"){
    good_mods = which.min(res$cvm)
  } else if (out == "1se.error"){
    good_mods = which(res$cvm <= (min_err + se_err))
  }
  params = data.frame(res[good_mods,])
  params = params[which(params$nzero == min(params$nzero)),]
  params = params[which.min(params$cvm),]
  # Refit models to get selected parameters (the coef() function output for caret is confusing)
  a = params$alpha
  l = params$lambda
  mod = glmnet(y = Y,x = X,alpha = a,lambda = l,family = model_type, penalty.factor = penalty.factor)
  selected = as.matrix(coef(mod))
  selected = rownames(selected)[selected[,1] != 0]
  selected = selected[selected != "(Intercept)"]
  #selected = predictors[match(selected,preds)]
  # Remove variables from global environment, just in case
  rm(X,Y,n_alphas,n_lambdas,model_type,folds,p,envir = .GlobalEnv)
  # Return selected variables
  return(selected)
}
