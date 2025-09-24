iAUC_glucose <- function(data,thresh,type = c("linear","spline"),
                         times = c(0,10,20,30,60,90,120,150,180),
                         cols = c("lab_ogtt_fasting",
                                  "min_10_glucose","min_20_glucose","min_30_glucose",
                                  "lab_ogtt_1_hour_glucose","min_90_glucose",
                                  "lab_ogtt_2_hour_glucose","min_150_glucose",
                                  "min_180_glucose")) {
  require(MESS)
  ogtt <- data[,cols] - thresh
  if (T %in% (rowSums(!is.na(ogtt)) < 2)) {
    stop("Too much missing data, all rows require at least two numeric values.")
  }
  aucs <- apply(ogtt,1, function(x) auc(times,x,type = type))
  return(aucs)
}

iAUC_insulin <- function(data,thresh,type = c("linear","spline"),
                         times = c(0,10,20,30,60,90,120,150,180),
                         cols = c("min_0_insulin","min_10_insulin","min_20_insulin","min_30_insulin","min_60_insulin",
                                  "min_90_insulin","min_120_insulin","min_150_insulin","min_180_insulin")) {
  require(MESS)
  ogtt <- data[,cols] - thresh
  if (T %in% (rowSums(!is.na(ogtt)) < 2)) {
    stop("Too much missing data, all rows require at least two numeric values.")
  }
  aucs <- apply(ogtt,1, function(x) auc(times,x,type = type))
  return(aucs)
}

iAUC_cpep <- function(data,thresh,type = c("linear","spline"),
                      times = c(0,10,20,30,60,90,120,150,180),
                      cols = c("min_0_c_peptide","min_10_c_peptide","min_20_c_peptide","min_30_c_peptide","min_60_c_peptide",
                               "min_90_c_peptide","min_120_c_peptide","min_150_c_peptide","min_180_c_peptide")) {
  require(MESS)
  ogtt <- data[,cols] - thresh
  if (T %in% (rowSums(!is.na(ogtt)) < 2)) {
    stop("Too much missing data, all rows require at least two numeric values.")
  }
  aucs <- apply(ogtt,1, function(x) auc(times,x,type = type))
  return(aucs)
}

#

iAUC_glucagon <- function(data,thresh,type = c("linear","spline"),
                          times = c(0,10,20,30,60,90,120,150,180),
                          cols = c("min_0_glucagon","min_10_glucagon","min_20_glucagon","min_30_glucagon","min_60_glucagon",
                                   "min_90_glucagon","min_120_glucagon","min_150_glucagon","min_180_glucagon")) {
  require(MESS)
  ogtt <- data[,cols] - thresh
  if (T %in% (rowSums(!is.na(ogtt)) < 2)) {
    stop("Too much missing data, all rows require at least two numeric values.")
  }
  aucs <- apply(ogtt,1, function(x) auc(times,x,type = type))
  return(aucs)
}

#spline <- iAUC(ogttdat,180,type = "spline")
#ogttdat$iauc180 <- iAUC(ogttdat,180,type = "linear")
#ogttdat$iauc0 <- iAUC(ogttdat,0,type = "linear")
#ogttdat$iaucbase <- iAUC(ogttdat,ogttdat$lab_ogtt_fasting,type = "linear")
