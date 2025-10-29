#Functions
render.geometric <- function(x, ...) {
  gm <- exp(mean(log(x[x > 0]), na.rm = TRUE))  # geometric mean
  gsd <- exp(sd(log(x[x > 0]), na.rm = TRUE))    # geometric SD
  c("", 
    " " = sprintf("%.2f (%.2f)", gm, gsd))
}