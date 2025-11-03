 # code from Moerbeek

nr.cohorts <- 4 ### number of cohorts
var.e <- 5 ### variance random error
D <- matrix(0, 2, 2) ### covariance matrix
D[1, 1] <- 1 ### variance random intercept
D[2, 2] <- 1 ### variance random slope
D[1, 2] <- -0.1 ### covariance random intercept and slope
D[2, 1] <- -0.1
delta <- 0.9 ### standardized effect size

age1 <- c(1,2) ### age stages of measurement cohort 1
age2 <- c(2,3) ### age stages of measurement cohort 2
age3 <- c(3,4) ### age stages of measurement cohort 3
age4 <- c(4,5) ### age stages of measurement cohort 4

X1 <- cbind(rep(1,length(age1)),age1) ### design matrix fixed part cohort 1
Z1 <- X1 ### design matrix random part cohort 1
X2 <- cbind(rep(1,length(age2)),age2) ### design matrix fixed part cohort 2
Z2 <- X2 ### design matrix random part cohort 2
X3 <- cbind(rep(1,length(age3)),age3) ### design matrix fixed part cohort 3
Z3 <- X3 ### design matrix random part cohort 3
X4 <- cbind(rep(1,length(age4)),age4) ### design matrix fixed part cohort 4
Z4 <- X4 ### design matrix random part cohort 4

I1 <- diag(length(age1))
V1 <- Z1 %*% D %*% t(Z1) + var.e * I1 ### covariance matrix of responses cohort 1
I2 <- diag(length(age2))
V2 <- Z2 %*% D %*% t(Z2) + var.e * I2 ### covariance matrix of responses cohort 2
I3 <- diag(length(age3)) 
V3 <- Z3 %*% D %*% t(Z3) + var.e * I3 ### covariance matrix of responses cohort 3
I4 <- diag(length(age4)) 
V4 <- Z4 %*% D %*% t(Z4) + var.e * I4 ### covariance matrix of responses cohort 4

### calculate and plot power as a function of the total number of subjects
plot.new()
nr.subjects <- c(54)
power <- rep(0, length(nr.subjects))
for(ii in 1:length(nr.subjects)) {
  XVX <- ((t(X1) %*% solve(V1) %*% X1 + t(X2) %*%solve(V2) %*% X2
    + t(X3) %*% solve(V3) %*% X3 + t(X4) %*% solve(V4) %*% X4) * nr.subjects[ii])/nr.cohorts
  var.beta <- solve(XVX)
  power[ii] <- pnorm((delta * sqrt(D[2, 2]))/sqrt(var.beta[2, 2])- qnorm(0.975))
}
plot(nr.subjects, power, type = "l", ylim=c(0,1))
temp <- cbind(nr.subjects,power)
print(temp)

### calculate and plot power as a function of the total number of measurements
nr.measurements <- c(50, 100, 150, 200, 250, 300, 350, 400, 450, 500)
nr.subjects <- nr.measurements/(length(c(age1, age2, age3))/nr.cohorts)
power <- rep(0, length(nr.measurements))
for(ii in 1:length(nr.measurements)) {
  XVX <- ((t(X1) %*% solve(V1) %*% X1 + t(X2) %*%solve(V2) %*% X2
  + t(X3) %*% solve(V3) %*% X3 + t(X4) %*% solve(V4) %*% X4) * nr.subjects[ii])/nr.cohorts
  var.beta <- solve(XVX)
  power[ii] <- pnorm((delta  * sqrt(D[2, 2]))/sqrt(var.beta[2, 2])- qnorm(0.975))
}
plot(nr.measurements, power, type = "l", ylim=c(0,1))
