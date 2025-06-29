library(nlme)
library(lme4)
library(ggplot2)
library(dplyr)

###########################
# 100 cSFK vs. 40 control #
###########################

# read in simulated data, for lean vs. overweight T1D group comparison (N = 24 per group)
# changed SD of random intercept to 5
data <- read.csv("/Users/pylell/Library/CloudStorage/OneDrive-SharedLibraries-UW/Bjornstad/Biostatistics Core Shared Drive/Grants/LEOPARD power and stats/Revision/ald_sim_100_vs_40_per_group.csv")

# model checking and plots
ctrl <- lmeControl(maxIter = 1000)
quick <- data[data$sim==1,]
mod1 <- lmer(y ~ group*age + (1 + age|id), data=quick)
mod1 <- lme(y ~ group*age, random = ~1|id, data=quick)
quick %>% ggplot(aes(x=age,y=y,group=group)) + geom_line()

# power calcs
p <- NULL
for (i in 1:1000) {
  quick <- data[data$sim==i,]
  try(mod1 <- lme(y ~ group*age, random = ~1|id, data=quick, control=ctrl))
  pcur <- summary(mod1)$tTable[4,5]
  p <- c(p,pcur)
}
length(p[p<0.05])

