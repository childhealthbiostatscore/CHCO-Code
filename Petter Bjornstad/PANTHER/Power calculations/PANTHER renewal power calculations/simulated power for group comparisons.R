library(nlme)
library(lme4)
library(ggplot2)
library(dplyr)

###########################
# 40 T1D vs. 60 control  #
###########################

# read in simulated data
# changed SD of random intercept to 5
data <- read.csv('/Users/lpyle/Library/CloudStorage/OneDrive-UW/Bjornstad/Biostatistics Core Shared Drive/Grants/Puberty and kidney structure and function R01 (PANTHER)/Renewal application/Power calculations/ald_sim_40_vs_60_per_group.csv')

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

###########################
# 40 T1D vs. 20 control  #
###########################

# read in simulated data
# changed SD of random intercept to 5
data <- read.csv('/Users/lpyle/Library/CloudStorage/OneDrive-UW/Bjornstad/Biostatistics Core Shared Drive/Grants/Puberty and kidney structure and function R01 (PANTHER)/Renewal application/Power calculations/ald_sim_40_vs_20_per_group.csv')

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

