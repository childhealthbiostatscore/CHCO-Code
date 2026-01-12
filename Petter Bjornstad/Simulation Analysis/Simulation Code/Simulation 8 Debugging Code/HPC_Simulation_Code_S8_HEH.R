#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse)
})

opt_list <- list(
  make_option("--iter", type="integer"),
  make_option("--seed", type="integer"),
  make_option("--out",  type="character")
)
opt <- parse_args(OptionParser(option_list = opt_list))

stopifnot(!is.na(opt$iter), !is.na(opt$seed), !is.na(opt$out))

set.seed(opt$seed)

# Load Libraries ----
library(usethis)
library(rjags)
library(R2jags)
library(dplyr)
library(pscl)
library(purrr) 
library(qgcomp)
library(MASS)
library(BaHZING)
# EAH once the new version has been pushed, add the version here

use_git()

#Set working directory----
# setwd("/project2/jagoodri_1060/BaHZING_code")
dir <- c(here::here(),"HPC Simulation Code")
setwd(dir)

#Load functions----
source("/Users/hhampson/CHCO-Code/Petter Bjornstad/Simulation Analysis/Simulation Code/HPC Simulation Code/HPC_Simulation_Functions_New.R")
# use package functions, except for the ZING model

#Set scenario (1-9)
simulation_scenario <- 8
iter <- 475

#Load data----
Object <- readRDS("/Users/hhampson/Library/CloudStorage/OneDrive-SharedLibraries-UW/Laura Pyle - Biostatistics Core Shared Drive/Hailey Hampson/1_Ongoing Projects/Microbiome Manuscript/Simulation Data 01_08/HPC Data/FormattedObject.RDS")
Table <- data.frame(Object$Table)
Table <- Table[-149]
P.s.causal.all <- readRDS("/Users/hhampson/Library/CloudStorage/OneDrive-SharedLibraries-UW/Laura Pyle - Biostatistics Core Shared Drive/Hailey Hampson/1_Ongoing Projects/Microbiome Manuscript/Simulation Data 01_08/HPC Data/Simulation_Scenarios_03_12_25.rds")
parameters <- read.csv("/Users/hhampson/Library/CloudStorage/OneDrive-SharedLibraries-UW/Laura Pyle - Biostatistics Core Shared Drive/Hailey Hampson/1_Ongoing Projects/Microbiome Manuscript/Simulation Data 01_08/HPC Data/Simulation_Parameters.csv")
sim.par <- parameters %>% 
  filter(scenario==simulation_scenario)

#Set Parameters----
P.s <- 206
P.g <- 92
P.f <- 33
P.o <- 22
P.c <- 11
P.p <- 7
P.e <- sim.par$P.e
p_e_causal_value <- sim.par$P.e.causal
p_s_causal_value <- P.s.causal.all[[sim.par$P.s.causal.all]]
current_n <- sim.par$N
OR_exposure <- sim.par$OR.exposure
Corr <- sim.par$Corr

# Load Z matrices----
# Genus
Z.s.g <- readRDS("/Users/hhampson/Library/CloudStorage/OneDrive-SharedLibraries-UW/Laura Pyle - Biostatistics Core Shared Drive/Hailey Hampson/1_Ongoing Projects/Microbiome Manuscript/Simulation Data 01_08/HPC Data/Z.s.g.RDS")
Genus.R <- ncol(Z.s.g)
#Family
Z.g.f <- readRDS("/Users/hhampson/Library/CloudStorage/OneDrive-SharedLibraries-UW/Laura Pyle - Biostatistics Core Shared Drive/Hailey Hampson/1_Ongoing Projects/Microbiome Manuscript/Simulation Data 01_08/HPC Data/Z.g.f.RDS")
Family.R <- ncol(Z.g.f)
#Order
Z.f.o <- readRDS("/Users/hhampson/Library/CloudStorage/OneDrive-SharedLibraries-UW/Laura Pyle - Biostatistics Core Shared Drive/Hailey Hampson/1_Ongoing Projects/Microbiome Manuscript/Simulation Data 01_08/HPC Data/Z.f.o.RDS")
Order.R <- ncol(Z.f.o)
#Class
Z.o.c <- readRDS("/Users/hhampson/Library/CloudStorage/OneDrive-SharedLibraries-UW/Laura Pyle - Biostatistics Core Shared Drive/Hailey Hampson/1_Ongoing Projects/Microbiome Manuscript/Simulation Data 01_08/HPC Data/Z.o.c.RDS")
Class.R <- ncol(Z.o.c)
#Phylum
Z.c.p <- readRDS("/Users/hhampson/Library/CloudStorage/OneDrive-SharedLibraries-UW/Laura Pyle - Biostatistics Core Shared Drive/Hailey Hampson/1_Ongoing Projects/Microbiome Manuscript/Simulation Data 01_08/HPC Data/Z.c.p.RDS")
Phylum.R <- ncol(Z.c.p)

#Create data for simulation
Y.cont <- Table[, grep("d__", names(Table))]
Y.names <- colnames(Y.cont)
Y.obs <- apply(Y.cont, 2, FUN=function(v) { ifelse(v==0, 0, 1)})
Y.freq <- as.numeric(apply(Y.obs, 2, mean))
Y.mean <- apply(Y.cont, 2, FUN=function(v) { mean(v[v!=0])})
Y.mean <- ifelse(is.na(Y.mean), 0.001, Y.mean)
Y.freq <- Y.freq[Y.freq !=0 & Y.freq !=1]
Y.freq <- Y.freq[1:P.s] # frequency of species presence vs. absence
phi <- 5 

#Correlation Matrix function
create_corr_matrix <- function(P.e, Corr) {
  # Create a matrix filled with Corr
  corr_matrix <- matrix(Corr, nrow = P.e, ncol = P.e)
  
  # Set diagonal elements to 1
  diag(corr_matrix) <- 1
  
  return(corr_matrix)
}
corr_matrix <- create_corr_matrix(P.e, Corr)

#Run Bayesian, Ridge and Zero-Inflated Models ----
  i <- opt$iter
  
  #Simulate Data----
  X <- mvrnorm(n = current_n, mu = rep(0, P.e), Sigma = corr_matrix)
  beta.X <- c(rep(log(OR_exposure), p_e_causal_value), rep(0, (P.e - p_e_causal_value)))
  profiles <- rbind(rep(-0.5, P.e), rep(0.5, P.e))  # for counterfactual contrasts for mixture effects
  
  Y <- do.call(cbind, lapply(1:P.s, FUN=function(v) {
    if(v %in% p_s_causal_value) { etaY <- log(Y.freq[v]/(1-Y.freq[v])) + scale(X, scale = TRUE, center = TRUE)%*%beta.X}
    if(!v %in% p_s_causal_value) { etaY <- log(Y.freq[v]/(1-Y.freq[v])) }
    ProbI_nonzero <- exp(etaY)/(1+exp(etaY))
    ProbI_zero=1-ProbI_nonzero
    I_zero <- rbinom(current_n, 1, ProbI_zero)
    # simulate the negative binomial part
    if(v %in% p_s_causal_value) { eta_mu <- log(Y.mean[v]) + scale(X, scale = TRUE, center = TRUE)%*%beta.X }
    if(!v %in% p_s_causal_value) { eta_mu <- log(Y.mean[v]) }
    
    Y <- rnbinom(current_n, size=phi, mu=exp(eta_mu))
    # combine the two parts
    Y[I_zero==1] <- 0
    return(Y)
  }))  
  
  Y2 <- data.frame(Y)
  colnames(Y2) <- paste0("species",1:length(colnames(Y2)))
  
# format data for the ZING model

  Y.g <- as.data.frame(as.matrix(Y2)%*%as.matrix(Z.s.g))
  Y.f <- as.data.frame(as.matrix(Y.g)%*%as.matrix(Z.g.f))
  Y.o <- as.data.frame(as.matrix(Y.f)%*%as.matrix(Z.f.o))
  Y.c <- as.data.frame(as.matrix(Y.o)%*%as.matrix(Z.o.c))
  Y.p <- as.data.frame(as.matrix(Y.c)%*%as.matrix(Z.c.p))
  Y.g$id <- rownames(Y.g)
  Y.f$id <- rownames(Y.f)
  Y.o$id <- rownames(Y.o)
  Y.c$id <- rownames(Y.c)
  Y.p$id <- rownames(Y.p)
  Y2$id <- rownames(Y2)
  Y2 <- full_join(Y2,Y.g,by="id")
  Y2 <- full_join(Y2,Y.f,by="id")
  Y2 <- full_join(Y2,Y.o,by="id")
  Y2 <- full_join(Y2,Y.c,by="id")
  Y2 <- full_join(Y2,Y.p,by="id")
  X2 <- as.data.frame(X)
  X2$id <- rownames(X2)
  colnames(X2) <- c(paste0("X.",1:(length(colnames(X2))-1)),"id")
  Y2 <- full_join(Y2,X2,by="id")
  id_col <- which(colnames(X2)=="id")
  exposure.names <- colnames(X2)[-id_col]
   
  #Names vector
  Taxa.names <- c(rownames(Z.s.g),colnames(Z.s.g),colnames(Z.g.f),colnames(Z.f.o),colnames(Z.o.c),colnames(Z.c.p))
  # HW added: test iter 475
  # saveRDS(Y2, "simulated_data_iter_475.rds")
  Y2 <- readRDS("/Users/hhampson/Library/CloudStorage/OneDrive-SharedLibraries-UW/Laura Pyle - Biostatistics Core Shared Drive/Hailey Hampson/1_Ongoing Projects/Microbiome Manuscript/Simulation Data 01_08/Simulation 8 Debugging/simulated_data_iter_475.rds")
  Y2_test <- Y2
  #Run ZING Model ----
  results_list <- map(Taxa.names,dat=Y2,exposure=exposure.names,ZING_Model)
  # results_list <- map(Taxa.names, ~ ZING_Model(.x, dat = Y2, exposure = exposure.names))

  ZING.results <- bind_rows(results_list) 
  ZING.results$sig <- ifelse(ZING.results$P_Value < 0.05, "*", "N.S.")
  print(colnames(ZING.results))
  ZING.results$non_convergence_rate <- paste0(round((1-((length(ZING.results$Taxa))/(length(Taxa.names)*(P.e+1)*2)))*100,1),"%")
  ZING.results$iteration <- i
  
  # update the names so that the merge will work
  colnames(ZING.results) <- c("taxa_full", "exposure", "estimate", "SE", "pval",
                              "component", "model", "sig", "non_convergence_rate",
                              "iteration")
  # add columns that we need to have in BaHZING and RBaHZING even though we 
  # don't need them for ZING
  ZING.results$taxa_name <- ZING.results$taxa_full
  ZING.results$bci_lcl <- NA
  ZING.results$bci_ucl <- NA
  ZING.results$pdir <- NA
  ZING.results$prope <- NA
  ZING.results$pmap <- NA
  ZING.results$domain <- ifelse(substr(ZING.results$taxa_name, 1, 1) == "p", "Phylum",
                              ifelse(substr(ZING.results$taxa_name, 1, 1) == "c", "Class",
                                     ifelse(substr(ZING.results$taxa_name, 1, 1) == "o", "Order",
                                            ifelse(substr(ZING.results$taxa_name,1,1) == "f", "Family",
                                                   ifelse(substr(ZING.results$taxa_name,1,1) == "g", "Genus", "Species")))))
  
  # Format for BaHZING ----
  
  # make a file using the simulated data that has the same structure as
  # Format_BaHZING()
  
  # first, make a Y3 with only the species level data that start with k__
  Y3 <- data.frame(Y)
  colnames(Y3) <- paste0("species",1:length(colnames(Y3)))
  colnames(Y3) <- paste0("k__", colnames(Y3))
  Y3$id <- row.names(Y3)
  Y3 <- full_join(Y3, X2, by = "id")
  
  # must transpose all of the Z matrices
  z.s.g.t <- t(Z.s.g)
  z.g.f.t <- t(Z.g.f)
  z.f.o.t <- t(Z.f.o)
  z.o.c.t <- t(Z.o.c)
  z.c.p.t <- t(Z.c.p)
  
  # now add into a list
  
  sim_input <- list(
    Table = Y3,
    Species.Genus.Matrix = z.s.g.t,
    Genus.Family.Matrix = z.g.f.t,
    Family.Order.Matrix = z.f.o.t,
    Order.Class.Matrix = z.o.c.t,
    Class.Phylum.Matrix = z.c.p.t
  )

  #Run BaH-ZING Model----
  
  # try to run BaHZING with testing conditions
    r <- BaHZING::BaHZING_Model(sim_input, 
                       x = exposure.names, 
                       exposure_standardization = "none", # already standardized above
                       n.chains = 3,
                       n.adapt = 100,
                       n.iter.burnin = 1000,
                       n.iter.sample = 5000,
                       counterfactual_profiles = profiles
                       )
  results <- r %>%
    mutate(sig = ifelse((bci_lcl<0 & bci_ucl<0) | (bci_lcl>0 & bci_ucl>0), "*","N.S."))
  
  BaHZING.results <- results #5664 - 6348
  rm(r,results)
 
  BaHZING.results$non_convergence_rate <- paste0(round((1-((length(BaHZING.results$taxa_full))/(length(Taxa.names)*(P.e+1)*2)))*100,1),"%")
  BaHZING.results$iteration <- i
  BaHZING.results$SE <- NA
  BaHZING.results$pval <- NA
  BaHZING.results$model <- "BaHZING"
  
  #Run Ridge BaHZING Model ----
  # run the ridge model as implemented in the BaHZING package
  
  # testing parameters
  r2 <- BaHZING::Ridge_BaHZING_Model(sim_input, 
                            x = exposure.names, 
                            exposure_standardization = "none", # already standardized above
                            n.chains = 3,
                            n.adapt = 100,
                            n.iter.burnin = 1000,
                            n.iter.sample = 5000,
                            counterfactual_profiles = profiles
                            )
 
  results <- r2 %>%
    mutate(sig = ifelse((bci_lcl<0 & bci_ucl<0) | (bci_lcl>0 & bci_ucl>0), "*","N.S."))

  RBaHZING.results <- results #6348

  RBaHZING.results$non_convergence_rate <- paste0(round((1-((length(RBaHZING.results$taxa_full))/(length(Taxa.names)*(P.e+1)*2)))*100,1),"%")
  RBaHZING.results$iteration <- i
  RBaHZING.results$SE <- NA
  RBaHZING.results$pval <- NA
  RBaHZING.results$model <- "RBaHZING"
  rm(r2,results)
  # return(RBaHZING.results)
  Results.Output <- rbind(ZING.results,BaHZING.results)
  Results.Output <- rbind(Results.Output,RBaHZING.results)
  Results.Output$Scenario <- simulation_scenario

  # manually clean up
  rm(Y, Y2, Y.g, Y.f, Y.o, Y.c, Y.p, X2, ZING.results, BaHZING.results, RBaHZING.results)
  gc()
  
  # save the result
  saveRDS(Results.Output, file = opt$out)








