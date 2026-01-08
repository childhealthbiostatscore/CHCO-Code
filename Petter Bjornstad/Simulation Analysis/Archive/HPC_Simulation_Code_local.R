# Load Libraries ----
# Load Libraries ----
library(doParallel)
library(foreach)
library(rjags)
library(R2jags)
library(dplyr)
library(pscl)
library(purrr) 
library(qgcomp)
library(parallel)
library(fs)
library(stringr)
library(foreach)
library(parallel)
library(doRNG)
library(doParallel)
library(lme4)
library(lmerTest)

#Set number of cores for parallellization
#maxCores <- detectCores()
#numCores <- maxCores-1
#cl <- makeCluster(numCores)  # Create a cluster with the desired number of cores
#registerDoParallel(cl) 

#Local file path
#dir.dat <- c("/Volumes/Peds Endo/Petter Bjornstad")
#dir.dat2 <- c("/Volumes/Peds Endo/Petter Bjornstad/scRNA/data_clean")
#dir.code <- c("/Users/hhampson/Documents/CHCO-Code/Petter Bjornstad/Liver analysis/Liver scRNAseq")
#dir.results <- c("/Volumes/Peds Endo/Petter Bjornstad/Kidney Project/Results")

# #Lambda file path
# dir.dat <- c("/run/user/1026/gvfs/smb-share:server=ucdenver.pvt,share=som/PEDS/RI Biostatistics Core/Shared/Shared Projects/Laura/Peds Endo/Petter Bjornstad")
# dir.code <- c("/home/Github_Repo/CHCO-Code/Petter Bjornstad/Kidney scRNA/Kidney scRNA")
# dir.results <- c(fs::path(dir.dat,"Kidney Project/Results"))

# #Mac Studio File Path
# dir.dat <- c("/Users/hhampson/Library/CloudStorage/OneDrive-UW/Biostatistics Core Shared Drive")
# dir.results <- c("/Users/hhampson/Library/CloudStorage/OneDrive-UW/Biostatistics Core Shared Drive/Kidney scRNAseq Project/Results")
# dir.ipa <- c("/Users/hhampson/Documents/IPA/Results")

#Load functions----
source("HPC_Simulation_Functions_03_12.R")

#Set location to store results ----
dir.results <- c("/Users/hhampson/Library/CloudStorage/OneDrive-UW/Biostatistics Core Shared Drive/Simulation Analysis")
dir.data <- c("/Users/hhampson/Library/CloudStorage/OneDrive-UW/Biostatistics Core Shared Drive/Simulation Analysis/Simulation_Data")

#Load data----
#Formatted Taxonomy Information
# Object <- readRDS("FormattedObject.RDS")
invisible(gc())
# bucket <- "simulation" # bucket name in Kopah
# temp_file <- tempfile(fileext = ".RDS") # need to create a temporary file
# s3$download_file(bucket, "Simulation_Data/FormattedObject.RDS", temp_file)
Object <- readRDS(fs::path(dir.data,"FormattedObject.RDS"))
invisible(gc())
Table <- data.frame(Object$Table)
Table <- Table[-149]

#Simulation Scenarios & Parameters
# P.s.causal.all <- readRDS("Simulation_Scenarios_03_12_25.rds")
# invisible(gc())
# temp_file <- tempfile(fileext = ".rds") # need to create a temporary file
# s3$download_file(bucket, "Simulation_Data/Simulation_Scenarios_03_12_25.rds", temp_file)
P.s.causal.all <- readRDS(fs::path(dir.data,"Simulation_Scenarios_03_12_25.rds"))
# invisible(gc())

#Set Parameters----
P.s <- 206
P.g <- 92
P.f <- 33
P.o <- 22
P.c <- 11
P.p <- 7

# Load Z matrices----
# Genus
# Z.s.g <- readRDS("/Users/hhampson/Dropbox (USC Lab)/Chatzi Projects (Active)/Env Chem SOL-CHS/Analysis/2_ongoing/CHS PFAS Microbiome (Hailey)/BHRM_microbiome/Formatted Data/Z.s.g.RDS")
Z.s.g <- readRDS(fs::path(dir.data,"Z.s.g.RDS"))
# temp_file <- tempfile(fileext = ".RDS") # need to create a temporary file
# s3$download_file(bucket, paste0("Simulation_Data/Z.s.g.RDS"), temp_file)
# Z.s.g <- readRDS(temp_file)
# invisible(gc())
# Genus.R <- ncol(Z.s.g)
#Family
# Z.g.f <- readRDS("/Users/hhampson/Dropbox (USC Lab)/Chatzi Projects (Active)/Env Chem SOL-CHS/Analysis/2_ongoing/CHS PFAS Microbiome (Hailey)/BHRM_microbiome/Formatted Data/Z.g.f.RDS")
Z.g.f <- readRDS(fs::path(dir.data,"Z.g.f.RDS"))
# temp_file <- tempfile(fileext = ".RDS") # need to create a temporary file
# s3$download_file(bucket, paste0("Simulation_Data/Z.g.f.RDS"), temp_file)
# Z.g.f <- readRDS(temp_file)
# invisible(gc())
# Family.R <- ncol(Z.g.f)
#Order
# Z.f.o <- readRDS("/Users/hhampson/Dropbox (USC Lab)/Chatzi Projects (Active)/Env Chem SOL-CHS/Analysis/2_ongoing/CHS PFAS Microbiome (Hailey)/BHRM_microbiome/Formatted Data/Z.f.o.RDS")
Z.f.o <- readRDS(fs::path(dir.data,"Z.f.o.RDS"))
# temp_file <- tempfile(fileext = ".RDS") # need to create a temporary file
# s3$download_file(bucket, paste0("Simulation_Data/Z.f.o.RDS"), temp_file)
# Z.f.o <- readRDS(temp_file)
# invisible(gc())
# Order.R <- ncol(Z.f.o)
#Class
# Z.o.c <- readRDS("/Users/hhampson/Dropbox (USC Lab)/Chatzi Projects (Active)/Env Chem SOL-CHS/Analysis/2_ongoing/CHS PFAS Microbiome (Hailey)/BHRM_microbiome/Formatted Data/Z.o.c.RDS")
Z.o.c <- readRDS(fs::path(dir.data,"Z.o.c.RDS"))
# temp_file <- tempfile(fileext = ".RDS") # need to create a temporary file
# s3$download_file(bucket, paste0("Simulation_Data/Z.o.c.RDS"), temp_file)
# Z.o.c <- readRDS(temp_file)
# invisible(gc())
# Class.R <- ncol(Z.o.c)
#Phylum
# Z.c.p <- readRDS("/Users/hhampson/Dropbox (USC Lab)/Chatzi Projects (Active)/Env Chem SOL-CHS/Analysis/2_ongoing/CHS PFAS Microbiome (Hailey)/BHRM_microbiome/Formatted Data/Z.c.p.RDS")
Z.c.p <- readRDS(fs::path(dir.data,"Z.c.p.RDS"))
# temp_file <- tempfile(fileext = ".RDS") # need to create a temporary file
# s3$download_file(bucket, paste0("Simulation_Data/Z.c.p.RDS"), temp_file)
# Z.c.p <- readRDS(temp_file)
# invisible(gc())
# Phylum.R <- ncol(Z.c.p)

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

model.fxn <- function(i) {
  # i = 1
  #Set seed for reproducibility
  set.seed(i)
  
  #Simulate Data----
  sim.output <- sim.fxn(i)
  
  # Save the simulated dataset for future use
  # saveRDS(sim.output, file = paste0("Dataset_",scenario,"_scenario_", k, "_", i, ".RDS"))
  
  #Extract outputs from simulated data
  X <- sim.output$X
  Y <- sim.output$Y
  profiles <- sim.output$profiles
  OR_exposure <- sim.par$OR.exposure[i]
  N <- sim.par$N[i]
  Y2 <- data.frame(Y)
  colnames(Y2) <- paste0("species",1:length(colnames(Y2)))
  
  #Format Data for ZING Model 
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
  
  #Run ZING Model ----
  results_list <- map(Taxa.names,dat=Y2,exposure=exposure.names,ZING_Model)
  ZING.results <- bind_rows(results_list) 
  # NonConvergence <- paste0(round((1-(length(ZING.results$Taxa)/(371*4*2)))*100,0),"%")
  # return(ZING.results)
  
  #Run BaH-ZING Model----
  jdata <- list(N=N, Y=Y, P.s=P.s, X=X, P.e=P.e,
                Z.s.g=Z.s.g, P.g=P.g,
                P.f=P.f, Z.g.f=Z.g.f,
                P.o=P.o, Z.f.o=Z.f.o,
                P.c=P.c, Z.o.c=Z.o.c,
                P.p=P.p, Z.c.p=Z.c.p,
                profiles=profiles)
  var.s <- c("species.beta", "genus.beta", "family.beta", "order.beta", "class.beta", "phylum.beta","species.psi","genus.psi","family.psi","order.psi","class.psi","phylum.psi","species.beta.zero", "genus.beta.zero", "family.beta.zero", "order.beta.zero", "class.beta.zero", "phylum.beta.zero","species.psi.zero","genus.psi.zero","family.psi.zero","order.psi.zero","class.psi.zero","phylum.psi.zero","omega","disp")
  # model.fit <- jags.model(file=textConnection(BHRM.microbiome), data=jdata, n.chains=1, n.adapt=100, quiet=F)
  # update(model.fit, n.iter=10)
  # model.fit <- coda.samples(model=model.fit, variable.names=var.s, n.iter=50, thin=1, progress.bar="none")
  model.fit <- jags.model(file=textConnection(BHRM.microbiome), data=jdata, n.chains=3, n.adapt=100, quiet=F)
  update(model.fit, n.iter=1000)
  model.fit <- coda.samples(model=model.fit, variable.names=var.s, n.iter=5000, thin=1, progress.bar="none")
  # # summarize results
  r <- summary(model.fit)
  results <- data.frame(round(r$statistics[,1:2],3), round(r$quantiles[,c(1,5)],3))
  results <- results %>%
    mutate(sig = ifelse((X2.5.<0 & X97.5.<0) | (X2.5.>0 & X97.5.>0), "*","N.S."))
  # mutate(Odds.Ratio = exp(Mean),
  #        lcl=exp(X2.5.),
  #        ucl=exp(X97.5.))
  BaHZING.results <- results #5664 - 6348
  rm(model.fit,r,results)
  # return(BaHZING.results)
  # BaHZING.results <- results
  BaHZING.results <- BaHZING.results %>% 
    mutate(Taxa=rownames(BaHZING.results),
           Exposure=rownames(BaHZING.results)) %>%
    mutate(SE=SD/sqrt(N)) %>% 
    rename(P_Value=sig) %>% 
    mutate(Component=ifelse(grepl("zero",Taxa),"Probability","Means")) %>% 
    mutate(Model="BaH-ZING") %>% 
    # rename(OddsRatio=Odds.Ratio) %>% 
    dplyr::select(Taxa,Exposure,Mean,SE,P_Value,Component,Model)
  rownames(BaHZING.results) <- NULL
  
  #Run Ridge BaHZING Model ----
  jdata <- list(N=N, Y=Y, P.s=P.s, X=X, P.e=P.e,
                P.g=P.g,
                P.f=P.f,
                P.o=P.o,
                P.c=P.c,
                P.p=P.p,
                profiles=profiles)
  var.s <- c("species.beta", "genus.beta", "family.beta", "order.beta", "class.beta", "phylum.beta","species.psi","genus.psi","family.psi","order.psi","class.psi","phylum.psi","species.beta.zero", "genus.beta.zero", "family.beta.zero", "order.beta.zero", "class.beta.zero", "phylum.beta.zero","species.psi.zero","genus.psi.zero","family.psi.zero","order.psi.zero","class.psi.zero","phylum.psi.zero","omega","disp")
  # model.fit <- jags.model(file=textConnection(RBHRM.microbiome), data=jdata, n.chains=1, n.adapt=100, quiet=F)
  # update(model.fit, n.iter=10)
  # model.fit <- coda.samples(model=model.fit, variable.names=var.s, n.iter=50, thin=1, progress.bar="none")
  model.fit <- jags.model(file=textConnection(RBHRM.microbiome), data=jdata, n.chains=3, n.adapt=100, quiet=F)
  update(model.fit, n.iter=1000)
  model.fit <- coda.samples(model=model.fit, variable.names=var.s, n.iter=5000, thin=1, progress.bar="none")
  # summarize results
  r <- summary(model.fit)
  results <- data.frame(round(r$statistics[,1:2],3), round(r$quantiles[,c(1,5)],3))
  results <- results %>%
    filter(grepl("species",rownames(results)))
  results <- results %>%
    mutate(sig = ifelse((X2.5.<0 & X97.5.<0) | (X2.5.>0 & X97.5.>0), "*","N.S."))
  # mutate(Odds.Ratio = exp(Mean),
  #        lcl=exp(X2.5.),
  #        ucl=exp(X97.5.))
  RBaHZING.results <- results #6348
  RBaHZING.results <- RBaHZING.results %>%
    mutate(Taxa=rownames(RBaHZING.results),
           Exposure=rownames(RBaHZING.results)) %>%
    # mutate(SE=SD) %>%
    mutate(SE=SD/sqrt(N)) %>%
    # select(-c(X2.5.,X97.5.,lcl,ucl)) %>%
    rename(P_Value=sig) %>%
    mutate(Component=ifelse(grepl("zero",Taxa),"Probability","Means")) %>%
    mutate(Model="RBaH-ZING") %>%
    # rename(OddsRatio=Odds.Ratio) %>%
    dplyr::select(Taxa,Exposure,Mean,SE,P_Value,Component,Model)
  rownames(RBaHZING.results) <- NULL
  rm(model.fit,r,results)
  # return(RBaHZING.results)
  Results.Output <- rbind(ZING.results,BaHZING.results)
  Results.Output <- rbind(Results.Output,RBaHZING.results)
  return(Results.Output)
}

# Get the number of cores your machine has
num_cores <- 4  # Reserve one core for the system, adjust as needed

# Function to run the model for each k and scenario
run_simulation <- function(k, scenario, dir.data, dir.results) {
  # Load simulation parameter data
  sim.par <- readRDS(fs::path(dir.data, paste0("Simulation_Parameter_", scenario, ".rds")))
  P.e <- sim.par$P.e[[1]]
  
  # Run the model 100 times
  full_results <- list()  # to store results
  for (i in 1:100) {  # You can adjust the number of iterations here
    result <- model.fxn(k)  # Run the model function
    full_results[[i]] <- result
  }
  
  # Combine all results into a single data frame or object
  full_results_combined <- do.call(rbind, full_results)
  
  # Ensure the results directory exists, create it if not
  if (!dir.exists(dir.results)) {
    dir.create(dir.results, recursive = TRUE)
  }
  
  # Define the output file path
  output_file <- fs::path(dir.results, paste0("Full_results_scenario_", scenario, "_iteration_", k, "_04_07_2025.csv"))
  
  # Save the results to the CSV file
  write.csv(full_results_combined, file = output_file)
  
  # Return a success message indicating the completion of this specific scenario and iteration
  return(paste("Simulation complete for k =", k, "and scenario =", scenario, ". Results saved to:", output_file))
}

# Set up parallel execution (using mclapply for Unix-based systems)
# If you're using Windows, you should replace mclapply with parLapply.
all_combinations <- expand.grid(k = 1:10, scenario = 1:40)  # Adjust the ranges as needed

# Parallel execution on all combinations of k and scenario
results <- mclapply(1:nrow(all_combinations), function(i) {
  k <- all_combinations$k[i]
  scenario <- all_combinations$scenario[i]
  
  # Call the run_simulation function for each combination
  run_simulation(k, scenario, dir.data, dir.results)
}, mc.cores = num_cores)

# Print results to check that it worked
print(results)
