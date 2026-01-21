#Simulation Results ----
library(reticulate)
use_python("/home/hhampson/miniconda3/bin/python", required = TRUE)

#1. Set Directories ----
# dir.dat <- c("/Users/hhampson/Library/Mobile Documents/com~apple~CloudDocs/BHRM_microbiome/Simulation Code 11_04/Simularion Results/Final")
# dir.code <- c("/Users/hhampson/Library/Mobile Documents/com~apple~CloudDocs/BHRM_microbiome/Simulation Code 11_04/HPC Code")
# dir.lib <- c("/Users/hhampson/Library/Mobile Documents/com~apple~CloudDocs/BHRM_microbiome/Simulation Code 11_04")
# dir.results <- c("/Users/hhampson/Library/CloudStorage/OneDrive-UW/Biostatistics Core Shared Drive/Hailey Hampson/1_Ongoing Projects/Microbiome Manuscript/Simulation Results/Processed Simulation Results")


#2. Load Libraries ----
#Load Libraries
# source("Libraries.R")

#3. Set up to load data ----
## Load boto3 and pandas
boto3 <- reticulate::import("boto3")
pd <- reticulate::import("pandas")

## Create an S3 client
# install.packages("jsonlite")  # Install if not already installed
library(jsonlite)  # Load the package

keys <- fromJSON("/home/hhampson/keys.json") # replace with your Lambda username
session <- boto3$session$Session(
  aws_access_key_id = keys$MY_ACCESS_KEY,
  aws_secret_access_key = keys$MY_SECRET_KEY
)

## Create an S3 client with the session
s3 <- session$client("s3", endpoint_url = "https://s3.kopah.uw.edu")

#Bucket
bucket <- "simulation" # bucket name in Kopah

#Load Taxa Z Matrices
#Genus
temp_file <- tempfile(fileext = ".RDS") 
s3$download_file(bucket,"Z.s.g.RDS", temp_file)
Z.s.g <- readRDS(temp_file)
# Z.s.g <- readRDS(fs::path(dir.code,"Z.s.g.RDS"))

#Family
temp_file <- tempfile(fileext = ".RDS") 
s3$download_file(bucket,"Z.g.f.RDS", temp_file)
Z.g.f <- readRDS(temp_file)
# Z.g.f <- readRDS(fs::path(dir.code,"Z.g.f.RDS"))

#Order
temp_file <- tempfile(fileext = ".RDS") 
s3$download_file(bucket,"Z.f.o.RDS", temp_file)
Z.f.o <- readRDS(temp_file)
# Z.f.o <- readRDS(fs::path(dir.code,"Z.f.o.RDS"))

#Class
temp_file <- tempfile(fileext = ".RDS") 
s3$download_file(bucket,"Z.o.c.RDS", temp_file)
Z.o.c <- readRDS(temp_file)
# Z.o.c <- readRDS(fs::path(dir.code,"Z.o.c.RDS"))

#Phylum
temp_file <- tempfile(fileext = ".RDS") 
s3$download_file(bucket,"Z.c.p.RDS", temp_file)
Z.c.p <- readRDS(temp_file)
# Z.c.p <- readRDS(fs::path(dir.code,"Z.c.p.RDS"))

#3. Load Results & Simulation Parameters ----
#Load simulation scenarios
temp_file <- tempfile(fileext = ".rds") 
s3$download_file(bucket,"Simulation_Scenarios_03_12_25.rds", temp_file)
simulation.parameters <- readRDS(temp_file)
# simulation.parameters <- readRDS(fs::path(dir.code,"Simulation_Scenarios_03_12_25.rds"))

#Load all simulatin parameters
temp_file <- tempfile(fileext = ".csv") 
s3$download_file(bucket,"Simulation_Parameters.csv", temp_file)
sim.par.all <- read.csv(temp_file)
# sim.par.all <- read.csv(fs::path(dir.code,"Simulation_Parameters.csv")) %>% 
sim.par.all <- sim.par.all %>% 
  dplyr::select(-X) %>% 
  dplyr::select(c("P.e.causal","P.s.causal.all","N","OR.exposure")) %>% 
  dplyr::rename(P.s.scenario=P.s.causal.all)

#Read in & Format all results
all_formatted_results <- data.frame()
temp_file <- tempfile(fileext = ".RDS") 
# for (scenario in 1:9) {
  for (scenario in c(1)) {
  s3$download_file(bucket,paste0("scn",scenario,"_iters1_to_1000.RDS"), temp_file)
  results <- readRDS(temp_file)
  formatted.results <- format.fxn(results)
  all_formatted_results <- rbind(all_formatted_results,formatted.results)
}
rm(formatted.results)

#5. Plot Results ----
temp_file <- tempfile(fileext = ".RDS") 
s3$download_file(bucket,"taxa_structure.RDS", temp_file)
names <- readRDS(temp_file)