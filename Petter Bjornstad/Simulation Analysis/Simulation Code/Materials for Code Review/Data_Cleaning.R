
# 1. Load Data ----
#Load Taxa Z Matrices
#Genus
# temp_file <- tempfile(fileext = ".RDS") 
# s3$download_file(bucket,"Z.s.g.RDS", temp_file)
# Z.s.g <- readRDS(temp_file)
Z.s.g <- readRDS(fs::path(dir.code,"Z.s.g.RDS"))

#Family
# temp_file <- tempfile(fileext = ".RDS") 
# s3$download_file(bucket,"Z.g.f.RDS", temp_file)
# Z.g.f <- readRDS(temp_file)
Z.g.f <- readRDS(fs::path(dir.code,"Z.g.f.RDS"))

#Order
# temp_file <- tempfile(fileext = ".RDS") 
# s3$download_file(bucket,"Z.f.o.RDS", temp_file)
# Z.f.o <- readRDS(temp_file)
Z.f.o <- readRDS(fs::path(dir.code,"Z.f.o.RDS"))

#Class
# temp_file <- tempfile(fileext = ".RDS") 
# s3$download_file(bucket,"Z.o.c.RDS", temp_file)
# Z.o.c <- readRDS(temp_file)
Z.o.c <- readRDS(fs::path(dir.code,"Z.o.c.RDS"))

#Phylum
# temp_file <- tempfile(fileext = ".RDS") 
# s3$download_file(bucket,"Z.c.p.RDS", temp_file)
# Z.c.p <- readRDS(temp_file)
Z.c.p <- readRDS(fs::path(dir.code,"Z.c.p.RDS"))

#3. Load Results & Simulation Parameters ----
#Load simulation scenarios
# temp_file <- tempfile(fileext = ".rds") 
# s3$download_file(bucket,"Simulation_Scenarios_03_12_25.rds", temp_file)
# simulation.parameters <- readRDS(temp_file)
simulation.parameters <- readRDS(fs::path(dir.code,"Simulation_Scenarios_03_12_25.rds"))

#Load all simulatin parameters
# temp_file <- tempfile(fileext = ".csv") 
# s3$download_file(bucket,"Simulation_Parameters.csv", temp_file)
# sim.par.all <- read.csv(temp_file)
sim.par.all <- read.csv(fs::path(dir.code,"Simulation_Parameters.csv")) %>%
# sim.par.all <- sim.par.all %>% 
  dplyr::select(-X) %>% 
  dplyr::select(c("P.e.causal","P.s.causal.all","N","OR.exposure")) %>% 
  dplyr::rename(P.s.scenario=P.s.causal.all)

#Read in & Format all results
all_formatted_results <- data.frame()
# for (scenario in 1:9) {
  for (scenario in 1) {
  results <- readRDS(fs::path(dir.dat,paste0("scn",scenario,"_iters1_to_1000.RDS")))
  formatted.results <- format.fxn(results)
  all_formatted_results <- rbind(all_formatted_results,formatted.results)
}
rm(formatted.results)

#5. Plot Results ----
# temp_file <- tempfile(fileext = ".RDS") 
# s3$download_file(bucket,"taxa_structure.RDS", temp_file)
# names <- readRDS(temp_file)
names <- readRDS(fs::path(dir.dat,"taxa_structure.RDS"))


