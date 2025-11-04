#Process Simulation Results
#Load Libraries

#Set up Directories
dir.dat <- c("/Users/hhampson/Library/CloudStorage/OneDrive-UW/Biostatistics Core Shared Drive/Hailey Hampson/1_Ongoing Projects/Microbiome Manuscript/Simulation Study/Simulation Data")
dir.results <- c("/Users/hhampson/Library/CloudStorage/OneDrive-UW/Biostatistics Core Shared Drive/Hailey Hampson/1_Ongoing Projects/Microbiome Manuscript/Simulation Study/Simulation Results")

#Load Data
sim1a <- readRDS(fs::path(dir.dat,"all_iters_combined_scn1_1to500.RDS"))
sim1b <- readRDS(fs::path(dir.dat,"all_iters_combined_scn1_501to1000.RDS"))
sim1 <- rbind(sim1a,sim1b)
rm(sim1a,sim1b)

#