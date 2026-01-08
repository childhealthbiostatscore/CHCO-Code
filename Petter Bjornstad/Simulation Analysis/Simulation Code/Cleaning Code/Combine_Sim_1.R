#Combine simulation scenarios 1 part A and 1 part B
simA <- readRDS("/Users/hhampson/Library/Mobile Documents/com~apple~CloudDocs/BHRM_microbiome/Simulation Code 11_04/Simularion Results/all_iters_combined_scn1_1to500.RDS")
simB <- readRDS("/Users/hhampson/Library/Mobile Documents/com~apple~CloudDocs/BHRM_microbiome/Simulation Code 11_04/Simularion Results/all_iters_combined_scn1_501to1000.RDS")

#Combine
sim1 <- rbind(simA,simB)
rm(simA,simB)

#Format sim 1 before saving
sim1$taxa_full <- str_replace_all(sim1$taxa_full,"k__X","species")

#Save combined file
saveRDS(sim1,"/Users/hhampson/Library/Mobile Documents/com~apple~CloudDocs/BHRM_microbiome/Simulation Code 11_04/Simularion Results/Final/Simulation_1_Results.RDS")
