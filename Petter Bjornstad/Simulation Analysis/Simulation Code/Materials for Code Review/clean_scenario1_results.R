#Filter Scenario 2 out of Scenario 1 File
results <- readRDS(fs::path(dir.dat,paste0("scn",scenario,"_iters1_to_1000.RDS")))

#Check if there is more than scenario 1 in this file
unique(results$Scenario)
#[1] 1 2 
length(results$taxa_full) #17164744 before filtering

#Remove scenario 2 from this file and resave
results_1 <- results %>% 
  filter(Scenario==1)
length(results_1$taxa_full) #17147531 after filtering
# 17164744-17147531
# 17213 extra rows of results coded as Scenario 2 

#Save file
saveRDS(results_1,fs::path(dir.dat,paste0("scn1_iters1_to_1000_filtered.RDS")))
