#Process Simulation Results
#Load Libraries
source(fs::path(here::here(),"Simulation Code 11_04","Libraries.R"))

#Set up Directories
dir.dat <- c("/Users/hhampson/Library/Mobile Documents/com~apple~CloudDocs/BHRM_microbiome/Simulation Code 11_04/Simularion Results/Final")
dir.code <- c("/Users/hhampson/Library/Mobile Documents/com~apple~CloudDocs/BHRM_microbiome/Simulation Code 11_04/HPC Code")
dir.results <- c("/Users/hhampson/Library/CloudStorage/OneDrive-UW/Biostatistics Core Shared Drive/Hailey Hampson/1_Ongoing Projects/Microbiome Manuscript/Simulation Results/Processed Simulation Results")

#Load Taxa Z Matrices
#Genus
Z.s.g <- readRDS(fs::path(dir.code,"Z.s.g.RDS"))
#Family
Z.g.f <- readRDS(fs::path(dir.code,"Z.g.f.RDS"))
#Order
Z.f.o <- readRDS(fs::path(dir.code,"Z.f.o.RDS"))
#Class
Z.o.c <- readRDS(fs::path(dir.code,"Z.o.c.RDS"))
#Phylum
Z.c.p <- readRDS(fs::path(dir.code,"Z.c.p.RDS"))

#Load Results & Simulation Parameters ----
simulation.parameters <- readRDS(fs::path(dir.code,"Simulation_Scenarios_03_12_25.rds"))
simulation.scenarios <- read.csv(fs::path(dir.code,"Simulation_Parameters.csv")) %>% 
  dplyr::select(-X)

#Load simulation results
sim_num <- 3 #Set simulation number
# sim_result <- readRDS(fs::path(dir.dat,paste0("Simulation_",sim_num,"_Results.RDS")))
sim_result <- readRDS("/Users/hhampson/Library/Mobile Documents/com~apple~CloudDocs/BHRM_microbiome/Simulation Code 11_04/Simularion Results/all_iters_combined_scn2_1_to_123.RDS") %>% 
  mutate(component=ifelse(component=="Count model coefficients","Means",
                          ifelse(component=="Zero-inflation model coefficients","Probability",component)))

#Prepare causal taxa and exposure vectors
P.e <- simulation.scenarios$P.e[sim_num]
P.e.causal <- simulation.scenarios$P.e.causal[sim_num]
P.s.causal <- simulation.scenarios$P.s.causal[sim_num]
P.s.causal.all <- simulation.scenarios$P.s.causal.all[sim_num]
causal_species <- simulation.parameters[[P.s.causal.all]]
causal_species <- paste0("species",causal_species)
causal_exposures <- paste0("X.",1:P.e.causal)
individual_effect <- log(simulation.scenarios$OR.exposure[sim_num])
mixture_effect <- individual_effect*P.e.causal

# Get all unique genera containing at least one causal species
causal_genera <- colnames(Z.s.g)[colSums(Z.s.g[causal_species, , drop = FALSE]) > 0]
causal_families <- colnames(Z.g.f)[colSums(Z.g.f[causal_genera, , drop = FALSE]) > 0]
causal_orders <- colnames(Z.f.o)[colSums(Z.f.o[causal_families, , drop = FALSE]) > 0]
causal_classes <- colnames(Z.o.c)[colSums(Z.o.c[causal_orders, , drop = FALSE]) > 0]
causal_phyla <- colnames(Z.c.p)[colSums(Z.c.p[causal_classes, , drop = FALSE]) > 0]
causal_taxa <- c(causal_species,causal_genera,causal_families,causal_orders,causal_classes,causal_phyla)

#Separate ZINB, RBaHZING and BaHZING
# result_zinb <- sim_result %>% 
#   filter(model=="ZINB" | model=="Poisson")
# result_bahzing <- sim_result %>% 
#   filter(model=="BaHZING")
# result_rbahzing <- sim_result %>% 
#   filter(model=="RBaHZING")
# rm(sim_result) #remove full simulation dataframe to save memory 

# #Separate Means & ZI components
# zinb_m <- result_zinb %>% 
#   filter(component=="Means")
# zinb_p <- result_zinb %>% 
#   filter(component=="Probability")
# 
# 
# #Separate Mixture & Individual effects, Calculate the FDR & expected value of the means component for ZINB
# #Means
# zinb_m_mixture <- zinb_m %>% 
#   filter(exposure=="Mixture") %>% 
#   dplyr::group_by(domain) %>%
#   mutate(fdr = p.adjust(pval, method = "fdr")) %>%
#   ungroup() %>% 
#   mutate(causal=ifelse(taxa_full %in% causal_taxa & P.e.causal>0,1,0)) %>% 
#   mutate(exp_val=ifelse(causal==1,mixture_effect,0))
# 
# zinb_m_ind <- zinb_m %>% 
#   filter(exposure!="Mixture") %>% 
#   dplyr::group_by(domain) %>%
#   mutate(fdr = p.adjust(pval, method = "fdr")) %>%
#   ungroup() %>% 
#   mutate(causal=ifelse(taxa_full %in% causal_taxa & exposure %in% causal_exposures,1,0))  %>% 
#   mutate(exp_val=ifelse(causal==1,individual_effect,0))
# 
# #Calculate the FDR & expected value of the probability component for ZINB
# #Prob
# zinb_p_mixture <- zinb_p %>% 
#   filter(exposure=="Mixture") %>% 
#   dplyr::group_by(domain) %>%
#   mutate(fdr = p.adjust(pval, method = "fdr")) %>%
#   ungroup() %>% 
#   mutate(causal=ifelse(taxa_full %in% causal_taxa & P.e.causal>0,1,0)) %>% 
#   mutate(exp_val=ifelse(causal==1,mixture_effect,0))
# 
# zinb_p_ind <- zinb_p %>% 
#   filter(exposure!="Mixture") %>% 
#   dplyr::group_by(domain) %>%
#   mutate(fdr = p.adjust(pval, method = "fdr")) %>%
#   ungroup() %>% 
#   mutate(causal=ifelse(taxa_full %in% causal_taxa & exposure %in% causal_exposures,1,0))  %>% 
#   mutate(exp_val=ifelse(causal==1,individual_effect,0))

#Separate Means & ZI components
means <- sim_result %>% 
  filter(component=="Means")
prob <- sim_result %>% 
  filter(component=="Probability")


#Separate Mixture & Individual effects, Calculate the FDR & expected value of the means component for ZINB
#Means
means_mixture <- means %>% 
  filter(exposure=="Mixture") %>% 
  dplyr::group_by(domain) %>%
  mutate(fdr = p.adjust(pval, method = "fdr")) %>%
  ungroup() %>% 
  mutate(causal=ifelse(taxa_full %in% causal_taxa & P.e.causal>0,1,0)) %>% 
  mutate(exp_val=ifelse(causal==1,mixture_effect,0))

means_ind <- means %>% 
  filter(exposure!="Mixture") %>% 
  dplyr::group_by(domain) %>%
  mutate(fdr = p.adjust(pval, method = "fdr")) %>%
  ungroup() %>% 
  mutate(causal=ifelse(taxa_full %in% causal_taxa & exposure %in% causal_exposures,1,0))  %>% 
  mutate(exp_val=ifelse(causal==1,individual_effect,0))

#Calculate the FDR & expected value of the probability component for ZINB
#Prob
prob_mixture <- prob %>% 
  filter(exposure=="Mixture") %>% 
  dplyr::group_by(domain) %>%
  mutate(fdr = p.adjust(pval, method = "fdr")) %>%
  ungroup() %>% 
  mutate(causal=ifelse(taxa_full %in% causal_taxa & P.e.causal>0,1,0)) %>% 
  mutate(exp_val=ifelse(causal==1,mixture_effect,0))

prob_ind <- prob %>% 
  filter(exposure!="Mixture") %>% 
  dplyr::group_by(domain) %>%
  mutate(fdr = p.adjust(pval, method = "fdr")) %>%
  ungroup() %>% 
  mutate(causal=ifelse(taxa_full %in% causal_taxa & exposure %in% causal_exposures,1,0))  %>% 
  mutate(exp_val=ifelse(causal==1,individual_effect,0))




