#Functions
#1. Format Raw Simulation Results
# results.data <- sim1
# format.fxn <- function(results.data){

scenario_num <- 1
results <- zinb_means
calculate_expected_and_fdr <- function(results, scenario_num) {
  
  # 1. Extract scenario parameters
  params <- simulation.scenarios[simulation.scenarios$Scenario == scenario_num, ]
  
  n_causal_exposures <- params$P.e.causal
  n_causal_species_index <- params$P.s.causal.all
  or_exposure <- params$OR.exposure
  
  log_or <- log(or_exposure)
  
  # 2. Identify causal exposures
  if (n_causal_exposures > 0) {
    causal_exposures <- paste0("X.", 1:n_causal_exposures)
  } else {
    causal_exposures <- character(0)
  }
  
  # 3. Identify causal species
  causal_species_nums <- simulation.parameters[[n_causal_species_index]]
  
  if (!is.null(causal_species_nums)) {
    causal_species <- paste0("species", causal_species_nums)
  } else {
    causal_species <- character(0)
  }
  
  # 4. Identify causal taxa at higher levels
  get_causal_higher_taxa <- function(Z_matrix, causal_lower_taxa) {
    if (is.null(Z_matrix) || length(causal_lower_taxa) == 0) {
      return(character(0))
    }
    causal_rows <- rownames(Z_matrix) %in% causal_lower_taxa
    if (sum(causal_rows) == 0) {
      return(character(0))
    }
    causal_cols <- colSums(Z_matrix[causal_rows, , drop = FALSE]) > 0
    return(colnames(Z_matrix)[causal_cols])
  }
  
  causal_genera <- get_causal_higher_taxa(Z.s.g, causal_species)
  causal_families <- get_causal_higher_taxa(Z.g.f, causal_genera)
  causal_orders <- get_causal_higher_taxa(Z.f.o, causal_families)
  causal_classes <- get_causal_higher_taxa(Z.o.c, causal_orders)
  causal_phyla <- get_causal_higher_taxa(Z.c.p, causal_classes)
  
  # 5. Calculate expected values with proportion weighting
  
  # Calculate base expected value (mixture-level) for a higher taxon
  calc_base_expected <- function(taxa_name, Z_matrix, causal_children, log_or) {
    if (is.null(Z_matrix) || !taxa_name %in% colnames(Z_matrix)) {
      return(0)
    }
    total <- sum(Z_matrix[, taxa_name] == 1)
    if (total == 0) return(0)
    
    children <- rownames(Z_matrix)[Z_matrix[, taxa_name] == 1]
    n_causal <- sum(children %in% causal_children)
    
    return(log_or * (n_causal / total))
  }
  
  # Initialize expected_value column
  results$expected_value <- 0
  
  # Only calculate if there are causal exposures and causal species
  if (n_causal_exposures > 0 && length(causal_species) > 0) {
    
    # Get unique taxa at each level
    species_names <- unique(results$taxa_name[results$domain == "Species"])
    genus_names <- unique(results$taxa_name[results$domain == "Genus"])
    family_names <- unique(results$taxa_name[results$domain == "Family"])
    order_names <- unique(results$taxa_name[results$domain == "Order"])
    class_names <- unique(results$taxa_name[results$domain == "Class"])
    phylum_names <- unique(results$taxa_name[results$domain == "Phylum"])
    
    # Calculate BASE expected values (= mixture expected) for each taxon
    # Species level - binary (causal or not)
    expected_species <- data.frame(
      taxa_name = species_names,
      base_expected = ifelse(species_names %in% causal_species, log_or, 0)
    )
    
    # Genus level - proportion of causal species
    expected_genus <- data.frame(
      taxa_name = genus_names,
      base_expected = sapply(genus_names, function(g) {
        calc_base_expected(g, Z.s.g, causal_species, log_or)
      })
    )
    
    # Family level - proportion of causal genera
    expected_family <- data.frame(
      taxa_name = family_names,
      base_expected = sapply(family_names, function(f) {
        calc_base_expected(f, Z.g.f, causal_genera, log_or)
      })
    )
    
    # Order level - proportion of causal families
    expected_order <- data.frame(
      taxa_name = order_names,
      base_expected = sapply(order_names, function(o) {
        calc_base_expected(o, Z.f.o, causal_families, log_or)
      })
    )
    
    # Class level - proportion of causal orders
    expected_class <- data.frame(
      taxa_name = class_names,
      base_expected = sapply(class_names, function(c) {
        calc_base_expected(c, Z.o.c, causal_orders, log_or)
      })
    )
    
    # Phylum level - proportion of causal classes
    expected_phylum <- data.frame(
      taxa_name = phylum_names,
      base_expected = sapply(phylum_names, function(p) {
        calc_base_expected(p, Z.c.p, causal_classes, log_or)
      })
    )
    
    # Combine all expected values
    all_expected <- rbind(expected_species, expected_genus, expected_family,
                          expected_order, expected_class, expected_phylum)
    
    # Join to results and calculate expected values
    results <- results %>%
      left_join(all_expected, by = "taxa_name") %>%
      mutate(
        expected_value = case_when(
          # Mixture gets the full base expected value
          exposure == "Mixture" ~ base_expected,
          # Individual causal exposures get base_expected / n_causal_exposures
          exposure %in% causal_exposures ~ base_expected / n_causal_exposures,
          # Non-causal exposures get 0
          TRUE ~ 0
        )
      ) %>%
      select(-base_expected)
  }
  
  # 6. Calculate within-level FDR
  results <- results %>%
    group_by(domain) %>%
    mutate(
      pval_adj_fdr = p.adjust(pval, method = "BH"),
      sig_fdr = pval_adj_fdr < 0.05
    ) %>%
    ungroup()
  
  # 7. Add indicator columns
  results$is_causal <- results$expected_value != 0
  results$sig_raw <- results$pval < 0.05
  
  results$TP <- results$is_causal & results$sig_fdr
  results$FP <- !results$is_causal & results$sig_fdr
  results$TN <- !results$is_causal & !results$sig_fdr
  results$FN <- results$is_causal & !results$sig_fdr
  
  # 8. Calculate summary statistics by domain
  summary_stats <- results %>%
    group_by(domain) %>%
    summarise(
      n_total = n(),
      n_causal = sum(is_causal),
      n_null = sum(!is_causal),
      
      TP = sum(TP),
      FP = sum(FP),
      TN = sum(TN),
      FN = sum(FN),
      
      MSE = mean((estimate - expected_value)^2, na.rm = TRUE),
      
      precision = ifelse((TP + FP) > 0, TP / (TP + FP), NA_real_),
      sensitivity = ifelse((TP + FN) > 0, TP / (TP + FN), NA_real_),
      specificity = ifelse((TN + FP) > 0, TN / (TN + FP), NA_real_),
      FDR = ifelse((FP + TP) > 0, FP / (FP + TP), NA_real_),
      
      .groups = "drop"
    ) %>%
    mutate(scenario = scenario_num)
  
  # 9. Return both results and summary
  return(list(
    results = results,
    summary = summary_stats
  ))
}


###############################################################################

params <- simulation.scenarios[simulation.scenarios$Scenario == scenario_num, ]

n_causal_exposures <- params$P.e.causal
n_causal_species_index <- params$P.s.causal.all  # Index into simulation.parameters list
or_exposure <- params$OR.exposure

# Expected effect on log scale for ZINB/Poisson means model
expected_effect <- log(or_exposure)


zinb_means <- results.data %>%
  filter(model=="ZINB" | model =="Poisson") %>%
  filter(component=="Means") %>%
  dplyr::select(-c("bci_lcl","bci_ucl","pdir","prope","pmap")) %>% 
  
  
  
  # mutate(P_Value = as.numeric(P_Value)) %>%
  mutate(fdr.p=case_when(grepl("species",Taxa) ~ P_Value*206,
                         grepl("genus",Taxa) ~ P_Value*92,
                         grepl("family",Taxa) ~ P_Value*33,
                         grepl("order",Taxa) ~ P_Value*22,
                         grepl("class",Taxa) ~ P_Value*11,
                         grepl("phylum",Taxa) ~ P_Value*7)) %>%
  mutate(P_Value=ifelse(P_Value<0.05,"*","N.S")) %>%
  mutate(fdr.p=ifelse(fdr.p<0.05,"*","N.S.")) %>%
  mutate(Model="ZINB") 
# mutate(SE=SD) %>%
# mutate(SD=SE*sqrt(100))

zinb_prob <- results.data %>%
  filter(model=="ZINB" | model =="Poisson") %>%
  filter(component=="Probability")

# formatted.data <- results.data %>%
#   filter(Component=="Means") %>%
#   filter(Model=="BaH-ZING" | Model =="RBaH-ZING") %>%
#   mutate(Exposure=str_replace(Exposure,"^[^,]*,",""),
#          Exposure=str_replace(Exposure,"]",""),
#          Exposure=ifelse(grepl("Mixture",Exposure) | grepl("X.",Exposure),Exposure,paste0("X.",Exposure)),
#          Exposure=ifelse(grepl("psi",Exposure),"Mixture",Exposure),
#          Exposure=ifelse(grepl("omega",Exposure),"Omega",Exposure),
#          Exposure=ifelse(grepl("disp",Exposure),"Dispersion",Exposure)) %>%
#   mutate(Taxa=str_replace(Taxa,",.*$",""),
#          Level=str_replace(Taxa,"\\[.*$",""),
#          Level=str_replace(Level,"\\..*$",""),
#          Taxa=str_replace(Taxa,".*?\\[",""),
#          Taxa=str_replace(Taxa,"]",""),
#          Taxa=paste0(Level,Taxa)) %>%
#   select(-Level) %>%
#   filter(Exposure!="Omega") %>%
#   filter(Exposure!="Dispersion") %>%
#   mutate(fdr.p = P_Value) 
# # mutate(SE=SD/sqrt(100))
# 
# formatted.data <- rbind(zing.data,formatted.data)
# formatted.data <- formatted.data %>%
#   mutate(Taxa.Level=case_when(grepl("species",Taxa) ~ "Species",
#                               grepl("genus",Taxa) ~ "Genus",
#                               grepl("family",Taxa) ~ "Family",
#                               grepl("order",Taxa) ~ "Order",
#                               grepl("class",Taxa) ~ "Class",
#                               grepl("phylum",Taxa) ~ "Phylum"))
# #Remove columns not needed
# formatted.data <- formatted.data %>%
#   select(-Component)
# 
# # Define a function to assign replication numbers
# assign_replication_numbers <- function(group) {
#   group$Replication <- seq_len(nrow(group))
#   return(group)
# }
# 
# # Apply the function to each group
# formatted.data <- formatted.data %>%
#   tidylog::group_by(Taxa, Exposure, Model) %>%
#   group_modify(~ assign_replication_numbers(.x))
# 
# #Group number for plotting
# formatted.data <- formatted.data %>%
#   tidylog::group_by(Taxa,Exposure,Replication) %>%
#   tidylog::mutate(group_number = cur_group_id()) %>%
#   ungroup()
# 
# sim.par <- sim.par.all[unique(formatted.data$Scenario),]
# par <- sim.par$P.s.scenario
# #Causal Exposures
# if (sim.par$P.e.causal==0){
#   causal.e <- NULL
# } else {
#   causal.e <- paste0("X.",rep(1:sim.par$P.e.causal))
# }
# 
# #Species level indicator
# if(par==1){
#   causal.s <- NULL
# }
# if(par!=1){
#   causal.s <- paste0("species",simulation.parameters[[par]])
# }
# 
# #Genus level indicator
# is_in_causal_s <- which(rownames(Z.s.g) %in% causal.s)
# subset_matrix <- Z.s.g[is_in_causal_s, ]
# genus_occurrences <- colSums(subset_matrix)
# causal.g <- names(which(genus_occurrences > 0))
# 
# #Family level indicator
# is_in_causal_g <- which(rownames(Z.g.f) %in% causal.g)
# subset_matrix <- Z.g.f[is_in_causal_g, ]
# family_occurrences <- colSums(subset_matrix)
# causal.f <- names(which(family_occurrences > 0))
# #Order level indicator
# is_in_causal_f <- which(rownames(Z.f.o) %in% causal.f)
# subset_matrix <- Z.f.o[is_in_causal_f, ]
# order_occurrences <- colSums(subset_matrix)
# causal.o <- names(which(order_occurrences > 0))
# # Class level indicator
# is_in_causal_c <- which(rownames(Z.o.c) %in% causal.o)
# subset_matrix <- Z.o.c[is_in_causal_c, ]
# class_occurrences <- colSums(subset_matrix)
# causal.c <- names(which(class_occurrences > 0))
# #Phylum level indicator
# is_in_causal_p<- which(rownames(Z.c.p) %in% causal.c)
# subset_matrix <- Z.c.p[is_in_causal_p, ]
# phylum_occurrences <- colSums(subset_matrix)
# causal.p <- names(which(phylum_occurrences > 0))
# 
# #Calculate expected values for Phyla
# ExpectedValues.p <- data.frame()
# phylum.names <- unique(formatted.data$Taxa[which(grepl("phylum",formatted.data$Taxa))])
# for (p in phylum.names) {
#   causal <- length(which(paste0("class",which(Z.c.p[,p]==1)) %in% causal.c))
#   total <- length(which(Z.c.p[,p] == 1))
#   Expected.p <- data.frame(Taxa=p,ExpectedLogOdds=log(1.5) *(causal / total),Index=as.numeric(str_remove(p,"phylum")))
#   ExpectedValues.p <- rbind(ExpectedValues.p,Expected.p)
# }
# 
# #Calculate expected values for Class
# ExpectedValues.c <- data.frame()
# class.names <- unique(formatted.data$Taxa[which(grepl("class",formatted.data$Taxa))])
# for (c in class.names) {
#   causal <- length(which(paste0("order",which(Z.o.c[,c]==1)) %in% causal.o))
#   total <- length(which(Z.o.c[,c] == 1))
#   Expected.c <- data.frame(Taxa=c,ExpectedLogOdds=log(1.5) *(causal / total),Index=as.numeric(str_remove(c,"class")))
#   ExpectedValues.c <- rbind(ExpectedValues.c,Expected.c)
# }
# 
# #Calculate expected values for Order
# ExpectedValues.o <- data.frame()
# order.names <- unique(formatted.data$Taxa[which(grepl("order",formatted.data$Taxa))])
# for (o in order.names) {
#   causal <- length(which(paste0("family",which(Z.f.o[,o]==1)) %in% causal.f))
#   total <- length(which(Z.f.o[,o] == 1))
#   Expected.o <- data.frame(Taxa=o,ExpectedLogOdds= log(1.5) *(causal / total),Index=as.numeric(str_remove(o,"order")))
#   ExpectedValues.o <- rbind(ExpectedValues.o,Expected.o)
# }
# #
# #Calculate expected values for Family
# ExpectedValues.f <- data.frame()
# family.names <- unique(formatted.data$Taxa[which(grepl("family",formatted.data$Taxa))])
# for (f in family.names) {
#   causal <- length(which(paste0("genus",which(Z.g.f[,f]==1)) %in% causal.g))
#   total <- length(which(Z.g.f[,f] == 1))
#   Expected.f <- data.frame(Taxa=f,ExpectedLogOdds=log(1.5) *(causal / total),Index=as.numeric(str_remove(f,"family")))
#   ExpectedValues.f <- rbind(ExpectedValues.f,Expected.f)
# }
# 
# #Calculate expected values for Genus
# ExpectedValues.g <- data.frame()
# genus.names <- unique(formatted.data$Taxa[which(grepl("genus",formatted.data$Taxa))])
# for (g in genus.names) {
#   causal <- length(which(paste0("species",which(Z.s.g[,g]==1)) %in% causal.s))
#   total <- length(which(Z.s.g[,g] == 1))
#   Expected.g <- data.frame(Taxa=g,ExpectedLogOdds=log(1.5) *(causal / total),Index=as.numeric(str_remove(g,"genus")))
#   ExpectedValues.g <- rbind(ExpectedValues.g,Expected.g)
# }
# 
# #Calculate expected values for Species
# ExpectedValues.s <- data.frame()
# species.names <- unique(formatted.data$Taxa[which(grepl("species",formatted.data$Taxa))])
# for (s in species.names) {
#   # causal <- length(which(paste0("species",which(Z.s.g[,s]==1)) %in% causal.s))
#   # total <- length(which(Z.s.g[,g] == 1))
#   Expected.s <- data.frame(Taxa=s,ExpectedLogOdds=ifelse(s %in% causal.s,log(1.5),0),Index=as.numeric(str_remove(s,"species")))
#   ExpectedValues.s <- rbind(ExpectedValues.s,Expected.s)
# }
# 
# ExpectedValues <- rbind(ExpectedValues.p,ExpectedValues.c)
# ExpectedValues <- rbind(ExpectedValues,ExpectedValues.o)
# ExpectedValues <- rbind(ExpectedValues,ExpectedValues.f)
# ExpectedValues <- rbind(ExpectedValues,ExpectedValues.g)
# ExpectedValues <- rbind(ExpectedValues,ExpectedValues.s)
# formatted.data <- tidylog::left_join(formatted.data,ExpectedValues,by="Taxa")
# 
# #Fix index
# formatted.data <- formatted.data %>%
#   tidylog::mutate(index=case_when(Taxa.Level=="Species" ~ str_remove(Taxa,"species"),
#                                   Taxa.Level=="Genus" ~ str_remove(Taxa,"genus"),
#                                   Taxa.Level=="Family" ~ str_remove(Taxa,"family"),
#                                   Taxa.Level=="Order" ~ str_remove(Taxa,"order"),
#                                   Taxa.Level=="Class" ~ str_remove(Taxa,"class"),
#                                   Taxa.Level=="Phylum" ~ str_remove(Taxa,"phylum")))
# # tidylog::mutate(index=as.numeric(index))
# formatted.data <- formatted.data %>%
#   tidylog::mutate(expectedLogOdds = case_when(Exposure=="Mixture" ~ NA,
#                                               Exposure %in% causal.e ~ ExpectedLogOdds,
#                                               !Exposure %in% causal.e ~ 0))
# 
# 
# # Expected value for mixture =sum of individual effects
# formatted.data <-  formatted.data %>%
#   # filter(Exposure!="Mixture") %>% 
#   tidylog::group_by(Taxa,Model,Replication) %>%
#   tidylog::mutate(mixture_sum = sum(expectedLogOdds,na.rm=T)) %>%
#   tidylog::mutate(expectedLogOdds=ifelse(is.na(expectedLogOdds),mixture_sum,expectedLogOdds)) %>%
#   ungroup() %>% 
#   select(-mixture_sum)
# 
# #Create Causal/Non-causal indicator: If expected value > 0 C, else NC
# formatted.data <- formatted.data %>%
#   mutate(indicator=ifelse(expectedLogOdds>0,"Causal","Non-Causal"))
# 
# causal.taxa <- c(causal.s,causal.g,causal.f,causal.o,causal.c,causal.p)
# 
# return(formatted.data)
# 
# }
# 
# formatted.scenario2 <- format.fxn(results.scenario2)
# formatted.data <- rbind(formatted.scenario2)
# rm(formatted.scenario2,results.scenario2)