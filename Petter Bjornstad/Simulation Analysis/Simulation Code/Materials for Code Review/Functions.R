#Format Simulation Results Function ----
format.fxn <- function(results.data){
  results.data <- results.data %>% 
    dplyr::mutate(component=case_when(component=="Count model coefficients" ~ "Means",
                                      component=="Zero-inflation model coefficients" ~ "Probability",
                                      component=="Means" ~ "Means",
                                      component=="Probability" ~ "Probability")) %>% 
    mutate(sig=case_when(pdir > 0.975 ~ "*",
                         pdir <= 0.975 ~ "N.S.",
                         is.na(pdir) ~ sig))
  
  zing.data <- results.data %>%
    filter(model=="ZINB" | model =="Poisson") %>%
    filter(component=="Means") %>%
    mutate(pval = as.numeric(pval)) %>%
    dplyr::group_by(domain) %>%
    mutate(fdr = p.adjust(pval, method = "fdr")) %>%
    ungroup() %>%
    mutate(fdr.sig=ifelse(fdr<0.05,"*","N.S.")) %>%
    mutate(model="ZING") %>% 
    mutate(lcl = estimate - (1.96 * SE),ucl = estimate + (1.96 * SE)
    ) %>% 
    dplyr::select(c("taxa_full","exposure","estimate","lcl","ucl","component","model","domain","Scenario","sig","fdr.sig","iteration"))
  
  # zing.data.prob <- results.data %>%
  #   filter(model=="ZINB" | model =="Poisson") %>%
  #   filter(component=="Probability") %>%
  #   mutate(pval = as.numeric(pval)) %>%
  #   dplyr::group_by(domain) %>%
  #   mutate(fdr = p.adjust(pval, method = "fdr")) %>%
  #   ungroup() %>%
  #   mutate(fdr.sig=ifelse(fdr<0.05,"*","N.S.")) %>%
  #   mutate(model="ZING") 
  
  bayesian <- results.data %>%
    filter(component=="Means") %>%
    filter(model=="BaHZING" | model =="RBaHZING") %>%
    mutate(taxa_full=str_replace_all(taxa_full,"k__species","species")) %>% 
    mutate(taxa_full=str_replace_all(taxa_full,"k__X","species")) %>% 
    # mutate(fdr=p) %>% 
    # mutate(fdr.sig=ifelse(fdr<0.05,"*","N.S.")) 
    mutate(fdr.sig=NA) %>% 
    dplyr::rename(lcl=bci_lcl,
                  ucl=bci_ucl) %>% 
    dplyr::select(c("taxa_full","exposure","estimate","lcl","ucl","component","model","domain","Scenario","sig","fdr.sig","iteration"))
  
  # formatted.data.prob <- results.data %>%
  #   filter(component=="Probability") %>%
  #   filter(model=="BaHZING" | model =="RBaHZING") %>%
  #   mutate(taxa_full=str_replace_all(taxa_full,"k__species","species")) %>% 
  #   mutate(fdr=pval) %>% 
  #   mutate(fdr.sig=ifelse(fdr<0.05,"*","N.S.")) 
  
  formatted.data <- rbind(zing.data,  bayesian)
  # formatted.data <- formatted.data %>% 
  # dplyr::rename(sig=fdr.sig)
  # formatted.data.prob <- rbind(zing.data.prob,formatted.data.prob)
  # rm(zing.data,zing.data.prob)
  rm(zing.data,bayesian)
  
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
  
  # # Define a function to assign replication numbers
  # assign_replication_numbers <- function(group) {
  #   group$Replication <- seq_len(nrow(group))
  #   return(group)
  # }
  
  # # Apply the function to each group
  # formatted.data <- formatted.data %>%
  #   tidylog::group_by(Taxa, Exposure, Model,iteration) %>%
  #   group_modify(~ assign_replication_numbers(.x))
  
  #Group number for plotting
  formatted.data <- formatted.data %>%
    tidylog::group_by(taxa_full,exposure,iteration) %>%
    tidylog::mutate(group_number = cur_group_id()) %>%
    ungroup()
  
  # formatted.data.prob <- formatted.data.prob %>%
  #   tidylog::group_by(taxa_full,exposure,iteration) %>%
  #   tidylog::mutate(group_number = cur_group_id()) %>%
  #   ungroup()
  
  sim.par <- sim.par.all[unique(formatted.data$Scenario),]
  par <- sim.par$P.s.scenario
  
  #Causal Exposures
  if (sim.par$P.e.causal==0){
    causal.e <- NULL
  } else {
    causal.e <- paste0("X.",rep(1:sim.par$P.e.causal))
  }
  
  #Species level indicator
  if(par==1){
    causal.s <- NULL
  }
  if(par!=1){
    causal.s <- paste0("species",simulation.parameters[[par]])
  }
  
  #Genus level indicator
  is_in_causal_s <- which(rownames(Z.s.g) %in% causal.s)
  subset_matrix <- Z.s.g[is_in_causal_s, ]
  genus_occurrences <- colSums(subset_matrix)
  causal.g <- names(which(genus_occurrences > 0))
  
  #Family level indicator
  is_in_causal_g <- which(rownames(Z.g.f) %in% causal.g)
  subset_matrix <- Z.g.f[is_in_causal_g, ]
  family_occurrences <- colSums(subset_matrix)
  causal.f <- names(which(family_occurrences > 0))
  
  #Order level indicator
  is_in_causal_f <- which(rownames(Z.f.o) %in% causal.f)
  subset_matrix <- Z.f.o[is_in_causal_f, ]
  order_occurrences <- colSums(subset_matrix)
  causal.o <- names(which(order_occurrences > 0))
  
  # Class level indicator
  is_in_causal_c <- which(rownames(Z.o.c) %in% causal.o)
  subset_matrix <- Z.o.c[is_in_causal_c, ]
  class_occurrences <- colSums(subset_matrix)
  causal.c <- names(which(class_occurrences > 0))
  
  #Phylum level indicator
  is_in_causal_p<- which(rownames(Z.c.p) %in% causal.c)
  subset_matrix <- Z.c.p[is_in_causal_p, ]
  phylum_occurrences <- colSums(subset_matrix)
  causal.p <- names(which(phylum_occurrences > 0))
  
  #Set the intended OR for the scenario
  setOR <- sim.par$OR.exposure
  null <- 1
  
  #Calculate expected values for Phyla
  ExpectedValues.p <- data.frame()
  phylum.names <- unique(formatted.data$taxa_full[which(grepl("phylum",formatted.data$taxa_full))])
  for (p in phylum.names) {
    causal <- length(which(paste0("class",which(Z.c.p[,p]==1)) %in% causal.c))
    total <- length(which(Z.c.p[,p] == 1))
    Expected.p <- data.frame(taxa_full=p,ExpectedLogOdds=log(setOR) *(causal / total),Index=as.numeric(str_remove(p,"phylum")))
    ExpectedValues.p <- rbind(ExpectedValues.p,Expected.p)
  }
  
  #Calculate expected values for Class
  ExpectedValues.c <- data.frame()
  class.names <- unique(formatted.data$taxa_full[which(grepl("class",formatted.data$taxa_full))])
  for (c in class.names) {
    causal <- length(which(paste0("order",which(Z.o.c[,c]==1)) %in% causal.o))
    total <- length(which(Z.o.c[,c] == 1))
    Expected.c <- data.frame(taxa_full=c,ExpectedLogOdds=log(setOR) *(causal / total),Index=as.numeric(str_remove(c,"class")))
    ExpectedValues.c <- rbind(ExpectedValues.c,Expected.c)
  }
  
  #Calculate expected values for Order
  ExpectedValues.o <- data.frame()
  order.names <- unique(formatted.data$taxa_full[which(grepl("order",formatted.data$taxa_full))])
  for (o in order.names) {
    causal <- length(which(paste0("family",which(Z.f.o[,o]==1)) %in% causal.f))
    total <- length(which(Z.f.o[,o] == 1))
    Expected.o <- data.frame(taxa_full=o,ExpectedLogOdds= log(setOR) *(causal / total),Index=as.numeric(str_remove(o,"order")))
    ExpectedValues.o <- rbind(ExpectedValues.o,Expected.o)
  }
  
  #Calculate expected values for Family
  ExpectedValues.f <- data.frame()
  family.names <- unique(formatted.data$taxa_full[which(grepl("family",formatted.data$taxa_full))])
  for (f in family.names) {
    causal <- length(which(paste0("genus",which(Z.g.f[,f]==1)) %in% causal.g))
    total <- length(which(Z.g.f[,f] == 1))
    Expected.f <- data.frame(taxa_full=f,ExpectedLogOdds=log(setOR) *(causal / total),Index=as.numeric(str_remove(f,"family")))
    ExpectedValues.f <- rbind(ExpectedValues.f,Expected.f)
  }
  
  #Calculate expected values for Genus
  ExpectedValues.g <- data.frame()
  genus.names <- unique(formatted.data$taxa_full[which(grepl("genus",formatted.data$taxa_full))])
  for (g in genus.names) {
    causal <- length(which(paste0("species",which(Z.s.g[,g]==1)) %in% causal.s))
    total <- length(which(Z.s.g[,g] == 1))
    Expected.g <- data.frame(taxa_full=g,ExpectedLogOdds=log(setOR) *(causal / total),Index=as.numeric(str_remove(g,"genus")))
    ExpectedValues.g <- rbind(ExpectedValues.g,Expected.g)
  }
  
  #Calculate expected values for Species - should this be unaffect by other levels?? Check this!!!
  ExpectedValues.s <- data.frame()
  species.names <- unique(formatted.data$taxa_full[which(grepl("species",formatted.data$taxa_full))])
  for (s in species.names) {
    # causal <- length(which(paste0("species",which(Z.s.g[,s]==1)) %in% causal.s))
    # total <- length(which(Z.s.g[,g] == 1))
    Expected.s <- data.frame(taxa_full=s,ExpectedLogOdds=ifelse(s %in% causal.s,log(setOR),0),Index=as.numeric(str_remove(s,"species")))
    ExpectedValues.s <- rbind(ExpectedValues.s,Expected.s)
  }
  
  ExpectedValues <- rbind(ExpectedValues.p,ExpectedValues.c)
  ExpectedValues <- rbind(ExpectedValues,ExpectedValues.o)
  ExpectedValues <- rbind(ExpectedValues,ExpectedValues.f)
  ExpectedValues <- rbind(ExpectedValues,ExpectedValues.g)
  ExpectedValues <- rbind(ExpectedValues,ExpectedValues.s)
  formatted.data <- tidylog::left_join(formatted.data,ExpectedValues,by="taxa_full")
  
  #Now make only the causal exposures causal
  formatted.data <- formatted.data %>%
    tidylog::mutate(expectedLogOdds = case_when(exposure=="Mixture" ~ NA,
                                                # exposure!="Mixture" ~ ExpectedLogOdds))
                                                exposure %in% causal.e ~ ExpectedLogOdds,
                                                !exposure %in% causal.e ~ log(null)))
  #Remove the intermediate ExpectedLogOdds variable
  formatted.data <- formatted.data %>% 
    dplyr::select(-c("ExpectedLogOdds"))
  
  
  # Expected value for mixture =sum of individual effects
  formatted.data <-  formatted.data %>%
    # filter(Exposure!="Mixture") %>% 
    tidylog::group_by(taxa_full,model,iteration) %>%
    tidylog::mutate(mixture_sum = sum(expectedLogOdds,na.rm=T)) %>%
    tidylog::mutate(expectedLogOdds=ifelse(is.na(expectedLogOdds),mixture_sum,expectedLogOdds)) %>%
    ungroup() %>% 
    select(-mixture_sum)
  
  #Create Causal/Non-causal indicator: If expected value > 0 C, else NC
  formatted.data <- formatted.data %>%
    mutate(indicator=ifelse(expectedLogOdds>0,"Associated","Not Associated"))
  
  associated.taxa <- c(causal.s,causal.g,causal.f,causal.o,causal.c,causal.p)
  
  return(formatted.data)
  
}

#Plot Function ----
create_estimate_plots <- function(data, 
                                  analysis_type = "mixture",
                                  output_dir = "/home/hhampson/Results/Microbiome Results",
                                  save_plots = FALSE,
                                  n_sample_iterations = 25,
                                  seed = 123) {
  
  require(tidyverse)
  require(ggplot2)
  
  # Validate analysis_type
  if (!analysis_type %in% c("mixture", "individual")) {
    stop("analysis_type must be either 'mixture' or 'individual'")
  }
  
  # Set seed and sample iterations
  set.seed(seed)
  sampled_iterations <- sample(unique(data$iteration), n_sample_iterations)
  data_sample <- data %>% filter(iteration %in% sampled_iterations)
  
  # Filter for analysis type
  if (analysis_type == "mixture") {
    plot_data <- data_sample %>% filter(exposure == "Mixture")
    file_prefix <- "Mixture_Plot"
  } else {
    plot_data <- data_sample %>% filter(exposure != "Mixture")
    file_prefix <- "Individual_Plot"
  }
  
  # Format factors
  plot_data <- plot_data %>%
    mutate(
      model = factor(model, levels = c("BaHZING", "RBaHZING", "ZING")),
      Scenario = factor(paste0("Scenario ", Scenario), 
                        levels = paste0("Scenario ", 1:9)),
      domain = factor(domain, 
                      levels = c("Species", "Genus", "Family", "Order", "Class", "Phylum"))
    )
  
  # Split by association status
  NA_data <- plot_data %>% filter(indicator == "Not Associated")
  A_data <- plot_data %>% filter(indicator == "Associated")
  
  # Base theme
  base_theme <- theme_bw() +
    theme(
      text = element_text(family = "Times", size = 20, color = "black"),
      axis.text = element_text(family = "Times", size = 15),
      axis.title.x = element_blank(),
      axis.title.y = element_text(family = "Times", size = 20, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      strip.text.x = element_text(size = 15, face = "bold"),
      strip.text.y = element_text(size = 15, face = "bold")
    )
  
  # Helper function for error bars
  mean_sd <- function(x) {
    mean_val <- mean(x)
    sd_val <- sd(x)
    data.frame(y = mean_val, ymin = mean_val - sd_val, ymax = mean_val + sd_val)
  }
  
  # Color by association status
  combined_data <- bind_rows(
    NA_data %>% mutate(status = "True Null"),
    A_data %>% mutate(status = "True Association")
  )
  
  plot <- ggplot(combined_data, aes(model, estimate, color = status)) +
    geom_line(aes(group = group_number), alpha = 0.15, linewidth = 0.3) +
    geom_point(aes(group = group_number), alpha = 0.3, size = 0.8) +
    geom_errorbar(stat = "summary", fun.data = mean_sd,
                  width = 0.2, linewidth = 0.8, alpha = 0.8) +
    geom_point(stat = "summary", fun = mean, size = 3, 
               shape = 21, fill = "white", stroke = 1.5) +
    facet_grid(Scenario ~ domain, scales = "free_y") +
    base_theme +
    theme(legend.position = "right",
          legend.title = element_text(size = 16, face = "bold"),
          legend.text = element_text(size = 14),
          title = element_text(size=20,face="bold")) +
    scale_color_manual(
      name = "Simulated Association",
      values = c("True Association" = "#540b0e", "True Null" = "gray")) +
    labs(y = "Model Estimate",
         title = paste0(str_to_title(analysis_type)))
  
  # Save plot if requested
  if (save_plots) {
    png(file.path(output_dir, paste0(file_prefix, ".png")),
        res = 300, height = 4000, width = 4000)
    print(plot)
    dev.off()
  }
  
  # Return plot
  return(plot)
}

# Sensitivity Measures ----
# Function to calculate sensitivity measures and create heatmaps
calculate_sensitivity_measures <- function(data,
                                           analysis_type = "individual",
                                           output_dir = "/home/hhampson/Results/Microbiome Results",
                                           save_plots = FALSE,
                                           save_csv = FALSE) {
  
  require(tidyverse)
  require(reshape2)
  
  # Filter for analysis type
  if (analysis_type == "individual") {
    filtered_data <- data %>% filter(exposure != "Mixture")
    file_prefix <- "Individual"
  } else if (analysis_type == "mixture") {
    filtered_data <- data %>% filter(exposure == "Mixture")
    file_prefix <- "Mixture"
  } else {
    stop("analysis_type must be either 'individual' or 'mixture'")
  }
  
  # Helper function to calculate metrics for one model
  calc_metrics <- function(model_data, model_name, sig_var = "sig") {
    model_data %>%
      mutate(PosNeg = case_when(
        indicator == "Associated" & !!sym(sig_var) == "*" ~ "True Positive",
        indicator == "Associated" & !!sym(sig_var) != "*" ~ "False Negative",
        indicator == "Not Associated" & !!sym(sig_var) == "*" ~ "False Positive",
        indicator == "Not Associated" & !!sym(sig_var) != "*" ~ "True Negative")) %>%
      group_by(domain, Scenario) %>%
      summarise(
        Specificity = length(which(PosNeg == "True Negative")) / 
          (length(which(PosNeg == "True Negative")) + length(which(PosNeg == "False Positive"))),
        FDR = length(which(PosNeg == "False Positive")) / 
          (length(which(PosNeg == "True Positive")) + length(which(PosNeg == "False Positive"))),
        .groups = 'drop'
      ) %>%
      mutate(Model = model_name)
  }
  
  # Calculate metrics for each model
  BaH_ZING <- filtered_data %>%
    filter(model == "BaHZING") %>%
    calc_metrics("BaHZING")
  
  RBaH_ZING <- filtered_data %>%
    filter(model == "RBaHZING") %>%
    calc_metrics("RBaHZING")
  
  ZING <- filtered_data %>%
    filter(model == "ZING") %>%
    calc_metrics("ZING")
  
  ZING_fdr <- filtered_data %>%
    filter(model == "ZING") %>%
    calc_metrics("Adj.ZING", sig_var = "fdr.sig")
  
  # Combine all results
  SensSpec <- rbind(BaH_ZING, RBaH_ZING, ZING, ZING_fdr)
  
  # Reshape to long format
  melted_data <- melt(SensSpec, id.vars = c('Model', 'domain', 'Scenario'))
  
  # Factor ordering
  model_order <- c('BaHZING', 'RBaHZING', 'ZING', 'Adj.ZING')
  melted_data$Model <- factor(melted_data$Model, levels = model_order)
  
  taxa_order <- c("Species", "Genus", "Family", "Order", "Class", "Phylum")
  melted_data$domain <- factor(melted_data$domain, levels = taxa_order)
  
  # Add scenario labels
  melted_data <- melted_data %>%
    mutate(Scenario = paste0("Scenario ", Scenario))
  
  scenario_order <- paste0("Scenario ", 1:9)
  melted_data$Scenario <- factor(melted_data$Scenario, levels = scenario_order)
  
  # Create normalized value and mark NAs
  melted_data <- melted_data %>%
    mutate(
      is_na = is.na(value) | is.nan(value) | is.infinite(value),
      normalized_value = case_when(
        is_na ~ NA_real_,
        variable == "FDR" ~ 1 - value,
        TRUE ~ value
      )
    )
  
  # Create wide format for CSV
  formatted_wide <- melted_data %>%
    select(Model, domain, Scenario, variable, value) %>%
    pivot_wider(names_from = domain, values_from = value) %>%
    select(Model, Scenario, variable, Species, Genus, Family, Order, Class, Phylum)
  
  # Save CSV if requested
  if (save_csv) {
    write.csv(formatted_wide, 
              file.path(output_dir, paste0("Sensitivity_Measures_", file_prefix, ".csv")),
              row.names = FALSE)
  }
  
  # Base theme for heatmaps
  base_heatmap_theme <- theme_bw() +
    theme(
      text = element_text(family = "Times", size = 12),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      axis.text.y = element_text(size = 10),
      axis.title = element_text(size = 12, face = "bold"),
      strip.text.x = element_text(size = 10, face = "bold"),
      strip.text.y = element_text(size = 10, face = "bold"),
      legend.position = "right",
      legend.title = element_text(size = 11, face = "bold"),
      legend.text = element_text(size = 10))
  
  # Create heatmap for Specificity (HIGH IS GOOD)
  spec_data <- melted_data %>%
    filter(variable == "Specificity")
  
  heatmap_specificity <- ggplot(spec_data, 
                                aes(x = domain, y = Model, fill = value)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_tile(data = filter(spec_data, is_na), 
              aes(x = domain, y = Model), 
              fill = "#f0f0f0", color = "white", linewidth = 0.5) +
    facet_grid(. ~ Scenario) +
    base_heatmap_theme +
    scale_fill_gradientn(
      colours = c("#9b2226", "#ee9b00", "#e9d8a6", "#94d2bd", "#005f73"),
      values = scales::rescale(c(0, 0.25, 0.5, 0.75, 1)),
      name = "Specificity\n(Higher = Better)",
      limits = c(0, 1),
      na.value = "#f0f0f0") +
    labs(x = "Taxonomic Level", y = "Model",
         title = "Specificity")
  
  # Create heatmap for FDR (LOW IS GOOD)
  fdr_data <- melted_data %>%
    filter(variable == "FDR")
  
  heatmap_fdr <- ggplot(fdr_data, 
                        aes(x = domain, y = Model, fill = normalized_value)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_tile(data = filter(fdr_data, is_na), 
              aes(x = domain, y = Model), 
              fill = "#f0f0f0", color = "white", linewidth = 0.5) +
    facet_grid(. ~ Scenario) +
    base_heatmap_theme +
    scale_fill_gradientn(
      colours = c("#005f73", "#94d2bd", "#e9d8a6", "#ee9b00", "#9b2226"),
      values = scales::rescale(c(0, 0.25, 0.5, 0.75, 1)),
      name = "FDR\n(Lower = Better)",
      limits = c(0, 1),
      na.value = "#f0f0f0")+
    labs(x = "Taxonomic Level", y = "Model",
         title = "False Discovery Rate")
  
  # Save plots if requested
  if (save_plots) {
    png(file.path(output_dir, paste0(file_prefix, "_Specificity_Heatmap.png")),
        res = 300, width = 5000, height = 1500)
    print(heatmap_specificity)
    dev.off()
    
    png(file.path(output_dir, paste0(file_prefix, "_FDR_Heatmap.png")),
        res = 300, width = 5000, height = 1500)
    print(heatmap_fdr)
    dev.off()
  }
  
  # Return results
  list(
    long_format = melted_data,
    wide_format = formatted_wide,
    heatmap_specificity = heatmap_specificity,
    heatmap_fdr = heatmap_fdr
  )
}

#Bias ----
calculate_bias_species <- function(data,
                                   analysis_type = "individual",
                                   output_dir = "/home/hhampson/Results/Microbiome Results",
                                   save_plots = FALSE,
                                   save_csv = FALSE) {
  
  require(tidyverse)
  
  # Filter for analysis type and species level only
  if (analysis_type == "individual") {
    filtered_data <- data %>% 
      filter(exposure != "Mixture", domain == "Species")
    file_prefix <- "Individual"
  } else if (analysis_type == "mixture") {
    filtered_data <- data %>% 
      filter(exposure == "Mixture", domain == "Species")
    file_prefix <- "Mixture"
  } else {
    stop("analysis_type must be either 'individual' or 'mixture'")
  }
  
  # Calculate OVERALL BIAS (all taxa combined)
  overall_bias <- filtered_data %>%
    group_by(model, Scenario) %>%
    summarise(
      Overall_Bias = mean(estimate - expectedLogOdds, na.rm = TRUE),
      .groups = 'drop'
    )
  
  # Calculate BIAS FOR ASSOCIATED TAXA (systematic over/underestimation)
  associated_bias <- filtered_data %>%
    filter(indicator == "Associated") %>%
    group_by(model, Scenario) %>%
    summarise(
      Associated_Bias = mean(estimate - expectedLogOdds, na.rm = TRUE),
      .groups = 'drop'
    )
  
  # Calculate ABSOLUTE BIAS FOR NOT ASSOCIATED TAXA (magnitude of false signals)
  not_associated_bias <- filtered_data %>%
    filter(indicator == "Not Associated") %>%
    group_by(model, Scenario) %>%
    summarise(
      NotAssociated_AbsBias = mean(abs(estimate - expectedLogOdds), na.rm = TRUE),
      NotAssociated_MSE = mean((estimate - expectedLogOdds)^2, na.rm = TRUE),
      .groups = 'drop'
    )
  
  # Combine all bias metrics
  bias_results <- overall_bias %>%
    left_join(associated_bias, by = c("model", "Scenario")) %>%
    left_join(not_associated_bias, by = c("model", "Scenario"))
  
  # Factor ordering
  model_order <- c('BaHZING', 'RBaHZING', 'ZING')
  bias_results$model <- factor(bias_results$model, levels = model_order)
  
  # Add scenario labels
  bias_results <- bias_results %>%
    mutate(Scenario = paste0("Scenario ", Scenario))
  
  scenario_order <- paste0("Scenario ", 1:9)
  bias_results$Scenario <- factor(bias_results$Scenario, levels = scenario_order)
  
  # Save CSV if requested
  if (save_csv) {
    write.csv(bias_results, 
              file.path(output_dir, paste0("Bias_Species_", file_prefix, ".csv")),
              row.names = FALSE)
  }
  
  # Base theme for heatmaps
  base_theme <- theme_bw() +
    theme(
      text = element_text(family = "Times", size = 12),
      axis.text.y = element_text(size = 10),
      axis.title = element_text(size = 12, face = "bold"),
      strip.text.x = element_text(size = 10, face = "bold"),
      strip.text.y = element_text(size = 10, face = "bold"),
      legend.position = "right",
      legend.title = element_text(size = 11, face = "bold"),
      legend.text = element_text(size = 10))
  
  # HEATMAP 1: Overall Bias
  max_overall <- max(abs(bias_results$Overall_Bias), na.rm = TRUE)
  
  heatmap_overall <- ggplot(bias_results, 
                            aes(x = Scenario, y = model, fill = Overall_Bias)) +
    geom_tile(color = "white", linewidth = 0.5) +
    base_theme +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10)) +
    scale_fill_gradientn(
      colours = c("#005f73", "#94d2bd", "#e9d8a6", "#ee9b00", "#9b2226"),
      values = scales::rescale(c(-max_overall, -max_overall/2, 0, max_overall/2, max_overall)),
      name = "Overall Bias\n(- = Underestimate\n+ = Overestimate)",
      limits = c(-max_overall, max_overall)) +
    labs(x = "Scenario", y = "Model",
         title = "Overall Bias at Species Level (All Taxa)")
  
  # HEATMAP 2: Associated Taxa Bias
  max_assoc <- max(abs(bias_results$Associated_Bias), na.rm = TRUE)
  
  heatmap_associated <- ggplot(bias_results, 
                               aes(x = Scenario, y = model, fill = Associated_Bias)) +
    geom_tile(color = "white", linewidth = 0.5) +
    base_theme +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10)) +
    scale_fill_gradientn(
      colours = c("#005f73", "#94d2bd", "#e9d8a6", "#ee9b00", "#9b2226"),
      values = scales::rescale(c(-max_assoc, -max_assoc/2, 0, max_assoc/2, max_assoc)),
      name = "Associated Bias\n(- = Underestimate\n+ = Overestimate)",
      limits = c(-max_assoc, max_assoc)) +
    labs(x = "Scenario", y = "Model",
         title = "Bias for Associated Taxa (Systematic Error in Effect Estimates)")
  
  # HEATMAP 3: Not Associated Taxa Absolute Bias
  heatmap_not_associated <- ggplot(bias_results, 
                                   aes(x = Scenario, y = model, fill = NotAssociated_AbsBias)) +
    geom_tile(color = "white", linewidth = 0.5) +
    base_theme +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10)) +
    scale_fill_gradientn(
      colours = c("#005f73", "#94d2bd", "#e9d8a6", "#ee9b00", "#9b2226"),
      values = scales::rescale(c(0, 0.25, 0.5, 0.75, 1)),
      name = "Absolute Bias\n(Lower = Better)",
      limits = c(0, max(bias_results$NotAssociated_AbsBias, na.rm = TRUE))) +
    labs(x = "Scenario", y = "Model",
         title = "Absolute Bias for Not Associated Taxa (Magnitude of False Signals)")
  
  # Save plots if requested
  if (save_plots) {
    png(file.path(output_dir, paste0(file_prefix, "_Bias_Overall_Heatmap.png")),
        res = 300, width = 4000, height = 1500)
    print(heatmap_overall)
    dev.off()
    
    png(file.path(output_dir, paste0(file_prefix, "_Bias_Associated_Heatmap.png")),
        res = 300, width = 4000, height = 1500)
    print(heatmap_associated)
    dev.off()
    
    png(file.path(output_dir, paste0(file_prefix, "_Bias_NotAssociated_Heatmap.png")),
        res = 300, width = 4000, height = 1500)
    print(heatmap_not_associated)
    dev.off()
  }
  
  # Return results
  list(
    bias_summary = bias_results,
    heatmap_overall = heatmap_overall,
    heatmap_associated = heatmap_associated,
    heatmap_not_associated = heatmap_not_associated
  )
}
