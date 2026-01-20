plot_function <- function(data, 
                          analysis_type = c("mixture", "individual"),
                          output_dir = "/home/hhampson/Results/Microbiome Results",
                          save_plots = FALSE,
                          n_sample_iterations = 25,
                          seed = 123) {
  
  analysis_type <- match.arg(analysis_type)
  
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
  
  # OPTION 1: Gradient color by estimate
  plot1 <- ggplot() +
    geom_line(data = NA_data, alpha = 0.5, 
              aes(model, estimate, group = group_number, color = estimate)) +
    geom_point(data = NA_data, aes(model, estimate, color = estimate), alpha = 0.5) +
    geom_line(data = A_data, alpha = 0.5, 
              aes(model, estimate, group = group_number, color = estimate)) +
    geom_point(data = A_data, aes(model, estimate, color = estimate), alpha = 0.5) +
    geom_errorbar(data = NA_data, aes(model, estimate, group = model),
                  stat = "summary", fun.data = mean_sd, width = 0.1, linewidth = 0.5, color = "black") +
    geom_errorbar(data = A_data, aes(model, estimate, group = model),
                  stat = "summary", fun.data = mean_sd, width = 0.1, linewidth = 0.5, color = "black") +
    geom_point(data = NA_data, aes(model, estimate, group = model),
               stat = "summary", fun = mean, size = 2, shape = 16, color = "black") +
    geom_point(data = A_data, aes(model, estimate, group = model),
               stat = "summary", fun = mean, size = 2, shape = 16, color = "black") +
    facet_grid(Scenario ~ domain, scales = "free_y") +
    base_theme +
    theme(legend.position = "left") +
    scale_colour_gradientn(
      colours = c("darkblue", "#2166ac", "#f7f7f7", "#b2182b", "darkred"),
      values = scales::rescale(c(-1, -0.5, 0, 0.5, 1)),
      name = "Estimate") +
    labs(y = "Estimate")
  
  # OPTION 2: Color by association status
  combined_data <- bind_rows(
    NA_data %>% mutate(status = "Not Associated"),
    A_data %>% mutate(status = "Associated")
  )
  
  plot2 <- ggplot(combined_data, aes(model, estimate, color = status)) +
    geom_line(aes(group = group_number), alpha = 0.4) +
    geom_errorbar(stat = "summary", fun.data = mean_sd,
                  width = 0.2, linewidth = 0.8, alpha = 0.8) +
    geom_point(stat = "summary", fun = mean, size = 3, 
               shape = 21, fill = "white", stroke = 1.5) +
    facet_grid(Scenario ~ domain, scales = "free_y") +
    base_theme +
    theme(legend.position = "right",
          legend.title = element_text(size = 16, face = "bold"),
          legend.text = element_text(size = 14)) +
    scale_color_manual(
      name = "Taxa Status",
      values = c("Associated" = "firebrick4", "Not Associated" = "gray60")) +
    labs(y = "Estimate")
  
  # OPTION 3: Categorical effect sizes
  NA_cat <- NA_data %>%
    mutate(effect_cat = factor(
      case_when(
        estimate >= 0.5 ~ "Strong Positive",
        estimate >= 0.25 ~ "Moderate Positive",
        estimate >= 0 ~ "Weak Positive",
        estimate >= -0.25 ~ "Weak Negative",
        TRUE ~ "Negative"),
      levels = c("Strong Positive", "Moderate Positive", "Weak Positive", 
                 "Weak Negative", "Negative")))
  
  A_cat <- A_data %>%
    mutate(effect_cat = factor(
      case_when(
        estimate >= 0.5 ~ "Strong Positive",
        estimate >= 0.25 ~ "Moderate Positive",
        estimate >= 0 ~ "Weak Positive",
        estimate >= -0.25 ~ "Weak Negative",
        TRUE ~ "Negative"),
      levels = c("Strong Positive", "Moderate Positive", "Weak Positive", 
                 "Weak Negative", "Negative")))
  
  plot3 <- ggplot() +
    geom_line(data = NA_cat, alpha = 0.5, 
              aes(model, estimate, group = group_number, color = effect_cat)) +
    geom_line(data = A_cat, alpha = 0.5, 
              aes(model, estimate, group = group_number, color = effect_cat)) +
    geom_errorbar(data = NA_cat, aes(model, estimate, group = model),
                  stat = "summary", fun.data = mean_sd, width = 0.1, linewidth = 0.5, color = "black") +
    geom_errorbar(data = A_cat, aes(model, estimate, group = model),
                  stat = "summary", fun.data = mean_sd, width = 0.1, linewidth = 0.5, color = "black") +
    geom_point(data = NA_cat, aes(model, estimate, group = model),
               stat = "summary", fun = mean, size = 2, shape = 16, color = "black") +
    geom_point(data = A_cat, aes(model, estimate, group = model),
               stat = "summary", fun = mean, size = 2, shape = 16, color = "black") +
    facet_grid(Scenario ~ domain, scales = "free_y") +
    base_theme +
    theme(legend.position = "right") +
    scale_color_manual(
      name = "Effect Size",
      values = c("Strong Positive" = "#660000", "Moderate Positive" = "#cc0000",
                 "Weak Positive" = "#e69f00", "Weak Negative" = "#cccccc",
                 "Negative" = "#9a9a9a")) +
    labs(y = "Estimate")
  
  # Save plots if requested
  if (save_plots) {
    png(file.path(output_dir, paste0(file_prefix, "_Option1.png")),
        res = 300, height = 4000, width = 4000)
    print(plot1)
    dev.off()
    
    png(file.path(output_dir, paste0(file_prefix, "_Option2.png")),
        res = 300, height = 4000, width = 4000)
    print(plot2)
    dev.off()
    
    png(file.path(output_dir, paste0(file_prefix, "_Option3.png")),
        res = 300, height = 4000, width = 4000)
    print(plot3)
    dev.off()
  }
  
  # Return list of plots
  list(option1 = plot1, option2 = plot2, option3 = plot3)
}

# Helper function for error bars
mean_sd <- function(x) {
  mean_val <- mean(x)
  sd_val <- sd(x)
  data.frame(y = mean_val, ymin = mean_val - sd_val, ymax = mean_val + sd_val)
}
  
  ##############
  
  # Filter your data to sampled iterations
  all_formatted_sample <- all_formatted_results %>% 
    filter(iteration %in% sampled_iterations)
  
  mixture <- all_formatted_sample %>% 
    filter(exposure=="Mixture") 
  
  mixture$model <- factor(mixture$model, levels = c("BaHZING","RBaHZING","ZING"))
  mixture$Scenario <- paste0("Scenario ",mixture$Scenario)
  scenario_order <- paste0("Scenario ",rep(1:9))
  mixture$Scenario <- factor(mixture$Scenario, levels = scenario_order)
  taxa_order <- c("Species","Genus","Family","Order","Class","Phylum")
  mixture$domain <- factor(mixture$domain, levels = taxa_order)
  
  # mixture <- mixture %>% 
  #   filter(Scenario=="Scenario 1")
  NA.mix <- mixture %>%
    filter(indicator=="Not Associated")
  A.mix <- mixture %>%
    filter(indicator=="Associated")
  
  mix.plot <- ggplot()+
    geom_line(data = NA.mix, alpha = 0.5, aes(model,estimate,group = group_number,color=estimate)) +
    geom_point(data=NA.mix,aes(model,estimate,color=estimate),alpha = 0.5) +
    geom_line(data = A.mix, alpha = 0.5, aes(model,estimate,group = group_number,color=estimate)) +
    geom_point(data=A.mix,aes(model,estimate,color=estimate),alpha = 0.5) +
    geom_errorbar(data=NA.mix,
                  aes(model,estimate,group = model),
                  stat = "summary",
                  fun.data = function(x) {
                    mean_val <- mean(x)
                    sd_val <- sd(x)
                    data.frame(y = mean_val, ymin = mean_val - sd_val, ymax = mean_val + sd_val)
                  },
                  width = 0.1,
                  linewidth = 0.5,
                  color = "black") +
    geom_errorbar(data=A.mix,
                  aes(model,estimate,group = model),
                  stat = "summary",
                  fun.data = function(x) {
                    mean_val <- mean(x)
                    sd_val <- sd(x)
                    data.frame(y = mean_val, ymin = mean_val - sd_val, ymax = mean_val + sd_val)
                  },
                  width = 0.1,
                  linewidth = 0.5,
                  color = "black") +
    geom_point(data=NA.mix,
               aes(model,estimate,group = model),
               stat = "summary",
               fun = mean,
               size = 2,
               shape = 16,
               color = "black",
               fill = "black") +
    geom_point(data=A.mix,
               aes(model,estimate,group = model),
               stat = "summary",
               fun = mean,
               size = 2,
               shape = 16,
               color = "black",
               fill = "black") +
    facet_grid(Scenario ~ domain,scales = "free_y")+
    theme_bw() +
    theme(
      text = element_text(family = "Times", size = 20,color="black"),
      axis.text = element_text(family = "Times", size = 15),
      plot.title = element_text(family = "Times", face = "bold", size = 16),
      plot.subtitle = element_text(family = "Times", size = 15),
      axis.title.x = element_blank(),
      axis.title.y = element_text(family = "Times", size = 20,face="bold"),
      axis.text.x=element_text(angle=45,hjust=1),
      strip.text.x = element_text(size=15,face="bold"),
      legend.position="left") +
    scale_colour_gradientn(
      colours = c("darkblue","#2166ac", "#f7f7f7", "#b2182b","darkred"),
      values = scales::rescale(c(-1, -0.5, 0, 0.5, 1)),
      name = "Estimate")
  
  # mix.plot
  png("/home/hhampson/Results/Microbiome Results/Mixture_Plot_Option1.png",res=300,height=4000,width=4000)
  plot(mix.plot)
  dev.off()
  
  mix.plot.edit <- ggplot()+
    geom_line(data = NA.mix, alpha = 0.3, 
              aes(model, estimate, group = group_number, color = "Not Associated")) +
    geom_line(data = A.mix, alpha = 0.5, 
              aes(model, estimate, group = group_number, color = "Associated")) +
    geom_errorbar(data = NA.mix,
                  aes(model, estimate, group = model),
                  stat = "summary",
                  fun.data = function(x) {
                    mean_val <- mean(x)
                    sd_val <- sd(x)
                    data.frame(y = mean_val, ymin = mean_val - sd_val, ymax = mean_val + sd_val)
                  },
                  width = 0.2,
                  linewidth = 0.8,
                  color = "gray30") +
    geom_errorbar(data = A.mix,
                  aes(model, estimate, group = model),
                  stat = "summary",
                  fun.data = function(x) {
                    mean_val <- mean(x)
                    sd_val <- sd(x)
                    data.frame(y = mean_val, ymin = mean_val - sd_val, ymax = mean_val + sd_val)
                  },
                  width = 0.2,
                  linewidth = 0.8,
                  color = "firebrick4") +
    geom_point(data = NA.mix,
               aes(model, estimate, group = model),
               stat = "summary",
               fun = mean,
               size = 3,
               shape = 21,
               color = "gray30",
               fill = "white",
               stroke = 1.5) +
    geom_point(data = A.mix,
               aes(model, estimate, group = model),
               stat = "summary",
               fun = mean,
               size = 3,
               shape = 21,
               color = "firebrick4",
               fill = "white",
               stroke = 1.5) +
    facet_grid(Scenario ~ domain, scales = "free_y") +
    theme_bw() +
    theme(
      text = element_text(family = "Times", size = 20, color = "black"),
      axis.text = element_text(family = "Times", size = 15),
      axis.title.x = element_blank(),
      axis.title.y = element_text(family = "Times", size = 20, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      strip.text.x = element_text(size = 15, face = "bold"),
      strip.text.y = element_text(size = 15, face = "bold")) +
    scale_color_manual(
      name = "Taxa Status",
      values = c("Associated" = "firebrick4", 
                 "Not Associated" = "gray60"),
      breaks = c("Associated", "Not Associated")) +
    theme(
      legend.position = "right",
      legend.title = element_text(size = 16, face = "bold"),
      legend.text = element_text(size = 14)) +
    labs(y = "Estimate")
  
  # mix.plot.edit
  
  # ggsave(plot=mix.plot,fs::path(dir.results,"Figure 3.3 Scenario 15.jpeg"),height=4,width=15)
  # png(fs::path(dir.results,"Scenario 4 Test.jpeg"),res=300,height=1000,width=3000)
  png("/home/hhampson/Results/Microbiome Results/Mixture_Plot_Option2.png",res=300,height=4000,width=4000)
  plot(mix.plot.edit)
  dev.off()
  
  # Add a categorical color variable
  # Add a categorical color variable
  NA.mix.cat <- NA.mix %>%
    mutate(effect_cat = case_when(
      estimate >= 0.5 ~ "Strong Positive",
      estimate >= 0.25 ~ "Moderate Positive",
      estimate >= 0 ~ "Weak Positive",
      estimate >= -0.25 ~ "Weak Negative",
      TRUE ~ "Negative"
    ))
  
  A.mix.cat <- A.mix %>%
    mutate(effect_cat = case_when(
      estimate >= 0.5 ~ "Strong Positive",
      estimate >= 0.25 ~ "Moderate Positive",
      estimate >= 0 ~ "Weak Positive",
      estimate >= -0.25 ~ "Weak Negative",
      TRUE ~ "Negative"
    ))
  
  # Create the plot with discrete colors
  mix.plot.edit3 <- ggplot() +
    geom_line(data = NA.mix.cat, alpha = 0.5, 
              aes(model, estimate, group = group_number, color = effect_cat)) +
    geom_line(data = A.mix.cat, alpha = 0.5, 
              aes(model, estimate, group = group_number, color = effect_cat)) +
    geom_errorbar(data = NA.mix.cat,
                  aes(model, estimate, group = model),
                  stat = "summary",
                  fun.data = function(x) {
                    mean_val <- mean(x)
                    sd_val <- sd(x)
                    data.frame(y = mean_val, ymin = mean_val - sd_val, ymax = mean_val + sd_val)
                  },
                  width = 0.1,
                  linewidth = 0.5,
                  color = "black") +
    geom_errorbar(data = A.mix.cat,
                  aes(model, estimate, group = model),
                  stat = "summary",
                  fun.data = function(x) {
                    mean_val <- mean(x)
                    sd_val <- sd(x)
                    data.frame(y = mean_val, ymin = mean_val - sd_val, ymax = mean_val + sd_val)
                  },
                  width = 0.1,
                  linewidth = 0.5,
                  color = "black") +
    geom_point(data = NA.mix.cat,
               aes(model, estimate, group = model),
               stat = "summary",
               fun = mean,
               size = 2,
               shape = 16,
               color = "black",
               fill = "black") +
    geom_point(data = A.mix.cat,
               aes(model, estimate, group = model),
               stat = "summary",
               fun = mean,
               size = 2,
               shape = 16,
               color = "black",
               fill = "black") +
    facet_grid(Scenario ~ domain, scales = "free_y") +
    theme_bw() +
    theme(
      text = element_text(family = "Times", size = 20, color = "black"),
      axis.text = element_text(family = "Times", size = 15),
      plot.title = element_text(family = "Times", face = "bold", size = 16),
      plot.subtitle = element_text(family = "Times", size = 15),
      axis.title.x = element_blank(),
      axis.title.y = element_text(family = "Times", size = 20, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      strip.text.x = element_text(size = 15, face = "bold"),
      legend.position = "right") +
    scale_color_manual(
      name = "Effect Size",
      values = c("Strong Positive" = "#660000",
                 "Moderate Positive" = "#cc0000",
                 "Weak Positive" = "#e69f00",
                 "Weak Negative" = "#cccccc",
                 "Null" = "#9a9a9a"),
      breaks = c("Strong Positive", "Moderate Positive", "Weak Positive", 
                 "Weak Negative", "Negative")) +
    labs(y = "Estimate")
  
  # mix.plot.edit3
  png("/home/hhampson/Results/Microbiome Results/Mixture_Plot_Option3.png",res=300,height=4000,width=4000)
  plot(mix.plot.edit3)
  dev.off()
  
  ##b. Individual Analysis ----
  # ind <- formatted.data %>%
  #   tidylog::filter(Exposure!="Mixture") 
  ind <- all_formatted_sample %>% 
    filter(exposure!="Mixture") 
  
  ind$model <- factor(ind$model, levels = c("BaHZING","RBaHZING","ZING"))
  ind$Scenario <- paste0("Scenario ",ind$Scenario)
  # scenario_order <- c("Scenario 4")
  ind$Scenario <- factor(ind$Scenario, levels = scenario_order)
  taxa_order <- c("Species","Genus","Family","Order","Class","Phylum")
  ind$domain <- factor(ind$domain, levels = taxa_order)
  
  NA.ind <- ind %>%
    filter(indicator=="Not Associated")
  A.ind <- ind %>%
    filter(indicator=="Associated")
  
  ind.plot <- ggplot() +
    geom_line(data = NA.ind, alpha = 0.5, 
              aes(model, estimate, group = group_number, color = estimate)) +
    geom_line(data = A.ind, alpha = 0.5, 
              aes(model, estimate, group = group_number, color = estimate)) +
    geom_errorbar(data = NA.ind,
                  aes(model, estimate, group = model),
                  stat = "summary",
                  fun.data = function(x) {
                    mean_val <- mean(x)
                    sd_val <- sd(x)
                    data.frame(y = mean_val, ymin = mean_val - sd_val, ymax = mean_val + sd_val)
                  },
                  width = 0.1,
                  linewidth = 0.5,
                  color = "black") +
    geom_errorbar(data = A.ind,
                  aes(model, estimate, group = model),
                  stat = "summary",
                  fun.data = function(x) {
                    mean_val <- mean(x)
                    sd_val <- sd(x)
                    data.frame(y = mean_val, ymin = mean_val - sd_val, ymax = mean_val + sd_val)
                  },
                  width = 0.1,
                  linewidth = 0.5,
                  color = "black") +
    geom_point(data = NA.ind,
               aes(model, estimate, group = model),
               stat = "summary",
               fun = mean,
               size = 2,
               shape = 16,
               color = "black",
               fill = "black") +
    geom_point(data = A.ind,
               aes(model, estimate, group = model),
               stat = "summary",
               fun = mean,
               size = 2,
               shape = 16,
               color = "black",
               fill = "black") +
    facet_grid(Scenario ~ domain, scales = "free_y") +
    theme_bw() +
    theme(
      text = element_text(family = "Times", size = 20, color = "black"),
      axis.text = element_text(family = "Times", size = 15),
      plot.title = element_text(family = "Times", face = "bold", size = 16),
      plot.subtitle = element_text(family = "Times", size = 15),
      axis.title.x = element_blank(),
      axis.title.y = element_text(family = "Times", size = 20, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      strip.text.x = element_text(size = 15, face = "bold"),
      legend.position = "left") +
    scale_colour_gradientn(
      colours = c("darkblue","#2166ac", "#f7f7f7", "#b2182b","darkred"),
      values = scales::rescale(c(-1, -0.5, 0, 0.5, 1)),
      name = "Estimate")+
    labs(y = "Estimate")
  
  # ind.plot
  png("/home/hhampson/Results/Microbiome Results/Individual_Plot_Option1.png",res=300,height=4000,width=4000)
  plot(ind.plot)
  dev.off()
  
  # First, combine the data with status indicator
  NA.ind_legend <- NA.ind %>% mutate(status = "Not Associated")
  A.ind_legend <- A.ind %>% mutate(status = "Associated")
  combined_ind_data <- bind_rows(NA.ind_legend, A.ind_legend)
  
  # Now create the plot with the combined data
  ind.plot.edit <- ggplot(combined_ind_data, aes(model, estimate, color = status)) +
    geom_line(aes(group = group_number), 
              alpha = 0.4) +
    geom_errorbar(stat = "summary",
                  fun.data = function(x) {
                    mean_val <- mean(x)
                    sd_val <- sd(x)
                    data.frame(y = mean_val, ymin = mean_val - sd_val, ymax = mean_val + sd_val)
                  },
                  width = 0.2,
                  linewidth = 0.8,
                  alpha = 0.8,
                  position = position_dodge(width = 0.3)) +
    geom_point(stat = "summary",
               fun = mean,
               size = 3,
               shape = 21,
               fill = "white",
               stroke = 1.5,
               position = position_dodge(width = 0.3)) +
    facet_grid(Scenario ~ domain, scales = "free_y") +
    theme_bw() +
    theme(
      text = element_text(family = "Times", size = 20, color = "black"),
      axis.text = element_text(family = "Times", size = 15),
      axis.title.x = element_blank(),
      axis.title.y = element_text(family = "Times", size = 20, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      strip.text.x = element_text(size = 15, face = "bold"),
      strip.text.y = element_text(size = 15, face = "bold"),
      legend.position = "right",
      legend.title = element_text(size = 16, face = "bold"),
      legend.text = element_text(size = 14)) +
    scale_color_manual(
      name = "Taxa Status",
      values = c("Associated" = "firebrick4", 
                 "Not Associated" = "gray60"),
      breaks = c("Associated", "Not Associated")) +
    guides(color = guide_legend(
      override.aes = list(
        shape = 21,
        size = 3,
        stroke = 1.5,
        fill = "white",
        linewidth = 1
      )
    )) +
    labs(y = "Estimate")
  # ind.plot.edit
  
  png("/home/hhampson/Results/Microbiome Results/Individual_Plot_Option2.png",res=300,height=4000,width=4000)
  plot(ind.plot.edit)
  dev.off()
  
  # Add categorical color variable
  NA.ind.cat <- NA.ind %>%
    mutate(effect_cat = case_when(
      estimate >= 0.5 ~ "Strong Positive",
      estimate >= 0.25 ~ "Moderate Positive",
      estimate >= 0 ~ "Weak Positive",
      estimate >= -0.25 ~ "Weak Negative",
      TRUE ~ "Negative"
    )) %>%
    mutate(effect_cat = factor(effect_cat, 
                               levels = c("Strong Positive", 
                                          "Moderate Positive", 
                                          "Weak Positive", 
                                          "Weak Negative", 
                                          "Negative")))
  
  A.ind.cat <- A.ind %>%
    mutate(effect_cat = case_when(
      estimate >= 0.5 ~ "Strong Positive",
      estimate >= 0.25 ~ "Moderate Positive",
      estimate >= 0 ~ "Weak Positive",
      estimate >= -0.25 ~ "Weak Negative",
      TRUE ~ "Negative"
    )) %>%
    mutate(effect_cat = factor(effect_cat, 
                               levels = c("Strong Positive", 
                                          "Moderate Positive", 
                                          "Weak Positive", 
                                          "Weak Negative", 
                                          "Negative")))
  
  ind.plot.v3 <- ggplot() +
    geom_line(data = NA.ind.cat, alpha = 0.5, 
              aes(model, estimate, group = group_number, color = effect_cat)) +
    geom_line(data = A.ind.cat, alpha = 0.5, 
              aes(model, estimate, group = group_number, color = effect_cat)) +
    geom_errorbar(data = NA.ind.cat,
                  aes(model, estimate, group = model),
                  stat = "summary",
                  fun.data = function(x) {
                    mean_val <- mean(x)
                    sd_val <- sd(x)
                    data.frame(y = mean_val, ymin = mean_val - sd_val, ymax = mean_val + sd_val)
                  },
                  width = 0.1,
                  linewidth = 0.5,
                  color = "black") +
    geom_errorbar(data = A.ind.cat,
                  aes(model, estimate, group = model),
                  stat = "summary",
                  fun.data = function(x) {
                    mean_val <- mean(x)
                    sd_val <- sd(x)
                    data.frame(y = mean_val, ymin = mean_val - sd_val, ymax = mean_val + sd_val)
                  },
                  width = 0.1,
                  linewidth = 0.5,
                  color = "black") +
    geom_point(data = NA.ind.cat,
               aes(model, estimate, group = model),
               stat = "summary",
               fun = mean,
               size = 2,
               shape = 16,
               color = "black",
               fill = "black") +
    geom_point(data = A.ind.cat,
               aes(model, estimate, group = model),
               stat = "summary",
               fun = mean,
               size = 2,
               shape = 16,
               color = "black",
               fill = "black") +
    facet_grid(Scenario ~ domain, scales = "free_y") +
    theme_bw() +
    theme(
      text = element_text(family = "Times", size = 20, color = "black"),
      axis.text = element_text(family = "Times", size = 15),
      plot.title = element_text(family = "Times", face = "bold", size = 16),
      plot.subtitle = element_text(family = "Times", size = 15),
      axis.title.x = element_blank(),
      axis.title.y = element_text(family = "Times", size = 20, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      strip.text.x = element_text(size = 15, face = "bold"),
      legend.position = "right") +
    scale_color_manual(
      name = "Effect Size",
      values = c("Strong Positive" = "#660000",
                 "Moderate Positive" = "#cc0000",
                 "Weak Positive" = "#e69f00",
                 "Weak Negative" = "#cccccc",
                 "Negative" = "#9a9a9a"),
      breaks = c("Strong Positive", "Moderate Positive", "Weak Positive", 
                 "Weak Negative", "Negative")) +
    labs(y = "Estimate")
  
  # ind.plot.v3
  png("/home/hhampson/Results/Microbiome Results/Individual_Plot_Option3.png",res=300,height=4000,width=4000)
  plot(ind.plot.v3)
  dev.off()
}

#7. Sensitivity Measures ----
#GO back and make bahzing and ridge pvalues based on bcis not pvals
##A. Individual----
individual <- all_formatted_results %>% 
  filter(exposure!="Mixture")
# data <- individual
#Update function
sens.fxn.old <- function(data){
  BaH_ZING <- data %>% 
    filter(Model=="BaH-ZING") %>% 
    filter(Exposure!="Mixture") %>%
    # filter(grepl("species",Taxa)) %>% 
    mutate(PosNeg = case_when(indicator=="Causal" & P_Value=="*" ~ "True Positive",
                              indicator=="Causal" & P_Value!="*" ~ "False Negative",
                              indicator=="Non-Causal" & P_Value=="*" ~ "False Positive",
                              indicator=="Non-Causal" & P_Value!="*" ~ "True Negative")) %>% 
    group_by(Taxa.Level,Scenario) %>%
    # group_by(Taxa.Level) %>%
    summarise(
      Specificity = (length(which(PosNeg=="True Negative")) / (length(which(PosNeg=="True Negative")) + length(which(PosNeg=="False Positive")))),
      Sensitivity = (length(which(PosNeg=="True Positive")) / (length(which(PosNeg=="True Positive")) + length(which(PosNeg=="False Negative")))),
      PPV = (length(which(PosNeg=="True Positive")) / (length(which(PosNeg=="True Positive")) + length(which(PosNeg=="False Positive")))),
      NPV = (length(which(PosNeg=="True Negative")) / (length(which(PosNeg=="False Negative")) + length(which(PosNeg=="True Negative")))),
      FDR =  (length(which(PosNeg=="False Positive")) / (length(which(PosNeg=="True Positive")) + length(which(PosNeg=="False Positive")))),
      FDR2 = (1 - PPV)*100,
      TypeI = (1-Specificity),
      TypeII = 1-Sensitivity
    ) %>%  
    mutate(Model="BaH-ZING") %>% 
    ungroup()
  # select(-Specificity)
  
  measures <- c("Specificity","Sensitivity","PPV","NPV","FDR","TypeI","TypeII")
  
  
  ZING <- data %>% 
    # filter(Model=="ZINB" | Model=="Poisson") %>% 
    filter(Model=="ZING") %>% 
    filter(Exposure!="Mixture") %>%
    # filter(grepl("species",Taxa)) %>% 
    mutate(PosNeg = case_when(indicator=="Causal" & P_Value=="*" ~ "True Positive",
                              indicator=="Causal" & P_Value!="*" ~ "False Negative",
                              indicator=="Non-Causal" & P_Value=="*" ~ "False Positive",
                              indicator=="Non-Causal" & P_Value!="*" ~ "True Negative")) %>% 
    group_by(Taxa.Level,Scenario) %>%
    # group_by(Taxa.Level) %>%
    summarise(
      Specificity = (length(which(PosNeg=="True Negative")) / (length(which(PosNeg=="True Negative")) + length(which(PosNeg=="False Positive")))),
      Sensitivity = (length(which(PosNeg=="True Positive")) / (length(which(PosNeg=="True Positive")) + length(which(PosNeg=="False Negative")))),
      PPV = (length(which(PosNeg=="True Positive")) / (length(which(PosNeg=="True Positive")) + length(which(PosNeg=="False Positive")))),
      NPV = (length(which(PosNeg=="True Negative")) / (length(which(PosNeg=="False Negative")) + length(which(PosNeg=="True Negative")))),
      FDR =  (length(which(PosNeg=="False Positive")) / (length(which(PosNeg=="True Positive")) + length(which(PosNeg=="False Positive")))),
      FDR2 = (1 - PPV)*100,
      TypeI = (1-Specificity),
      TypeII = 1-Sensitivity
    ) %>% 
    mutate(Model="ZING") %>% 
    ungroup()
  # select(-Specificity)
  
  ZING.fdr <- data %>% 
    # filter(Model=="ZINB" | Model=="Poisson") %>% 
    filter(Model=="ZING") %>% 
    filter(Exposure!="Mixture") %>%
    # filter(grepl("species",Taxa)) %>% 
    mutate(PosNeg = case_when(indicator=="Causal" & fdr.p=="*" ~ "True Positive",
                              indicator=="Causal" & fdr.p!="*" ~ "False Negative",
                              indicator=="Non-Causal" & fdr.p=="*" ~ "False Positive",
                              indicator=="Non-Causal" & fdr.p!="*" ~ "True Negative")) %>% 
    group_by(Taxa.Level,Scenario) %>%
    # group_by(Taxa.Level) %>%
    summarise(
      Specificity = (length(which(PosNeg=="True Negative")) / (length(which(PosNeg=="True Negative")) + length(which(PosNeg=="False Positive")))),
      Sensitivity = (length(which(PosNeg=="True Positive")) / (length(which(PosNeg=="True Positive")) + length(which(PosNeg=="False Negative")))),
      PPV = (length(which(PosNeg=="True Positive")) / (length(which(PosNeg=="True Positive")) + length(which(PosNeg=="False Positive")))),
      NPV = (length(which(PosNeg=="True Negative")) / (length(which(PosNeg=="False Negative")) + length(which(PosNeg=="True Negative")))),
      FDR =  (length(which(PosNeg=="False Positive")) / (length(which(PosNeg=="True Positive")) + length(which(PosNeg=="False Positive")))),
      FDR2 = (1 - PPV)*100,
      TypeI = (1-Specificity),
      TypeII = 1-Sensitivity
    ) %>% 
    mutate(Model="Adj.ZING") %>% 
    ungroup() 
  # select(-Specificity)
  
  # SensSpec <- rbind(BaH_ZING,RBaH_ZING)
  SensSpec <- rbind(BaH_ZING,ZING)
  # SensSpec <- rbind(SensSpec,ZING)
  SensSpec <- rbind(SensSpec,ZING.fdr)
  # SensSpec <- rbind(SensSpec,ZING.fdr)
  # SensSpec.s <- SensSpec %>%
  #   filter(Taxa.Level=="Species")
  # SensSpec.g <- SensSpec %>%
  #   filter(Taxa.Level=="Genus")
  # SensSpec.f <- SensSpec %>%
  #   filter(Taxa.Level=="Family")
  # SensSpec.o <- SensSpec %>%
  #   filter(Taxa.Level=="Order")
  # SensSpec.c <- SensSpec %>%
  #   filter(Taxa.Level=="Class")
  # SensSpec.p <- SensSpec %>%
  #   filter(Taxa.Level=="Phylum")
  
  # B <- SensSpec %>%
  #   filter(Model=="BaH-ZING")
  # mean(B$Specificity,na.rm=T)#0.998
  # range(B$Specificity,na.rm=T)
  # mean(B$PPV,na.rm=T)#0.75
  # mean(B$Sensitivity,na.rm=T) #0.37
  # mean(B$TypeI,na.rm=T)*100 #0.19%
  # mean(B$TypeII,na.rm=T)*100 #63%
  # 
  # R <- SensSpec.s %>% 
  #   filter(Model=="Ridge ZING")
  # mean(R$Specificity,na.rm=T) #0.995
  # range(R$Specificity,na.rm=T)
  # mean(R$PPV,na.rm=T) #0.85
  # mean(R$Sensitivity,na.rm=T) #100
  # mean(R$TypeI,na.rm=T)*100 #0.48% 
  # mean(R$TypeII,na.rm=T)*100 #0.024% 
  # 
  # Z <- SensSpec %>% 
  #   filter(Model=="ZING")
  # mean(Z$Specificity,na.rm=T)
  # range(Z$Specificity,na.rm=T)
  # mean(Z$PPV,na.rm=T)
  # mean(Z$Sensitivity,na.rm=T)
  # mean(Z$TypeI,na.rm=T)*100 #25% 
  # mean(Z$TypeII,na.rm=T)*100 #0.006%
  # 
  # Z.a <- SensSpec %>% 
  #   filter(Model=="Adj.ZING")
  # mean(Z.a$Specificity,na.rm=T)
  # range(Z.a$Specificity,na.rm=T)
  # mean(Z.a$PPV,na.rm=T)
  # mean(Z.a$Sensitivity,na.rm=T)
  # mean(Z.a$TypeI,na.rm=T)*100 #20% 
  # mean(Z.a$TypeII,na.rm=T)*100 #0.226%
  
  # taxa_order <- c("Species","Genus","Family","Order","Class","Phylum")
  # SensSpec$Taxa.Level <- factor(SensSpec$Taxa.Level, levels = taxa_order)
  # write.csv(SensSpec,fs::path(dir.results,"SensitivitySpecificity_All.csv"))
  measures <- c("Specificity","Sensitivity","PPV","NPV","FDR","TypeI","TypeII")
  # measures <- c("PPV","FDR","TypeI")
  # SensSpec.Sum <-SensSpec %>% 
  #   tidylog::select(-c(Taxa.Level,Scenario)) %>% 
  #   tidylog::group_by(Model) %>% 
  #   tidylog::summarize(across(all_of(measures),~mean(.,na.rm=T))) 
  
  
  #   return(SensSpec)
  # }
  
  #Bias non causal and causal...and add MSE? Variance of mean?
  # Melt the data for ggplot2
  # SensSpec2 <- SensSpec %>%
  #   select(-FDR2) 
  # select(-Scenario)
  # select(FDR,Model,Scenario:Type2ErrorRate) %>% 
  # mutate(Scenario=paste0("Scenario ",Scenario))
  # SensitivitySpecificity2$FDR[1:4] <- NA
  # scenario.order <- c("Scenario 1", "Scenario 2", "Scenario 3","Scenario 4",
  #                     "Scenario 5","Scenario 6","Scenario 7","Scenario 8","Scenario 11")
  # Convert Scenario to a factor with the desired order
  # SensitivitySpecificity2$Scenario <- factor(SensitivitySpecificity2$Scenario, levels = scenario.order)
  
  # melted_data <- reshape2::melt(SensSpec2, id.vars = c('Taxa.Level', 'Model','Scenario'))
  melted_data <- reshape2::melt(SensSpec, id.vars = c('Model','Taxa.Level','Scenario'))
  # melted_data <- melted_data %>% 
  # filter(Scenario!="Scenario 11")
  # model_order <- c('Ridge ZING', 'BaH-ZING', 'ZING','Adj.ZING')
  model_order <- c('BaH-ZING', 'ZING','Adj.ZING')
  melted_data$Model <- factor(melted_data$Model, levels = model_order)
  taxa_order <- c("Species","Genus","Family","Order","Class","Phylum")
  # melted_data$Taxa.Level <- factor(melted_data$Taxa.Level, levels = taxa_order)
  
  
  #Order Taxa levels
  # taxa_order
  melted_data$Taxa.Level <- factor(melted_data$Taxa.Level, levels = taxa_order)
  return(melted_data)
}
sens.fxn <- function(data){
  BaH_ZING <- data %>% 
    filter(model=="BaHZING") %>% 
    filter(exposure!="Mixture") %>%
    mutate(PosNeg = case_when(indicator=="Associated" & sig=="*" ~ "True Positive",
                              indicator=="Associated" & sig!="*" ~ "False Negative",
                              indicator=="Not Associated" & sig=="*" ~ "False Positive",
                              indicator=="Not Associated" & sig!="*" ~ "True Negative")) %>% 
    group_by(domain, Scenario) %>%
    summarise(
      Specificity = (length(which(PosNeg=="True Negative")) / (length(which(PosNeg=="True Negative")) + length(which(PosNeg=="False Positive")))),
      Sensitivity = (length(which(PosNeg=="True Positive")) / (length(which(PosNeg=="True Positive")) + length(which(PosNeg=="False Negative")))),
      PPV = (length(which(PosNeg=="True Positive")) / (length(which(PosNeg=="True Positive")) + length(which(PosNeg=="False Positive")))),
      NPV = (length(which(PosNeg=="True Negative")) / (length(which(PosNeg=="False Negative")) + length(which(PosNeg=="True Negative")))),
      FDR =  (length(which(PosNeg=="False Positive")) / (length(which(PosNeg=="True Positive")) + length(which(PosNeg=="False Positive")))),
      FDR2 = (1 - PPV)*100,
      TypeI = (1-Specificity),
      TypeII = 1-Sensitivity,
      .groups = 'drop'
    ) %>%  
    mutate(Model="BaHZING")
  
  RBaH_ZING <- data %>% 
    filter(model=="RBaHZING") %>% 
    filter(exposure!="Mixture") %>%
    mutate(PosNeg = case_when(indicator=="Associated" & sig=="*" ~ "True Positive",
                              indicator=="Associated" & sig!="*" ~ "False Negative",
                              indicator=="Not Associated" & sig=="*" ~ "False Positive",
                              indicator=="Not Associated" & sig!="*" ~ "True Negative")) %>% 
    group_by(domain, Scenario) %>%
    summarise(
      Specificity = (length(which(PosNeg=="True Negative")) / (length(which(PosNeg=="True Negative")) + length(which(PosNeg=="False Positive")))),
      Sensitivity = (length(which(PosNeg=="True Positive")) / (length(which(PosNeg=="True Positive")) + length(which(PosNeg=="False Negative")))),
      PPV = (length(which(PosNeg=="True Positive")) / (length(which(PosNeg=="True Positive")) + length(which(PosNeg=="False Positive")))),
      NPV = (length(which(PosNeg=="True Negative")) / (length(which(PosNeg=="False Negative")) + length(which(PosNeg=="True Negative")))),
      FDR =  (length(which(PosNeg=="False Positive")) / (length(which(PosNeg=="True Positive")) + length(which(PosNeg=="False Positive")))),
      FDR2 = (1 - PPV)*100,
      TypeI = (1-Specificity),
      TypeII = 1-Sensitivity,
      .groups = 'drop'
    ) %>% 
    mutate(Model="RBaHZING")
  
  ZING <- data %>% 
    filter(model=="ZING") %>% 
    filter(exposure!="Mixture") %>%
    mutate(PosNeg = case_when(indicator=="Associated" & sig=="*" ~ "True Positive",
                              indicator=="Associated" & sig!="*" ~ "False Negative",
                              indicator=="Not Associated" & sig=="*" ~ "False Positive",
                              indicator=="Not Associated" & sig!="*" ~ "True Negative")) %>% 
    group_by(domain, Scenario) %>%
    summarise(
      Specificity = (length(which(PosNeg=="True Negative")) / (length(which(PosNeg=="True Negative")) + length(which(PosNeg=="False Positive")))),
      Sensitivity = (length(which(PosNeg=="True Positive")) / (length(which(PosNeg=="True Positive")) + length(which(PosNeg=="False Negative")))),
      PPV = (length(which(PosNeg=="True Positive")) / (length(which(PosNeg=="True Positive")) + length(which(PosNeg=="False Positive")))),
      NPV = (length(which(PosNeg=="True Negative")) / (length(which(PosNeg=="False Negative")) + length(which(PosNeg=="True Negative")))),
      FDR =  (length(which(PosNeg=="False Positive")) / (length(which(PosNeg=="True Positive")) + length(which(PosNeg=="False Positive")))),
      FDR2 = (1 - PPV)*100,
      TypeI = (1-Specificity),
      TypeII = 1-Sensitivity,
      .groups = 'drop'
    ) %>% 
    mutate(Model="ZING")
  
  ZING.fdr <- data %>% 
    filter(model=="ZING") %>% 
    filter(exposure!="Mixture") %>%
    mutate(PosNeg = case_when(indicator=="Associated" & fdr.sig=="*" ~ "True Positive",
                              indicator=="Associated" & fdr.sig!="*" ~ "False Negative",
                              indicator=="Not Associated" & fdr.sig=="*" ~ "False Positive",
                              indicator=="Not Associated" & fdr.sig!="*" ~ "True Negative")) %>% 
    group_by(domain, Scenario) %>%
    summarise(
      Specificity = (length(which(PosNeg=="True Negative")) / (length(which(PosNeg=="True Negative")) + length(which(PosNeg=="False Positive")))),
      Sensitivity = (length(which(PosNeg=="True Positive")) / (length(which(PosNeg=="True Positive")) + length(which(PosNeg=="False Negative")))),
      PPV = (length(which(PosNeg=="True Positive")) / (length(which(PosNeg=="True Positive")) + length(which(PosNeg=="False Positive")))),
      NPV = (length(which(PosNeg=="True Negative")) / (length(which(PosNeg=="False Negative")) + length(which(PosNeg=="True Negative")))),
      FDR =  (length(which(PosNeg=="False Positive")) / (length(which(PosNeg=="True Positive")) + length(which(PosNeg=="False Positive")))),
      FDR2 = (1 - PPV)*100,
      TypeI = (1-Specificity),
      TypeII = 1-Sensitivity,
      .groups = 'drop'
    ) %>% 
    mutate(Model="Adj.ZING")
  
  SensSpec <- rbind(BaH_ZING, RBaH_ZING, ZING, ZING.fdr)
  
  measures <- c("Specificity","Sensitivity","PPV","NPV","FDR","TypeI","TypeII")
  
  melted_data <- reshape2::melt(SensSpec, id.vars = c('Model','domain','Scenario'))
  
  model_order <- c('BaHZING', 'RBaHZING', 'ZING','Adj.ZING')
  melted_data$Model <- factor(melted_data$Model, levels = model_order)
  
  taxa_order <- c("Species","Genus","Family","Order","Class","Phylum")
  melted_data$domain <- factor(melted_data$domain, levels = taxa_order)
  
  return(melted_data)
}
formatted.all <- sens.fxn(individual)
rm(individual)
formatted.all <- formatted.all %>% 
  mutate(Scenario=paste0("Scenario ",Scenario))
formatted.all$Scenario <- factor(formatted.all$Scenario,levels=scenario_order)

formatted.wide <- formatted.all %>%
  pivot_wider(names_from = domain,
              values_from=value) 
# select(-Scenario)
formatted.wide <- formatted.wide[c("Model","Scenario","variable","Species","Genus",
                                   "Family","Order","Class","Phylum")]
# write.csv(formatted.wide,fs::path(dir.results,"Sensitivity_Individual.csv"))
# write.csv(formatted.wide,"/home/hhampson/Results/Microbiome Results/Sensitivity_Measures_Individual.csv")

# Assuming your data is already in long format with columns:
# Model, Scenario, variable, Taxonomic_Level, Specificity

# Make sure Taxonomic_Level is a factor in the right order
specificity_long <- formatted.all %>%
  mutate(domain = factor(domain, 
                         levels = c("Species", "Genus", "Family", 
                                    "Order", "Class", "Phylum")))

# PLOT 1: BOXPLOT - All metrics by Model and Taxonomic Level
boxplot_all_metrics <- ggplot(specificity_long, 
                              aes(x = domain, y = value, fill = Model)) +
  geom_boxplot(alpha = 0.7, outlier.shape = 16, outlier.size = 1) +
  facet_grid(variable ~ ., scales = "free_y") +
  theme_bw() +
  theme(
    text = element_text(family = "Times", size = 14),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 14, face = "bold"),
    strip.text.y = element_text(size = 12, face = "bold"),
    legend.position = "right",
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11)) +
  scale_fill_manual(values = c("BaHZING" = "#2166ac", 
                               "RBaHZING" = "#92c5de", 
                               "ZING" = "#b2182b",
                               "Adj.ZING" = "gray")) +
  labs(x = "Taxonomic Level", 
       y = "Value",
       fill = "Model")

png("/home/hhampson/Results/Microbiome Results/Individual_Metrics_BoxPlots.png",res=300,width=2000,height=1000)
plot(boxplot_all_metrics)
dev.off()

##B. Mixture----
mixture.sens <- all_formatted_results %>% 
  filter(exposure=="Mixture")
data <- mixture.sens
sens.fxn.mix.old <- function(data){
  BaH_ZING <- data %>% 
    filter(Model=="BaH-ZING") %>% 
    filter(Exposure=="Mixture") %>%
    # filter(grepl("species",Taxa)) %>% 
    mutate(PosNeg = case_when(indicator=="Causal" & P_Value=="*" ~ "True Positive",
                              indicator=="Causal" & P_Value!="*" ~ "False Negative",
                              indicator=="Non-Causal" & P_Value=="*" ~ "False Positive",
                              indicator=="Non-Causal" & P_Value!="*" ~ "True Negative")) %>% 
    group_by(Taxa.Level,Scenario) %>%
    # group_by(Taxa.Level) %>%
    summarise(
      Specificity = (length(which(PosNeg=="True Negative")) / (length(which(PosNeg=="True Negative")) + length(which(PosNeg=="False Positive")))),
      Sensitivity = (length(which(PosNeg=="True Positive")) / (length(which(PosNeg=="True Positive")) + length(which(PosNeg=="False Negative")))),
      PPV = (length(which(PosNeg=="True Positive")) / (length(which(PosNeg=="True Positive")) + length(which(PosNeg=="False Positive")))),
      # NPV = (length(which(PosNeg=="True Negative")) / (length(which(PosNeg=="False Negative")) + length(which(PosNeg=="True Negative")))),
      FDR =  (length(which(PosNeg=="False Positive")) / (length(which(PosNeg=="True Positive")) + length(which(PosNeg=="False Positive")))),
      # FDR2 = (1 - PPV)*100,
      TypeI = (1-Specificity),
      TypeII = 1-Sensitivity
    ) %>%  
    mutate(Model="BaH-ZING") %>% 
    ungroup()
  
  measures <- c("Specificity","Sensitivity","PPV","NPV","FDR","TypeI","TypeII")
  # measures <- c("PPV","FDR")
  
  # RBaH_ZING <- data %>% 
  #   filter(Model=="Ridge ZING") %>%
  #   filter(Exposure!="Mixture") %>%
  #   # filter(grepl("species",Taxa)) %>% 
  #   mutate(PosNeg = case_when(indicator=="Causal" & P_Value=="*" ~ "True Positive",
  #                             indicator=="Causal" & P_Value!="*" ~ "False Negative",
  #                             indicator=="Non-Causal" & P_Value=="*" ~ "False Positive",
  #                             indicator=="Non-Causal" & P_Value!="*" ~ "True Negative")) %>% 
  #   group_by(Taxa.Level,Scenario) %>%
  #   # group_by(Taxa.Level) %>%
  #   summarise(
  #     Specificity = (length(which(PosNeg=="True Negative")) / (length(which(PosNeg=="True Negative")) + length(which(PosNeg=="False Positive")))),
  #     Sensitivity = (length(which(PosNeg=="True Positive")) / (length(which(PosNeg=="True Positive")) + length(which(PosNeg=="False Negative")))),
  #     PPV = (length(which(PosNeg=="True Positive")) / (length(which(PosNeg=="True Positive")) + length(which(PosNeg=="False Positive")))),
  #     # NPV = (length(which(PosNeg=="True Negative")) / (length(which(PosNeg=="False Negative")) + length(which(PosNeg=="True Negative")))),
  #     FDR =  (length(which(PosNeg=="False Positive")) / (length(which(PosNeg=="True Positive")) + length(which(PosNeg=="False Positive")))),
  #     # FDR2 = (1 - PPV)*100,
  #     TypeI = (1-Specificity),
  #     TypeII = 1-Sensitivity
  #   ) %>% 
  #   mutate(across(all_of(measures),~ifelse(Taxa.Level!="Species",NA,.))) %>% 
  #   mutate(Model="Ridge ZING") %>% 
  #   ungroup()
  
  ZING <- data %>% 
    # filter(Model=="ZINB" | Model=="Poisson") %>% 
    filter(Model=="ZING") %>% 
    filter(Exposure=="Mixture") %>%
    # filter(grepl("species",Taxa)) %>% 
    mutate(PosNeg = case_when(indicator=="Causal" & P_Value=="*" ~ "True Positive",
                              indicator=="Causal" & P_Value!="*" ~ "False Negative",
                              indicator=="Non-Causal" & P_Value=="*" ~ "False Positive",
                              indicator=="Non-Causal" & P_Value!="*" ~ "True Negative")) %>% 
    group_by(Taxa.Level,Scenario) %>%
    # group_by(Taxa.Level) %>%
    summarise(
      Specificity = (length(which(PosNeg=="True Negative")) / (length(which(PosNeg=="True Negative")) + length(which(PosNeg=="False Positive")))),
      Sensitivity = (length(which(PosNeg=="True Positive")) / (length(which(PosNeg=="True Positive")) + length(which(PosNeg=="False Negative")))),
      PPV = (length(which(PosNeg=="True Positive")) / (length(which(PosNeg=="True Positive")) + length(which(PosNeg=="False Positive")))),
      # NPV = (length(which(PosNeg=="True Negative")) / (length(which(PosNeg=="False Negative")) + length(which(PosNeg=="True Negative")))),
      FDR =  (length(which(PosNeg=="False Positive")) / (length(which(PosNeg=="True Positive")) + length(which(PosNeg=="False Positive")))),
      # FDR2 = (1 - PPV)*100,
      TypeI = (1-Specificity),
      TypeII = 1-Sensitivity
    ) %>% 
    mutate(Model="ZING") %>% 
    ungroup()
  
  ZING.fdr <- data %>% 
    # filter(Model=="ZINB" | Model=="Poisson") %>% 
    filter(Model=="ZING") %>% 
    filter(Exposure=="Mixture") %>%
    # filter(grepl("species",Taxa)) %>% 
    mutate(PosNeg = case_when(indicator=="Causal" & fdr.p=="*" ~ "True Positive",
                              indicator=="Causal" & fdr.p!="*" ~ "False Negative",
                              indicator=="Non-Causal" & fdr.p=="*" ~ "False Positive",
                              indicator=="Non-Causal" & fdr.p!="*" ~ "True Negative")) %>% 
    group_by(Taxa.Level,Scenario) %>%
    # group_by(Taxa.Level) %>%
    summarise(
      Specificity = (length(which(PosNeg=="True Negative")) / (length(which(PosNeg=="True Negative")) + length(which(PosNeg=="False Positive")))),
      Sensitivity = (length(which(PosNeg=="True Positive")) / (length(which(PosNeg=="True Positive")) + length(which(PosNeg=="False Negative")))),
      PPV = (length(which(PosNeg=="True Positive")) / (length(which(PosNeg=="True Positive")) + length(which(PosNeg=="False Positive")))),
      # NPV = (length(which(PosNeg=="True Negative")) / (length(which(PosNeg=="False Negative")) + length(which(PosNeg=="True Negative")))),
      FDR =  (length(which(PosNeg=="False Positive")) / (length(which(PosNeg=="True Positive")) + length(which(PosNeg=="False Positive")))),
      # FDR2 = (1 - PPV)*100,
      TypeI = (1-Specificity),
      TypeII = 1-Sensitivity
    ) %>% 
    mutate(Model="Adj.ZING") %>% 
    ungroup()
  
  # SensSpec <- rbind(BaH_ZING,RBaH_ZING)
  SensSpec <- rbind(BaH_ZING,ZING)
  # SensSpec <- rbind(SensSpec,ZING)
  SensSpec <- rbind(SensSpec,ZING.fdr)
  
  # melted_data <- reshape2::melt(SensSpec2, id.vars = c('Taxa.Level', 'Model','Scenario'))
  melted_data <- reshape2::melt(SensSpec, id.vars = c('Model','Taxa.Level','Scenario'))
  # melted_data <- melted_data %>% 
  # filter(Scenario!="Scenario 11")
  model_order <- c( 'BaH-ZING', 'ZING','Adj.ZING')
  # model_order <- c('BaH-ZING', 'ZING','Adj.ZING')
  melted_data$Model <- factor(melted_data$Model, levels = model_order)
  taxa_order <- c("Species","Genus","Family","Order","Class","Phylum")
  # melted_data$Taxa.Level <- factor(melted_data$Taxa.Level, levels = taxa_order)
  
  melted_data$Taxa.Level <- factor(melted_data$Taxa.Level, levels = taxa_order)
  return(melted_data)
}
sens.fxn.mixture <- function(data){
  BaH_ZING <- data %>% 
    filter(model=="BaHZING") %>% 
    filter(exposure=="Mixture") %>%
    mutate(PosNeg = case_when(indicator=="Associated" & sig=="*" ~ "True Positive",
                              indicator=="Associated" & sig!="*" ~ "False Negative",
                              indicator=="Not Associated" & sig=="*" ~ "False Positive",
                              indicator=="Not Associated" & sig!="*" ~ "True Negative")) %>% 
    group_by(domain, Scenario) %>%
    summarise(
      Specificity = (length(which(PosNeg=="True Negative")) / (length(which(PosNeg=="True Negative")) + length(which(PosNeg=="False Positive")))),
      Sensitivity = (length(which(PosNeg=="True Positive")) / (length(which(PosNeg=="True Positive")) + length(which(PosNeg=="False Negative")))),
      PPV = (length(which(PosNeg=="True Positive")) / (length(which(PosNeg=="True Positive")) + length(which(PosNeg=="False Positive")))),
      NPV = (length(which(PosNeg=="True Negative")) / (length(which(PosNeg=="False Negative")) + length(which(PosNeg=="True Negative")))),
      FDR =  (length(which(PosNeg=="False Positive")) / (length(which(PosNeg=="True Positive")) + length(which(PosNeg=="False Positive")))),
      FDR2 = (1 - PPV)*100,
      TypeI = (1-Specificity),
      TypeII = 1-Sensitivity,
      .groups = 'drop'
    ) %>%  
    mutate(Model="BaHZING")
  
  RBaH_ZING <- data %>% 
    filter(model=="RBaHZING") %>% 
    filter(exposure=="Mixture") %>%
    mutate(PosNeg = case_when(indicator=="Associated" & sig=="*" ~ "True Positive",
                              indicator=="Associated" & sig!="*" ~ "False Negative",
                              indicator=="Not Associated" & sig=="*" ~ "False Positive",
                              indicator=="Not Associated" & sig!="*" ~ "True Negative")) %>% 
    group_by(domain, Scenario) %>%
    summarise(
      Specificity = (length(which(PosNeg=="True Negative")) / (length(which(PosNeg=="True Negative")) + length(which(PosNeg=="False Positive")))),
      Sensitivity = (length(which(PosNeg=="True Positive")) / (length(which(PosNeg=="True Positive")) + length(which(PosNeg=="False Negative")))),
      PPV = (length(which(PosNeg=="True Positive")) / (length(which(PosNeg=="True Positive")) + length(which(PosNeg=="False Positive")))),
      NPV = (length(which(PosNeg=="True Negative")) / (length(which(PosNeg=="False Negative")) + length(which(PosNeg=="True Negative")))),
      FDR =  (length(which(PosNeg=="False Positive")) / (length(which(PosNeg=="True Positive")) + length(which(PosNeg=="False Positive")))),
      FDR2 = (1 - PPV)*100,
      TypeI = (1-Specificity),
      TypeII = 1-Sensitivity,
      .groups = 'drop'
    ) %>% 
    mutate(Model="RBaHZING")
  
  ZING <- data %>% 
    filter(model=="ZING") %>% 
    filter(exposure=="Mixture") %>%
    mutate(PosNeg = case_when(indicator=="Associated" & sig=="*" ~ "True Positive",
                              indicator=="Associated" & sig!="*" ~ "False Negative",
                              indicator=="Not Associated" & sig=="*" ~ "False Positive",
                              indicator=="Not Associated" & sig!="*" ~ "True Negative")) %>% 
    group_by(domain, Scenario) %>%
    summarise(
      Specificity = (length(which(PosNeg=="True Negative")) / (length(which(PosNeg=="True Negative")) + length(which(PosNeg=="False Positive")))),
      Sensitivity = (length(which(PosNeg=="True Positive")) / (length(which(PosNeg=="True Positive")) + length(which(PosNeg=="False Negative")))),
      PPV = (length(which(PosNeg=="True Positive")) / (length(which(PosNeg=="True Positive")) + length(which(PosNeg=="False Positive")))),
      NPV = (length(which(PosNeg=="True Negative")) / (length(which(PosNeg=="False Negative")) + length(which(PosNeg=="True Negative")))),
      FDR =  (length(which(PosNeg=="False Positive")) / (length(which(PosNeg=="True Positive")) + length(which(PosNeg=="False Positive")))),
      FDR2 = (1 - PPV)*100,
      TypeI = (1-Specificity),
      TypeII = 1-Sensitivity,
      .groups = 'drop'
    ) %>% 
    mutate(Model="ZING")
  
  # Note: No FDR-adjusted version for ZING mixture since there's only one test per domain
  
  SensSpec <- rbind(BaH_ZING, RBaH_ZING, ZING)
  
  measures <- c("Specificity","Sensitivity","PPV","NPV","FDR","TypeI","TypeII")
  
  melted_data <- reshape2::melt(SensSpec, id.vars = c('Model','domain','Scenario'))
  
  model_order <- c('BaHZING', 'RBaHZING', 'ZING')
  melted_data$Model <- factor(melted_data$Model, levels = model_order)
  
  taxa_order <- c("Species","Genus","Family","Order","Class","Phylum")
  melted_data$domain <- factor(melted_data$domain, levels = taxa_order)
  
  return(melted_data)
}
formatted.all.mix <- sens.fxn.mixture(mixture.sens)

formatted.all.mix <- formatted.all.mix %>% 
  mutate(Scenario=paste0("Scenario ",Scenario))
formatted.all.mix$Scenario <- factor(formatted.all.mix$Scenario,levels=scenario_order)

formatted.all.mix.wide <- formatted.all.mix %>%  
  pivot_wider(names_from = Taxa.Level,
              values_from=value) %>% 
  select(-Scenario)
formatted.all.mix.wide<- formatted.all.mix.wide[c("Model","variable","Species","Genus",
                                                  "Family","Order","Class","Phylum")]
write.csv(formatted.wide,fs::path(dir.results,"Sensitivity_Mixture_Scenario15.csv"))


#8. Dendrograms ----
##A. Individual ----
#Scenarios
means <- formatted.data %>%
  tidylog::filter(Exposure!="Mixture") %>%
  tidylog::filter(Model!="Ridge ZING")  
# tidylog::filter(Model!="ZING")
means2 <- means

source("Microbiome_Cleaning.R")

key <- species.key %>%
  select(species,sim.num) %>%
  rename(Taxa=sim.num) %>%
  rename(taxa.name=species)
genus.key <- genus.key %>%
  rename(Taxa=sim.num,
         taxa.name=genus)
key <- rbind(key,genus.key)
family.key <- family.key %>%
  rename(Taxa=sim.num,
         taxa.name=family)
key <- rbind(key,family.key)
order.key <- order.key %>%
  rename(Taxa=sim.num,
         taxa.name=order)
key <- rbind(key,order.key)
class.key <- class.key %>%
  rename(Taxa=sim.num,
         taxa.name=class)
key <- rbind(key,class.key)
phylum.key <- phylum.key %>%
  rename(Taxa=sim.num,
         taxa.name=phylum)
key <- rbind(key,phylum.key)


means2 <- tidylog::left_join(means2,key,by="Taxa")
#Format for tree
# species.key$species.new
means2$taxa.name.new <- str_replace(means2$taxa.name,"_p_","--p_")
means2$taxa.name.new <- str_replace(means2$taxa.name.new,"_c_","--c_")
means2$taxa.name.new <- str_replace(means2$taxa.name.new,"_o_","--o_")
means2$taxa.name.new <- str_replace(means2$taxa.name.new,"_f_","--f_")
means2$taxa.name.new <- str_replace(means2$taxa.name.new,"_g_","--g_")
means2$taxa.name.new <- str_replace(means2$taxa.name.new,"_s_","--s_")

means2$Scenario <- paste0("Scenario ",means2$Scenario)
scenario14.dataA <- means2 %>%
  filter(Scenario=="Scenario 14") %>% 
  filter(Exposure=="X.1" |Exposure=="X.2" | Exposure =="X.3") %>% 
  filter(Model=="BaH-ZING") %>% 
  # filter(Exposure=="X.1") %>% 
  group_by(taxa.name.new) %>% 
  tidylog::summarize(Mean.total = mean(Mean)) %>% 
  filter(grepl("p__Bacteroidota",taxa.name.new)) %>% 
  mutate(mean.color = case_when(Mean.total<=0 ~ "<= 0",
                                Mean.total>0 & Mean.total <=0.2 ~ "0 < Mean <= 0.2",
                                Mean.total>0.2 & Mean.total <0.4 ~ "0.2 < Mean < 0.4",
                                Mean.total>=0.4 ~ "Mean >= 0.4")) %>% 
  mutate(mean.size= case_when(Mean.total<=0 ~ "<= 0",
                              Mean.total>0 & Mean.total <=0.2 ~ "0 < Mean <= 0.2",
                              Mean.total>0.2 & Mean.total <0.4 ~ "0.2 < Mean < 0.4",
                              Mean.total>=0.4 ~ "Mean >= 0.4")) %>% 
  mutate(mean.color = factor(mean.color,levels=c("<= 0","0 < Mean <= 0.2","0.2 < Mean < 0.4","Mean >= 0.4"))) %>% 
  mutate(mean.size = factor(mean.size,levels=c("<= 0","0 < Mean <= 0.2","0.2 < Mean < 0.4","Mean >= 0.4"))) %>% 
  mutate(taxa.name.new = ifelse(taxa.name.new=="d__Bacteria--p__Bacteroidota","Bacteria",taxa.name.new))

scenario15.dataA <- means2 %>%
  filter(Scenario=="Scenario 15") %>% 
  filter(Model=="BaH-ZING") %>% 
  # filter(Exposure=="X.1") %>% 
  group_by(taxa.name.new) %>% 
  tidylog::summarize(Mean.total = mean(Mean)) %>% 
  filter(grepl("p__Bacteroidota",taxa.name.new)) %>% 
  mutate(mean.color = case_when(Mean.total<=0 ~ "<= 0",
                                Mean.total>0 & Mean.total <=0.2 ~ "0 < Mean <= 0.2",
                                Mean.total>0.2 & Mean.total <0.4 ~ "0.2 < Mean < 0.4",
                                Mean.total>=0.4 ~ "Mean >= 0.4")) %>% 
  mutate(mean.size= case_when(Mean.total<=0 ~ "<= 0",
                              Mean.total>0 & Mean.total <=0.2 ~ "0 < Mean <= 0.2",
                              Mean.total>0.2 & Mean.total <0.4 ~ "0.2 < Mean < 0.4",
                              Mean.total>=0.4 ~ "Mean >= 0.4")) %>% 
  mutate(mean.color = factor(mean.color,levels=c("<= 0","0 < Mean <= 0.2","0.2 < Mean < 0.4","Mean >= 0.4"))) %>% 
  mutate(mean.size = factor(mean.size,levels=c("<= 0","0 < Mean <= 0.2","0.2 < Mean < 0.4","Mean >= 0.4"))) %>% 
  # mutate(taxa.name.new = str_replace(taxa.name.new,"d__Bacteria--p__Bacteroidota","Bacteria")) 
  mutate(taxa.name.new = ifelse(taxa.name.new=="d__Bacteria--p__Bacteroidota","Bacteria",taxa.name.new))

scen14.b <- ggtree(tree2,layout="fan",open.angle=230) %<+% scenario14.dataA +
  geom_point(aes(color=mean.color, size = mean.size), alpha=1)+
  scale_size_manual(values = c("<= 0" = 1, "0 < Mean <= 0.2" = 2, "0.2 < Mean < 0.4" = 3, "Mean >= 0.4" = 4))+
  scale_color_manual(values = c("<= 0"="#9a9a9a", "0 < Mean <= 0.2"="#e69f00", "0.2 < Mean < 0.4"="#cc0000", "Mean >= 0.4"="#660000"))+
  # scale_size_manual(values = c("small" = 1, "small.medium" = 2, "medium" = 3, "large" = 4))+
  # scale_color_manual(values = c("#9a9a9a"="#9a9a9a", "#e69f00"="#e69f00", "#cc0000"="#cc0000", "#660000"="#660000"))
  # scale_color_gradient2(low = "gray", mid = c("gray","yellow","orange","#CC0066"), 
  #                       high = "#CC0066", midpoint = 0.2, limits = c(min, max), na.value = NA) +
  theme(legend.position = "none")


scen15.b <- ggtree(tree2,layout="fan",open.angle=230) %<+% scenario15.dataA +
  geom_point(aes(color=mean.color, size = mean.size), alpha=1)+
  scale_size_manual(values = c("<= 0" = 1, "0 < Mean <= 0.2" = 2, "0.2 < Mean < 0.4" = 3, "Mean >= 0.4" = 4))+
  scale_color_manual(values = c("<= 0"="#9a9a9a", "0 < Mean <= 0.2"="#e69f00", "0.2 < Mean < 0.4"="#cc0000", "Mean >= 0.4"="#660000"))+
  
  # scen15 <- ggtree(tree2,layout="fan",open.angle=230) %<+% scenario15.dataA +
  #   geom_point(aes(color=Mean.total, size = Mean.total), alpha=1)+
  #   scale_colour_gradientn(colours = c("#660000","#cc0000","#cc0000","#e69f00","#9a9a9a"),
  #                          values = c(1.0,0.8,0.6,0.4,0.2,0))
  # scale_color_gradient2(low = "gray", mid = c("gray","yellow","orange","#CC0066"),
  #                       high = "#CC0066", midpoint = 0.2, limits = c(min, max), na.value = NA) +
  theme(legend.position = "none")

scen15.b
ggsave(plot=scen15.b,fs::path(dir.results,"Dendrogram_Individual_Scenario15.jpeg"))

##B. Mixture----
mixture2 <- formatted.data %>% 
  filter(Exposure=="Mixture") %>% 
  filter(Model=="BaH-ZING")
# mixture2 <- mixture 
# filter(Replication==1)
means2 <- tidylog::left_join(mixture2,key,by="Taxa")
#Format for tree
# species.key$species.new
means2$taxa.name.new <- str_replace(means2$taxa.name,"_p_","--p_")
means2$taxa.name.new <- str_replace(means2$taxa.name.new,"_c_","--c_")
means2$taxa.name.new <- str_replace(means2$taxa.name.new,"_o_","--o_")
means2$taxa.name.new <- str_replace(means2$taxa.name.new,"_f_","--f_")
means2$taxa.name.new <- str_replace(means2$taxa.name.new,"_g_","--g_")
means2$taxa.name.new <- str_replace(means2$taxa.name.new,"_s_","--s_")

#Scenario 1
# scenario1.data <- scenario.fxn(scenario1)
#Scenario 2
# scenario2.data <- scenario.fxn(scenario2)
means2$Scenario <- paste0("Scenario ",means2$Scenario)
scenario12.dataA <- means2 %>%
  filter(Scenario=="Scenario 12") %>% 
  filter(Model=="BaH-ZING") %>% 
  # filter(Exposure=="X.1") %>% 
  group_by(taxa.name.new) %>%
  tidylog::summarize(Mean.total = mean(Mean)) %>% 
  mutate(mean.color = case_when(Mean.total<=0 ~ "<= 0",
                                Mean.total>0 & Mean.total <=1 ~ "0 < Mean <= 1.0",
                                Mean.total>1 & Mean.total <2 ~ "1.0 < Mean < 2.0",
                                Mean.total>=2 ~ "Mean >= 2.0")) %>% 
  mutate(mean.size= case_when(Mean.total<=0 ~ "<= 0",
                              Mean.total>0 & Mean.total <=1 ~ "0 < Mean <= 1.0",
                              Mean.total>1 & Mean.total <2 ~ "1.0 < Mean < 2.0",
                              Mean.total>=2 ~ "Mean >= 2.0")) %>% 
  mutate(mean.color = factor(mean.color,levels=c("<= 0","0 < Mean <= 1.0","1.0 < Mean < 2.0","Mean >= 2.0"))) %>% 
  mutate(mean.size = factor(mean.size,levels=c("<= 0","0 < Mean <= 1.0","1.0 < Mean < 2.0","Mean >= 2.0"))) %>% 
  mutate(taxa.name.new = ifelse(taxa.name.new=="d__Bacteria--p__Bacteroidota","Bacteria",taxa.name.new))

# select(taxa.name.new,Mean)
#Scenario 13
scenario13.dataA <- means2 %>%
  filter(Scenario=="Scenario 13") %>% 
  filter(Model=="BaH-ZING") %>% 
  # filter(Exposure=="X.1") %>% 
  group_by(taxa.name.new) %>%
  tidylog::summarize(Mean.total = mean(Mean)) %>% 
  mutate(mean.color = case_when(Mean.total<=0 ~ "<= 0",
                                Mean.total>0 & Mean.total <=1 ~ "0 < Mean <= 1.0",
                                Mean.total>1 & Mean.total <2 ~ "1.0 < Mean < 2.0",
                                Mean.total>=2 ~ "Mean >= 2.0")) %>% 
  mutate(mean.size= case_when(Mean.total<=0 ~ "<= 0",
                              Mean.total>0 & Mean.total <=1 ~ "0 < Mean <= 1.0",
                              Mean.total>1 & Mean.total <2 ~ "1.0 < Mean < 2.0",
                              Mean.total>=2 ~ "Mean >= 2.0")) %>% 
  mutate(mean.color = factor(mean.color,levels=c("<= 0","0 < Mean <= 1.0","1.0 < Mean < 2.0","Mean >= 2.0"))) %>% 
  mutate(mean.size = factor(mean.size,levels=c("<= 0","0 < Mean <= 1.0","1.0 < Mean < 2.0","Mean >= 2.0"))) %>% 
  mutate(taxa.name.new = ifelse(taxa.name.new=="d__Bacteria--p__Bacteroidota","Bacteria",taxa.name.new))

# select(taxa.name.new,Mean)
#Scenario 14
scenario14.dataA <- means2 %>%
  filter(Scenario=="Scenario 14") %>% 
  filter(Model=="BaH-ZING") %>% 
  # filter(Exposure=="X.1") %>% 
  group_by(taxa.name.new) %>%
  tidylog::summarize(Mean.total = mean(Mean)) %>% 
  mutate(mean.color = case_when(Mean.total<=0 ~ "<= 0",
                                Mean.total>0 & Mean.total <=1 ~ "0 < Mean <= 1.0",
                                Mean.total>1 & Mean.total <2 ~ "1.0 < Mean < 2.0",
                                Mean.total>=2 ~ "Mean >= 2.0")) %>% 
  mutate(mean.size= case_when(Mean.total<=0 ~ "<= 0",
                              Mean.total>0 & Mean.total <=1 ~ "0 < Mean <= 1.0",
                              Mean.total>1 & Mean.total <2 ~ "1.0 < Mean < 2.0",
                              Mean.total>=2 ~ "Mean >= 2.0")) %>% 
  mutate(mean.color = factor(mean.color,levels=c("<= 0","0 < Mean <= 1.0","1.0 < Mean < 2.0","Mean >= 2.0"))) %>% 
  mutate(mean.size = factor(mean.size,levels=c("<= 0","0 < Mean <= 1.0","1.0 < Mean < 2.0","Mean >= 2.0"))) %>% 
  mutate(taxa.name.new = ifelse(taxa.name.new=="d__Bacteria--p__Bacteroidota","Bacteria",taxa.name.new))

# select(taxa.name.new,Mean)
#Scenario 15
scenario15.dataA <- means2 %>%
  filter(Scenario=="Scenario 15") %>% 
  filter(Model=="BaH-ZING") %>% 
  group_by(taxa.name.new) %>%
  tidylog::summarize(Mean.total = mean(Mean)) %>% 
  mutate(mean.color = case_when(Mean.total<=0 ~ "<= 0",
                                Mean.total>0 & Mean.total <=1 ~ "0 < Mean <= 1.0",
                                Mean.total>1 & Mean.total <2 ~ "1.0 < Mean < 2.0",
                                Mean.total>=2 ~ "Mean >= 2.0")) %>% 
  mutate(mean.size= case_when(Mean.total<=0 ~ "<= 0",
                              Mean.total>0 & Mean.total <=1 ~ "0 < Mean <= 1.0",
                              Mean.total>1 & Mean.total <2 ~ "1.0 < Mean < 2.0",
                              Mean.total>=2 ~ "Mean >= 2.0")) %>% 
  mutate(mean.color = factor(mean.color,levels=c("<= 0","0 < Mean <= 1.0","1.0 < Mean < 2.0","Mean >= 2.0"))) %>% 
  mutate(mean.size = factor(mean.size,levels=c("<= 0","0 < Mean <= 1.0","1.0 < Mean < 2.0","Mean >= 2.0"))) %>% 
  mutate(taxa.name.new = ifelse(taxa.name.new=="d__Bacteria--p__Bacteroidota","Bacteria",taxa.name.new))

#Scenario12
scen12 <- ggtree(tree2,layout="fan",open.angle=230) %<+% scenario12.dataA +
  geom_point(aes(color=mean.color, size = mean.size), alpha=1)+
  scale_size_manual(values = c("<= 0" = 1, "0 < Mean <= 1.0" = 2, "1.0 < Mean < 2.0" = 3, "Mean >= 2.0" = 4))+
  scale_color_manual(values = c("<= 0"="#9a9a9a", "0 < Mean <= 1.0"="#e69f00", "1.0 < Mean < 2.0"="#cc0000", "Mean >= 2.0"="#660000"))+
  
  # geom_point(aes(color=Mean.total, size = Mean.total), alpha=1)+
  # scale_color_gradient2(low = "gray", mid = c("gray","yellow","orange","#CC0066"),
  #                       high = "#CC0066", midpoint = 1.2, limits = c(min, max), na.value = NA) +
  theme(legend.position = "none")

#Scenario13
scen13 <- ggtree(tree2,layout="fan",open.angle=230) %<+% scenario13.dataA +
  geom_point(aes(color=mean.color, size = mean.size), alpha=1)+
  scale_size_manual(values = c("<= 0" = 1, "0 < Mean <= 1.0" = 2, "1.0 < Mean < 2.0" = 3, "Mean >= 2.0" = 4))+
  scale_color_manual(values = c("<= 0"="#9a9a9a", "0 < Mean <= 1.0"="#e69f00", "1.0 < Mean < 2.0"="#cc0000", "Mean >= 2.0"="#660000"))+
  
  # geom_point(aes(color=Mean.total, size = Mean.total), alpha=1)+
  # scale_color_gradient2(low = "gray", mid = c("gray","yellow","orange","#CC0066"), 
  #                       high = "#CC0066", midpoint = 1.2, limits = c(min, max), na.value = NA) +
  theme(legend.position = "none")

#Scenario14
scen14 <- ggtree(tree2,layout="fan",open.angle=230) %<+% scenario14.dataA +
  geom_point(aes(color=mean.color, size = mean.size), alpha=1)+
  scale_size_manual(values = c("<= 0" = 1, "0 < Mean <= 1.0" = 2, "1.0 < Mean < 2.0" = 3, "Mean >= 2.0" = 4))+
  scale_color_manual(values = c("<= 0"="#9a9a9a", "0 < Mean <= 1.0"="#e69f00", "1.0 < Mean < 2.0"="#cc0000", "Mean >= 2.0"="#660000"))+
  
  # geom_point(aes(color=Mean.total, size = Mean.total), alpha=1)+
  # scale_color_gradient2(low = "gray", mid = c("gray","yellow","orange","#CC0066"),
  #                       high = "#CC0066", midpoint = 1.2, limits = c(min, max), na.value = NA) +
  theme(legend.position = "none")

#Scenario15
scen15 <- ggtree(tree2,layout="fan",open.angle=230) %<+% scenario15.dataA +
  geom_point(aes(color=mean.color, size = mean.size), alpha=1)+
  scale_size_manual(values = c("<= 0" = 1, "0 < Mean <= 1.0" = 2, "1.0 < Mean < 2.0" = 3, "Mean >= 2.0" = 4))+
  scale_color_manual(values = c("<= 0"="#9a9a9a", "0 < Mean <= 1.0"="#e69f00", "1.0 < Mean < 2.0"="#cc0000", "Mean >= 2.0"="#660000"))+
  
  # geom_point(aes(color=Mean.total, size = Mean.total), alpha=1)+
  # scale_color_gradient2(low = "gray", mid = c("gray","yellow","orange","#CC0066"), 
  #                       high = "#CC0066", midpoint = 1.2, limits = c(min, max), na.value = NA) +
  theme(legend.position = "none")
scen15
ggsave(plot=scen15,fs::path(dir.results,"Dendrogram_Mixture_Scenario15.jpeg"))
# 
# #9. Individual Sharing Plot ----
# #New Individual Analysis ----
# means <- formatted.data %>%
#   tidylog::filter(Exposure!="Mixture")
# # tidylog::filter(Model!="Ridge ZING")
# means$Model <- factor(means$Model, levels = c("BaH-ZING", "ZING"))
# means$Scenario <- paste0("Scenario ",means$Scenario)
# # scenario_order <- c("Scenario 1","Scenario 2","Scenario 3","Scenario 4","Scenario 5",
# #                     "Scenario 6","Scenario 7","Scenario 8","Scenario 9","Scenario 10","Scenario 11",
# #                     "Scenario 12","Scenario 13","Scenario 14",
# #                     "Scenario 15")
# # means$Scenario <- factor(means$Scenario,levels=scenario_order)
# 
# taxa_order <- c("Species","Genus","Family","Order","Class","Phylum")
# means$Taxa.Level <- factor(means$Taxa.Level, levels = taxa_order)
# 
# #Filter to just a few replications for simplicity
# # means.hold <- means
# 
# means<-means %>%
#   filter(Replication==1 | Replication ==50 |Replication==150 |Replication==200|
#            Replication==30 | Replication==60|Replication==80 |
#            Replication==250|Replication==300|Replication==350|
#            Replication==400|Replication==450|
#            Replication==500|Replication==550|
#            Replication==600)
# 
# #Species & Genus Level
# names2 <- names %>% 
#   select(species,genus) %>% 
#   rename(Taxa=species,
#          group=genus)
# 
# species <- means %>% 
#   select(Taxa,Exposure,Model,Mean,Scenario,Taxa.Level,Replication,indicator) %>% 
#   filter(Taxa.Level=="Species") %>% 
#   select(-Taxa.Level)
# species <- tidylog::left_join(species,names2) 
# species <- species %>% 
#   rename(Species=Taxa) %>%
#   rename(Genus=group)
# genus <- means %>% 
#   filter(Taxa.Level=="Genus") %>% 
#   select(Taxa,Exposure,Model,Mean,Scenario,Replication,indicator) %>% 
#   rename(Genus=Taxa)
# species.genus <- tidylog::left_join(species,genus,by=c("Genus","Exposure","Model","Scenario","Replication"))
# species.genus <- species.genus %>% 
#   rename(Mean.species=Mean.x,
#          Mean.genus=Mean.y,
#          Indicator.species=indicator.x,
#          Indicator.genus=indicator.y)
# 
# #Add Family 
# names3 <- names %>% 
#   select(genus,family) %>% 
#   rename(Taxa=genus,
#          group=family) %>% 
#   distinct()
# 
# genus <- means %>% 
#   select(Taxa,Exposure,Model,Mean,Scenario,Taxa.Level,Replication,indicator) %>% 
#   filter(Taxa.Level=="Genus") %>% 
#   select(-Taxa.Level)
# genus <- tidylog::left_join(genus,names3) 
# genus <- genus%>% 
#   rename(Genus=Taxa) %>%
#   rename(Family=group)
# family <- means %>% 
#   filter(Taxa.Level=="Family") %>% 
#   select(Taxa,Exposure,Model,Mean,Scenario,Replication,indicator) %>% 
#   rename(Family=Taxa)
# genus.family <- tidylog::left_join(genus,family,by=c("Family","Exposure","Model","Scenario","Replication"))
# genus.family <- genus.family %>% 
#   rename(Mean.genus=Mean.x,
#          Mean.family=Mean.y,
#          Indicator.genus=indicator.x,
#          Indicator.family=indicator.y)
# 
# #Add order
# names4 <- names %>% 
#   select(family,order) %>% 
#   rename(Taxa=family,
#          group=order) %>% 
#   distinct()
# 
# family <- means %>% 
#   select(Taxa,Exposure,Model,Mean,Scenario,Taxa.Level,Replication,indicator) %>% 
#   filter(Taxa.Level=="Family") %>% 
#   select(-Taxa.Level)
# family <- tidylog::left_join(family,names4) 
# family <- family%>% 
#   rename(Family=Taxa) %>%
#   rename(Order=group)
# order <- means %>% 
#   filter(Taxa.Level=="Order") %>% 
#   select(Taxa,Exposure,Model,Mean,Scenario,Replication,indicator) %>% 
#   rename(Order=Taxa)
# family.order <- tidylog::left_join(family,order,by=c("Order","Exposure","Model","Scenario","Replication"))
# family.order <- family.order %>% 
#   rename(Mean.family=Mean.x,
#          Mean.order=Mean.y,
#          Indicator.family=indicator.x,
#          Indicator.order=indicator.y)
# 
# #Add Class
# names5 <- names %>% 
#   select(order,class) %>% 
#   rename(Taxa=order,
#          group=class) %>% 
#   distinct()
# 
# order <- means %>% 
#   select(Taxa,Exposure,Model,Mean,Scenario,Taxa.Level,Replication,indicator) %>% 
#   filter(Taxa.Level=="Order") %>% 
#   select(-Taxa.Level)
# order <- tidylog::left_join(order,names5) 
# order <- order%>% 
#   rename(Order=Taxa) %>%
#   rename(Class=group)
# class <- means %>% 
#   filter(Taxa.Level=="Class") %>% 
#   select(Taxa,Exposure,Model,Mean,Scenario,Replication,indicator) %>% 
#   rename(Class=Taxa)
# order.class <- tidylog::left_join(order,class,by=c("Class","Exposure","Model","Scenario","Replication"))
# order.class <- order.class %>% 
#   rename(Mean.order=Mean.x,
#          Mean.class=Mean.y,
#          Indicator.order=indicator.x,
#          Indicator.class=indicator.y)
# 
# #Add phylum
# names6 <- names %>% 
#   select(class,phylum) %>% 
#   rename(Taxa=class,
#          group=phylum) %>% 
#   distinct()
# 
# class <- means %>% 
#   select(Taxa,Exposure,Model,Mean,Scenario,Taxa.Level,Replication,indicator) %>% 
#   filter(Taxa.Level=="Class") %>% 
#   select(-Taxa.Level)
# class <- tidylog::left_join(class,names6) 
# class <- class%>% 
#   rename(Class=Taxa) %>%
#   rename(Phylum=group)
# phylum <- means %>% 
#   filter(Taxa.Level=="Phylum") %>% 
#   select(Taxa,Exposure,Model,Mean,Scenario,Replication,indicator) %>% 
#   rename(Phylum=Taxa)
# class.phylum <- tidylog::left_join(class,phylum,by=c("Phylum","Exposure","Model","Scenario","Replication"))
# class.phylum <- class.phylum %>% 
#   rename(Mean.class=Mean.x,
#          Mean.phylum=Mean.y,
#          Indicator.class=indicator.x,
#          Indicator.phylum=indicator.y)
# #Combine species.genus with genus.family
# taxa.data <- tidylog::left_join(species.genus,genus.family,by=c("Genus","Exposure","Model","Scenario","Replication","Mean.genus","Indicator.genus"))
# taxa.data <- tidylog::left_join(taxa.data,family.order,by=c("Family","Exposure","Model","Scenario","Replication","Mean.family","Indicator.family"))
# taxa.data <- tidylog::left_join(taxa.data,order.class,by=c("Order","Exposure","Model","Scenario","Replication","Mean.order","Indicator.order"))
# taxa.data <- tidylog::left_join(taxa.data,class.phylum,by=c("Class","Exposure","Model","Scenario","Replication","Mean.class","Indicator.class"))
# taxa.data2 <- taxa.data 
# # filter(Scenario=="Scenario 1") %>% 
# # filter(Model=="BaH-ZING")
# NC <- taxa.data2 %>%
#   filter(Indicator.species=="Non-Causal")
# C <- taxa.data2 %>% 
#   filter(Indicator.species=="Causal")
# C <- C %>%
#   mutate(Species="Species",
#          Genus="Genus",
#          Family="Family",
#          Order="Order",
#          Class="Class",
#          Phylum="Phylum")
# NC <- NC %>%
#   mutate(Species="Species",
#          Genus="Genus",
#          Family="Family",
#          Order="Order",
#          Class="Class",
#          Phylum="Phylum")
# max <- max(C$Mean.species,C$Mean.genus,C$Mean.family,C$Mean.order,C$Mean.class,C$Mean.phylum)
# # min <- min(C$Mean.species,C$Mean.genus,C$Mean.family,C$Mean.order,C$Mean.class,C$Mean.phylum)
# min <- min(NC$Mean.species,NC$Mean.genus,NC$Mean.family,NC$Mean.order,NC$Mean.class,NC$Mean.phylum)
# 
# 
# #add error bars
# 
# individual.plot <- ggplot() +
#   geom_segment(data = NC, aes(x = Species, xend = Genus, y = Mean.species, yend = Mean.genus, color = Mean.species), linewidth = 0.25,alpha=0.5) +
#   geom_segment(data = C, aes(x = Species, xend = Genus, y = Mean.species, yend = Mean.genus, color = Mean.species), linewidth = 0.25,alpha=0.5)+
#   
#   geom_segment(data = NC, aes(x = Genus, xend = Family, y = Mean.genus, yend = Mean.family, color = Mean.genus), linewidth = 0.25,alpha=0.5) +
#   geom_segment(data = C, aes(x = Genus, xend = Family, y = Mean.genus, yend = Mean.family, color = Mean.genus), linewidth = 0.25,alpha=0.5)+
#   
#   geom_segment(data = NC, aes(x = Family, xend = Order, y = Mean.family, yend = Mean.order, color = Mean.family), linewidth = 0.25,alpha=0.5) +
#   geom_segment(data = C, aes(x = Family, xend = Order, y = Mean.family, yend = Mean.order, color = Mean.family), linewidth = 0.25,alpha=0.5)+
#   
#   geom_segment(data = NC, aes(x = Order, xend = Class, y = Mean.order, yend = Mean.class, color = Mean.order), linewidth = 0.25,alpha=0.5) +
#   geom_segment(data = C, aes(x = Order, xend = Class, y = Mean.order, yend = Mean.class, color = Mean.order), linewidth = 0.25,alpha=0.5)+
#   
#   geom_segment(data = NC, aes(x = Class, xend = Phylum, y = Mean.class, yend = Mean.phylum, color = Mean.class), linewidth = 0.25,alpha=0.5) +
#   geom_segment(data = C, aes(x = Class, xend = Phylum, y = Mean.class, yend = Mean.phylum, color = Mean.class), linewidth = 0.25,alpha=0.5)+
#   
#   
#   geom_point(data = NC, aes(x = Species, y = Mean.species, color = Mean.species), size=0.25) +
#   geom_point(data = NC, aes(x = Genus, y = Mean.genus, color = Mean.genus), size=0.25) +
#   geom_point(data = NC, aes(x = Family, y = Mean.family, color = Mean.family), size=0.25) +
#   geom_point(data = NC, aes(x = Order, y = Mean.order, color = Mean.order), size=0.25) +
#   geom_point(data = NC, aes(x = Class, y = Mean.class, color = Mean.class), size=0.25) +
#   geom_point(data = NC, aes(x = Phylum, y = Mean.phylum, color = Mean.phylum), size=0.25) +
#   
#   geom_point(data = C, aes(x = Species, y = Mean.species, color = Mean.species), size=0.25) +  
#   geom_point(data = C, aes(x = Genus, y = Mean.genus, color = Mean.genus), size=0.25) +
#   geom_point(data = C, aes(x = Family, y = Mean.family, color = Mean.family), size=0.25) +
#   geom_point(data = C, aes(x = Order, y = Mean.order, color = Mean.order), size=0.25) +
#   geom_point(data = C, aes(x = Class, y = Mean.class, color = Mean.class), size=0.25) +
#   geom_point(data = C, aes(x = Phylum, y = Mean.phylum, color = Mean.phylum), size=0.25) +
#   geom_errorbar(data=NC,
#                 aes(Species,Mean.species,group = Species),
#                 stat = "summary",
#                 fun.data = function(x) {
#                   mean_val <- mean(x)
#                   sd_val <- sd(x)
#                   data.frame(y = mean_val, ymin = mean_val - sd_val, ymax = mean_val + sd_val)
#                 },
#                 width = 0.1,
#                 linewidth = 0.35,
#                 color = "black"
#   ) +
#   geom_errorbar(data=C,
#                 aes(Species,Mean.species,group = Species),
#                 stat = "summary",
#                 fun.data = function(x) {
#                   mean_val <- mean(x)
#                   sd_val <- sd(x)
#                   data.frame(y = mean_val, ymin = mean_val - sd_val, ymax = mean_val + sd_val)
#                 },
#                 width = 0.1,
#                 linewidth = 0.35,
#                 color = "black"
#   ) +
#   geom_point(data=NC,
#              aes(Species,Mean.species,group = Species),
#              stat = "summary",
#              fun = mean,   # Use the mean function to calculate the summary statistic
#              size = 1,     # Size of the points
#              shape = 16,   # Shape of the points
#              color = "black",  # Color of the points
#              fill = "black") +
#   geom_point(data=C,
#              aes(Species,Mean.species,group = Species),
#              stat = "summary",
#              fun = mean,   # Use the mean function to calculate the summary statistic
#              size = 1,     # Size of the points
#              shape = 16,   # Shape of the points
#              color = "black",  # Color of the points
#              fill = "black") +
#   geom_errorbar(data=NC,
#                 aes(Genus,Mean.genus,group = Genus),
#                 stat = "summary",
#                 fun.data = function(x) {
#                   mean_val <- mean(x)
#                   sd_val <- sd(x)
#                   data.frame(y = mean_val, ymin = mean_val - sd_val, ymax = mean_val + sd_val)
#                 },
#                 width = 0.1,
#                 linewidth = 0.35,
#                 color = "black"
#   ) +
#   geom_errorbar(data=C,
#                 aes(Genus,Mean.genus,group = Genus),
#                 stat = "summary",
#                 fun.data = function(x) {
#                   mean_val <- mean(x)
#                   sd_val <- sd(x)
#                   data.frame(y = mean_val, ymin = mean_val - sd_val, ymax = mean_val + sd_val)
#                 },
#                 width = 0.1,
#                 linewidth = 0.35,
#                 color = "black"
#   ) +
#   geom_point(data=NC,
#              aes(Genus,Mean.genus,group = Genus),
#              stat = "summary",
#              fun = mean,   # Use the mean function to calculate the summary statistic
#              size = 1,     # Size of the points
#              shape = 16,   # Shape of the points
#              color = "black",  # Color of the points
#              fill = "black") +
#   geom_point(data=C,
#              aes(Genus,Mean.genus,group = Genus),
#              stat = "summary",
#              fun = mean,   # Use the mean function to calculate the summary statistic
#              size = 1,     # Size of the points
#              shape = 16,   # Shape of the points
#              color = "black",  # Color of the points
#              fill = "black") +
#   geom_errorbar(data=NC,
#                 aes(Family,Mean.family,group = Family),
#                 stat = "summary",
#                 fun.data = function(x) {
#                   mean_val <- mean(x)
#                   sd_val <- sd(x)
#                   data.frame(y = mean_val, ymin = mean_val - sd_val, ymax = mean_val + sd_val)
#                 },
#                 width = 0.1,
#                 linewidth = 0.35,
#                 color = "black"
#   ) +
#   geom_errorbar(data=C,
#                 aes(Family,Mean.family,group = Family),
#                 stat = "summary",
#                 fun.data = function(x) {
#                   mean_val <- mean(x)
#                   sd_val <- sd(x)
#                   data.frame(y = mean_val, ymin = mean_val - sd_val, ymax = mean_val + sd_val)
#                 },
#                 width = 0.1,
#                 linewidth = 0.35,
#                 color = "black"
#   ) +
#   geom_point(data=NC,
#              aes(Family,Mean.family,group = Family),
#              stat = "summary",
#              fun = mean,   # Use the mean function to calculate the summary statistic
#              size = 1,     # Size of the points
#              shape = 16,   # Shape of the points
#              color = "black",  # Color of the points
#              fill = "black") +
#   geom_point(data=C,
#              aes(Family,Mean.family,group = Family),
#              stat = "summary",
#              fun = mean,   # Use the mean function to calculate the summary statistic
#              size = 1,     # Size of the points
#              shape = 16,   # Shape of the points
#              color = "black",  # Color of the points
#              fill = "black") +
#   geom_errorbar(data=NC,
#                 aes(Order,Mean.order,group = Order),
#                 stat = "summary",
#                 fun.data = function(x) {
#                   mean_val <- mean(x)
#                   sd_val <- sd(x)
#                   data.frame(y = mean_val, ymin = mean_val - sd_val, ymax = mean_val + sd_val)
#                 },
#                 width = 0.1,
#                 linewidth = 0.35,
#                 color = "black"
#   ) +
#   geom_errorbar(data=C,
#                 aes(Order,Mean.order,group = Order),
#                 stat = "summary",
#                 fun.data = function(x) {
#                   mean_val <- mean(x)
#                   sd_val <- sd(x)
#                   data.frame(y = mean_val, ymin = mean_val - sd_val, ymax = mean_val + sd_val)
#                 },
#                 width = 0.1,
#                 linewidth = 0.35,
#                 color = "black"
#   ) +
#   geom_point(data=NC,
#              aes(Order,Mean.order,group = Order),
#              stat = "summary",
#              fun = mean,   # Use the mean function to calculate the summary statistic
#              size = 1,     # Size of the points
#              shape = 16,   # Shape of the points
#              color = "black",  # Color of the points
#              fill = "black") +
#   geom_point(data=C,
#              aes(Order,Mean.order,group = Order),
#              stat = "summary",
#              fun = mean,   # Use the mean function to calculate the summary statistic
#              size = 1,     # Size of the points
#              shape = 16,   # Shape of the points
#              color = "black",  # Color of the points
#              fill = "black") +
#   
#   geom_errorbar(data=NC,
#                 aes(Class,Mean.class,group = Class),
#                 stat = "summary",
#                 fun.data = function(x) {
#                   mean_val <- mean(x)
#                   sd_val <- sd(x)
#                   data.frame(y = mean_val, ymin = mean_val - sd_val, ymax = mean_val + sd_val)
#                 },
#                 width = 0.1,
#                 linewidth = 0.35,
#                 color = "black"
#   ) +
#   geom_errorbar(data=C,
#                 aes(Class,Mean.class,group = Class),
#                 stat = "summary",
#                 fun.data = function(x) {
#                   mean_val <- mean(x)
#                   sd_val <- sd(x)
#                   data.frame(y = mean_val, ymin = mean_val - sd_val, ymax = mean_val + sd_val)
#                 },
#                 width = 0.1,
#                 linewidth = 0.35,
#                 color = "black"
#   ) +
#   geom_point(data=NC,
#              aes(Class,Mean.class,group = Class),
#              stat = "summary",
#              fun = mean,   # Use the mean function to calculate the summary statistic
#              size = 1,     # Size of the points
#              shape = 16,   # Shape of the points
#              color = "black",  # Color of the points
#              fill = "black") +
#   geom_point(data=C,
#              aes(Class,Mean.class,group = Class),
#              stat = "summary",
#              fun = mean,   # Use the mean function to calculate the summary statistic
#              size = 1,     # Size of the points
#              shape = 16,   # Shape of the points
#              color = "black",  # Color of the points
#              fill = "black") +
#   geom_errorbar(data=NC,
#                 aes(Phylum,Mean.phylum,group = Phylum),
#                 stat = "summary",
#                 fun.data = function(x) {
#                   mean_val <- mean(x)
#                   sd_val <- sd(x)
#                   data.frame(y = mean_val, ymin = mean_val - sd_val, ymax = mean_val + sd_val)
#                 },
#                 width = 0.1,
#                 linewidth = 0.35,
#                 color = "black"
#   ) +
#   geom_errorbar(data=C,
#                 aes(Phylum,Mean.phylum,group = Phylum),
#                 stat = "summary",
#                 fun.data = function(x) {
#                   mean_val <- mean(x)
#                   sd_val <- sd(x)
#                   data.frame(y = mean_val, ymin = mean_val - sd_val, ymax = mean_val + sd_val)
#                 },
#                 width = 0.1,
#                 linewidth = 0.35,
#                 color = "black"
#   ) +
#   geom_point(data=NC,
#              aes(Phylum,Mean.phylum,group = Phylum),
#              stat = "summary",
#              fun = mean,   # Use the mean function to calculate the summary statistic
#              size = 1,     # Size of the points
#              shape = 16,   # Shape of the points
#              color = "black",  # Color of the points
#              fill = "black") +
#   geom_point(data=C,
#              aes(Phylum,Mean.phylum,group = Phylum),
#              stat = "summary",
#              fun = mean,   # Use the mean function to calculate the summary statistic
#              size = 1,     # Size of the points
#              shape = 16,   # Shape of the points
#              color = "black",  # Color of the points
#              fill = "black") +
#   
#   # scale_color_gradient2(low = "gray", mid = c("yellow","orange","#CC0066"), high = "#CC0066",
#   #                       midpoint = max/2, limits = c(min, max), na.value = NA)+
#   theme_bw()+
#   facet_wrap("Model",ncol=1)+
#   # scale_color_gradientn(colors = my_colors, limits = my_limits)
#   # scale_color_gradientn(colors = c("#9a9a9a","#e69f00","#cc0000","#660000"), 
#   #                       values = c(min, 0, 0.2, 0.4, max))
#   # scale_color_gradientn(colors = c("#9a9a9a","yellow","orange","#CC0066"), 
#   #                       values = c(-0.2489804, 0, 0.2, 0.4, 0.6550573))
#   # scale_colour_gradientn(colours = c("red","yellow","green","lightblue","darkblue"),
#   #                        values = c(1.0,0.8,0.6,0.4,0.2,0)) 
#   scale_colour_gradientn(colours = c("#660000","#660000","#cc0000","#cc0000","#e69f00","#e69f00","#9a9a9a","#9a9a9a","#9a9a9a","#9a9a9a"),
#                          values = c(1.0,0.8,0.6,0.4,0.2,0)) +
#   scale_x_discrete(name = "Taxonomic Level", limits = c("Species","Genus","Family","Order","Class","Phylum"))+
#   theme(text=element_text(family="Times",size=15),
#         strip.text.x = element_text(face="bold",size=15),
#         axis.text.x = element_text(angle=45,hjust=1,face="bold",size=15),
#         axis.title.y=element_text(face="bold"))+
#   ylab("Fold Change")+
#   labs(color = "Fold Change") +
#   theme(legend.position="left")+
#   facet_grid(Scenario ~ Model)
# individual.plot
# ggsave(plot=individual.plot,fs::path(dir.results,"Sharing Plot Scenario 7.jpeg"),height=4,width=8)
