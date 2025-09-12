#A. Early Integration ----
# Set Color Palettes 
col_pal <- RColorBrewer::brewer.pal(n = 8, name = "Dark2")

sankey_colors <- matrix(c("exposure", "darkgray",
                          "lc1",      col_pal[1],
                          "lc2",      col_pal[2],
                          "lc3",      col_pal[3],
                          "lc4",      col_pal[4],
                          "layer1",   col_pal[1],
                          "layer2",   "#7570B3",
                          "layer3",   col_pal[3],
                          "layer4",   col_pal[4],
                          "Outcome",  col_pal[8],
                          "TRUE",     "#6372e0", # Blue
                          "FALSE",    "#d1d4ff", # Light grey
                          "pos_clus_to_out", "red", 
                          "neg_clus_to_out", "#e4e5f2"), 
                        byrow = TRUE, nrow = 14)

# Change to dataframe
colnames(sankey_colors) <- c("domain", "range")
sankey_colors <- as.data.frame(sankey_colors)
# lucid_fit1 <- fit1
sankey_early_integration <- function(lucid_fit1, text_size = 15) {
  # Get sankey dataframe ----
  # x <- lucid_fit1
  get_sankey_df <- function(x,
                            G_color = "dimgray", 
                            X_color = "#eb8c30",
                            Z_color = "#2fa4da", 
                            Y_color = "#afa58e", 
                            pos_link_color = "#67928b", 
                            neg_link_color = "#d1e5eb", 
                            fontsize = 10) {
    K <- x$K
    var.names <- x$var.names
    dimG <- length(var.names$Gnames)
    dimZ <- length(var.names$Znames)
    valueGtoX <- as.vector(t(x$res_Beta[, -1]))
    valueXtoZ <- as.vector(t(x$res_Mu))
    valueXtoY <- as.vector(x$res_Gamma$beta)[1:K]
    
    # GtoX
    GtoX <- data.frame(
      source = rep(x$var.names$Gnames, K), 
      target = paste0("Latent Cluster", 
                      as.vector(sapply(1:K, function(x) rep(x, dimG)))), 
      value = abs(valueGtoX), 
      group = as.factor(valueGtoX > 0))
    
    # XtoZ
    XtoZ <- data.frame(
      source = paste0("Latent Cluster", 
                      as.vector(sapply(1:K, 
                                       function(x) rep(x, dimZ)))), 
      target = rep(var.names$Znames, 
                   K), value = abs(valueXtoZ),
      group = as.factor(valueXtoZ > 
                          0))
    
    # subset top 25% of each omics layer
    top25<- XtoZ %>%
      filter(source == "Latent Cluster1") %>%
      mutate(omics = case_when(grepl("metabolite", target) ~ "Metabolomics",
                               !grepl("metabolite", target) ~ "Microbiome")) %>% 
      # grepl("miR", target) ~ "miRNA")) %>%
      group_by(omics) %>%
      arrange(desc(value)) %>%
      slice_max(value, n = 10) %>% 
      ungroup()
    
    XtoZ_sub<- XtoZ %>%
      filter(target %in% top25$target)
    
    
    # XtoY
    XtoY <- data.frame(source = paste0("Latent Cluster", 1:K), 
                       target = rep(var.names$Ynames, K), value = abs(valueXtoY), 
                       group = as.factor(valueXtoY > 0))
    links <- rbind(GtoX, XtoZ_sub, XtoY)
    # links <- rbind(GtoX, XtoZ, XtoY)
    
    nodes <- data.frame(
      name = unique(c(as.character(links$source), 
                      as.character(links$target))), 
      group = as.factor(c(rep("exposure",
                              dimG), rep("lc", K), rep("biomarker", nrow(XtoZ_sub)/2), "outcome")))
    # group = as.factor(c(rep("exposure", 
    # dimG), rep("lc", K), rep("biomarker", dimZ), "outcome")))
    ## the following two lines were used to exclude covars from the plot
    links <- links %>% filter(!grepl("cohort", source) & 
                                !grepl("age", source) & 
                                !grepl("fish", source) &
                                !grepl("sex", source))
    nodes <- nodes %>% filter(!grepl("cohort", name) &
                                !grepl("age", name) & 
                                !grepl("fish", name) &
                                !grepl("sex", name)) 
    
    links$IDsource <- match(links$source, nodes$name) - 1
    links$IDtarget <- match(links$target, nodes$name) - 1
    
    color_scale <- data.frame(
      domain = c("exposure", "lc", "biomarker", 
                 "outcome", "TRUE", "FALSE"), 
      range = c(G_color, X_color, 
                Z_color, Y_color, pos_link_color, neg_link_color))
    
    sankey_df = list(links = links, 
                     nodes = nodes)
    return(sankey_df)
  }
  # 1. Get sankey dataframes ----
  sankey_dat <- get_sankey_df(lucid_fit1)
  n_omics <- length(lucid_fit1$var.names$Znames)
  # link data
  links <- sankey_dat[["links"]] 
  # node data
  nodes <- sankey_dat[["nodes"]] 
  
  nodes1 <- nodes %>% 
    mutate(group = case_when(str_detect(name,"Cluster") ~ "lc",
                             # str_detect(name, "cg") ~ "CpG",
                             str_detect(name, "eGFR_log_scaled") ~ "outcome",
                             # str_detect(name, "pro") ~ "Prot",
                             str_detect(name, "met") ~ "Metabolomics",
                             name %in% omics2 ~ "Microbiome",
                             # str_detect(name, "mic") ~ "Mic",
                             # str_detect(name, "tc") ~ "TC",
                             # str_detect(name, "miR") ~ "miRNA",
                             str_detect(name, "ma_") ~ "exposure"))
  # name = ifelse(group == "exposure", "Hg",name))
  links1 <- links 
  # mutate(source = ifelse(source == "G1", "Hg",source))
  # 6. Plotly Version ----
  
  ## 6.1 Set Node Color Scheme: ----
  color_pal_sankey <- matrix(
    c("exposure", sankey_colors$range[sankey_colors$domain == "exposure"],
      "lc",       "#b3d8ff",
      "Microbiome",     sankey_colors$range[sankey_colors$domain == "layer1"],
      "Metabolomics",      sankey_colors$range[sankey_colors$domain == "layer2"],
      # "miRNA", sankey_colors$range[sankey_colors$domain == "layer3"],
      "outcome",  sankey_colors$range[sankey_colors$domain == "Outcome"]), 
    ncol = 2, byrow = TRUE) %>%
    as_tibble(.name_repair = "unique") %>% 
    janitor::clean_names() %>%
    dplyr::rename(group = x1, color = x2)
  
  # Add color scheme to nodes
  nodes_new_plotly <- nodes1 %>% 
    left_join(color_pal_sankey) %>%
    mutate(
      x = case_when(
        group == "exposure" ~ 0,
        str_detect(name, "Cluster") ~ 1/3,
        str_detect(group, "Metabolomics")|
          str_detect(group, "Microbiome")|
          # str_detect(name, "miR")|
          str_detect(name, "eGFR_log_scaled")~ 2/3
      ))
  
  nodes_new_plotly1 <- nodes_new_plotly %>%
    # Modify names of features for plotting
    dplyr::select(group, color, x, name)%>% 
    mutate(name = case_when(name == "Latent Cluster1" ~ "<b>Joint Omics\nProfile 0</b>",
                            name == "Latent Cluster2" ~ "<b>Joint Omics\nProfile 1</b>",
                            TRUE ~ name))
  
  
  ## 6.2 Get links for Plotly, set color ----
  links_new <- links1  %>%
    mutate(
      link_color = case_when(
        # Ref link color
        value == 0 ~     "#f3f6f4",
        # # Cluster 
        # str_detect(source, "Cluster1") &  group == TRUE  ~  "#706C6C",
        # str_detect(source, "Cluster1") &  group == FALSE ~  "#D3D3D3",
        # str_detect(source, "Cluster2") &  group == TRUE  ~  "#706C6C",
        # str_detect(source, "Cluster2") &  group == FALSE ~  "#D3D3D3",
        ##############
        # Exposure
        str_detect(source, "ma_") &  group == TRUE  ~  "red",
        str_detect(source, "ma_") &  group == FALSE  ~  "#D3D3D3",
        # Outcome
        str_detect(target, "eGFR_log_scaled") &  group == TRUE  ~  "red",
        str_detect(target, "eGFR_log_scaled") &  group == FALSE  ~  "#D3D3D3",
        # Metabolomics
        str_detect(target, "met") &  group == TRUE  ~  "#a64d79",
        str_detect(target, "met") &  group == FALSE ~  "#ead1dc",
        # Transcriptome
        target %in% omics2 &  group == TRUE  ~  "#38761d",
        target %in% omics2 &  group == FALSE ~  "#b6d7a8",
        # proteome
        # str_detect(target, "miR") &  group == TRUE  ~  "#a64d79",
        # str_detect(target, "miR") &  group == FALSE ~  "#ead1dc",
        ##
        group == FALSE ~ "#D3D3D3", # Negative association
        group == TRUE ~  "#706C6C")) # Positive association
  
  links_new1<- links_new %>%
    dplyr::select(colnames(links_new), target)
  
  plotly_link <- list(
    source = links_new1$IDsource,
    target = links_new1$IDtarget,
    value = links_new1$value+.00000000000000000000001, 
    color = links_new1$link_color)  
  
  # Get list of nodes for Plotly
  plotly_node <- list(
    label = nodes_new_plotly1$name, 
    color = nodes_new_plotly1$color,
    pad = 15,
    thickness = 20,
    line = list(color = "black",width = 0.5),
    x = nodes_new_plotly1$x, 
    # y = c(0.01, 
    #       0.3, 0.7, # clusters
    #       seq(from = .01, to = 1, by = 0.04)[1:(dimZ * 0.25)], # biomaker
    #       .95
    y = c(0.01,
          0.1, 0.5, # clusters
          seq(from = .05, to = 1, by = 0.04)[1:21],
          # seq(from = (.01+0.06*7), to = 1, by = 0.08)[1:5],
          # 0.9,
          # biomaker
          0.98
    ))
  
  
  ## 6.3 Plot Figure ----
  (fig <- plot_ly(
    type = "sankey",
    domain = list(
      x =  c(0,1),
      y =  c(0,1)),
    orientation = "h",
    node = plotly_node,
    link = plotly_link))
  
  (fig <- fig %>% layout(
    # title = "Basic Sankey Diagram",
    font = list(
      size = text_size
    ))
  )
  return(fig)
}


#Plot omics profiles
# fit <- fit1
# integration_type <- "early"
# omics_lst_data <- omics_lst
plot_omics_profiles <- function(fit, integration_type, omics_lst_data) {
  
  # Combines omics data into one dataframe
  omics_lst_df <- purrr::map(omics_lst_data, ~tibble::as_tibble(.x, rownames = "name"))
  
  # Get metadata file
  meta_df <- imap_dfr(omics_lst_df,
                      ~tibble(omic_layer = .y, ftr_name = names(.x))) |>
    filter(ftr_name != "name") |>
    mutate(omic_num = case_when(str_detect(omic_layer, "microbiome") ~ 1,
                                # str_detect(omic_layer, "transc") ~ 2,
                                # str_detect(omic_layer, "miR") ~ 3,
                                # str_detect(omic_layer, "pro") ~ 4,
                                str_detect(omic_layer, "metabolomics") ~ 2))
  
  if(integration_type == "Early"){
    M_mean = as.data.frame(fit$res_Mu)
    M_mean$cluster = as.factor(1:2)
    # Reshape the data
    M_mean_melt <- M_mean %>% 
      pivot_longer(cols = -cluster, names_to = "variable", values_to = "value")
    
    M_mean_melt <- M_mean_melt %>% 
      mutate(cluster = paste0("Cluster ", cluster))
    # add color label for omics layer
    M_mean_melt = M_mean_melt %>%
      mutate(color_label = case_when(str_detect(variable,  "_scaled") ~ "1", 
                                     str_detect(variable, "met") ~ "2", 
                                     TRUE ~ "3"))
    fig <- ggplot(M_mean_melt, 
                  aes(fill = color_label, y = value, x = variable)) +
      geom_bar(position="dodge", stat="identity") +
      ggtitle("Omics profiles for the two latent clusters") +
      facet_grid(rows = vars(cluster), scales = "free_y") +
      theme(legend.position="none") +
      geom_hline(yintercept = 0) +
      xlab("") +
      theme(text = element_text(size=10),
            axis.text.x = element_text(angle = 90, vjust = 1,
                                       hjust = 1),
            plot.margin = margin(10, 10, 10, 80),
            panel.background = element_rect(fill="white"), 
            strip.background = element_rect(fill = "white"),
            axis.line.x = element_line(color = "black"),
            axis.line.y = element_line(color = "black"),) +
      scale_fill_manual(values = c("#2fa4da", "#e7b6c1","#A77E69"))
    
  } else if(integration_type == "Intermediate"){
    M_mean = as_tibble(fit$res_Mu[[1]], rownames = "variable") %>%
      bind_rows(as_tibble(fit$res_Mu[[2]], rownames = "variable")) %>%
      bind_rows(as_tibble(fit$res_Mu[[3]], rownames = "variable"))
    
    # Reorder results because mirna order is reversed
    M_mean1 <- M_mean %>% 
      left_join(meta_df, by = c("variable" = "ftr_name")) %>%
      mutate(`Low Risk`  =  if_else(omic_layer == "miRna", V2, V1), 
             `High Risk` =  if_else(omic_layer == "miRna", V1, V2)) %>%
      dplyr::select(-c("V1", "V2"))
    
    # Pivot longer for figure 
    M_mean_l <- M_mean1 %>% 
      pivot_longer(cols = c(`Low Risk`, `High Risk`),
                   names_to = "cluster",
                   values_to = "value")
    
    # add color label for omics layer
    M_mean2 = M_mean_l %>%
      mutate(color_label = case_when(omic_layer == "methylome" ~ "1", 
                                     omic_layer == "transcriptome" ~ "2", 
                                     omic_layer == "miRna" ~ "3"), 
             low_high = if_else(str_detect(cluster, "Low"), 0,1),
             omic = if_else(omic_layer == "miRna", 
                            "miR",
                            str_sub(omic_layer, end = 1) %>% toupper()),
             omic_cluster = str_c(omic, low_high))
    
    # Filter only the top ## differential expressed features 
    M_mean2_top <- M_mean2 %>% 
      group_by(variable) %>% 
      filter(abs(value) == max(abs(value))) %>% 
      ungroup() %>% 
      arrange(max(abs(value))) %>% 
      group_by(omic_layer) %>% 
      slice_head(n=12) %>%
      ungroup()
    
    # Plots top 12 features
    fig <- ggplot(M_mean2  %>% filter(variable %in% M_mean2_top$variable),
                  aes(fill = color_label, y = value, x = variable)) +
      geom_bar(position="dodge", stat="identity") +
      ggtitle("Omics profiles for 2 latent clusters - Lucid in Parallel") +
      facet_grid(rows = vars(cluster),
                 cols = vars(omic_layer), scales = "free_x", space = "free") +
      theme(legend.position="none") +
      geom_hline(yintercept = 0) +
      xlab("") +
      theme(text = element_text(size=10),
            axis.text.x = element_text(angle = 90, vjust = 1,
                                       hjust = 1),
            plot.margin = margin(10, 10, 10, 80),
            panel.background = element_rect(fill="white"),
            strip.background = element_rect(fill = "white"),
            axis.line.x = element_line(color = "black"),
            axis.line.y = element_line(color = "black"),) +
      scale_fill_manual(values = c("#2fa4da", "#A77E69", "#e7b6c1"))
  }
  
  return(fig)
}


#B. Intermediate Integration ----
reorder_lucid_parallel <- function(lucidus_fit,
                                   reference = NULL) {
  if(is.null(reference)) {
    warning("no reference specified, return the original model")
    return(lucidus_fit)
  }
  
  n_omic <- length(reference)
  
  # reorder beta
  GtoX <- lucidus_fit$res_Beta$Beta
  lucidus_fit$res_Beta$Beta <- lapply(1:n_omic, function(i) {
    (-1)^(reference[i] - 1) * GtoX[[i]] # if reference = 1, no changes; 
    # if reference = 2, flip the reference and negate the estimates
  })
  # reorder mu
  XtoZ <- lucidus_fit$res_Mu
  lucidus_fit$res_Mu <- lapply(1:n_omic, function(i) {
    x <- c(1, 2) # order of clusters
    if(reference[i] == 2) {
      x <- c(2, 1)
      XtoZ[[i]][, x]
    } else{
      XtoZ[[i]][, x]
    }
  }) 
  # reorder gamma
  XtoY <- lucidus_fit$res_Gamma$Gamma$mu
  XtoY[1] <- XtoY[1] + sum(XtoY[-1] * (reference - 1)) # reference level using the new reference
  XtoY[-1] <- (-1)^(reference - 1) * XtoY[-1] # if reference = 2, flip the estimates
  lucidus_fit$res_Gamma$Gamma$mu <- XtoY
  lucidus_fit$res_Gamma$fit$coefficients <- XtoY
  
  # return the object using the new reference
  return(lucidus_fit)
}


# Set Color Palettes 
col_pal <- RColorBrewer::brewer.pal(n = 8, name = "Dark2")

# Set Sankey Colors ----
# Color pallet for sankey diagrams
sankey_colors <- matrix(c("exposure", col_pal[6],
                          "lc1",      col_pal[1],
                          "lc2",      col_pal[2],
                          "lc3",      col_pal[3],
                          "lc4",      col_pal[4],
                          "layer1",   col_pal[1],
                          "layer2",   col_pal[2],
                          "layer3",   col_pal[3],
                          "layer4",   col_pal[4],
                          "Outcome",  col_pal[8],
                          "TRUE",     "#6372e0", # Blue
                          "FALSE",    "#d1d4ff", # Light grey
                          "pos_clus_to_out", "red",
                          "neg_clus_to_out", "#e4e5f2"),
                        byrow = TRUE, nrow = 14)

# Change to dataframe
colnames(sankey_colors) <- c("domain", "range")
sankey_colors <- as.data.frame(sankey_colors)



#' Plot Sankey Diagram for LUCID in parallel
# lucidus_fit <- fit_reordered
plot_lucid_in_parallel_plotly<- function(lucidus_fit,
                                         sankey_colors,
                                         text_size = 10, 
                                         n_z_ftrs_to_plot = NULL){
  # Get number of clusters, layers, etc.
  K <- lucidus_fit$K
  dimG <- lucidus_fit$res_Beta$Beta[[1]] %>% ncol()-1
  n_layers   <- length(lucidus_fit$res_Beta$Beta)
  
  # Get top omics features based on effect size
  if(!is.null(n_z_ftrs_to_plot)){
    top_ftrs <- vector("list", n_layers)
    for(i in seq_along(top_ftrs)){
      top_ftrs[[i]] <- lucidus_fit$res_Mu[[i]] %>%
        as.data.frame() %>%
        rownames_to_column("name") %>%
        mutate(effect_size = abs(V1) + abs(V2)) %>%
        arrange(desc(effect_size))
      top_ftr_nms <- top_ftrs[[i]]$name[1:n_z_ftrs_to_plot[i]]
      lucidus_fit$res_Mu[[i]] <- 
        lucidus_fit$res_Mu[[i]][
          rownames(lucidus_fit$res_Mu[[i]]) %in% top_ftr_nms, ]
    }
  }
  
  mu_lst <- purrr::map(lucidus_fit$res_Mu, 
                       ~as.data.frame(.x) %>%
                         rownames_to_column("name"))
  names(mu_lst) <- paste0("layer", c(1:n_layers))
  dimZ <- purrr::map(mu_lst, ncol) %>% as.numeric()-1
  n_features <- purrr::map(mu_lst, nrow) %>% as.numeric()
  names(n_features) <- paste0("layer", c(1:n_layers))
  # Names of features and set order of omics features
  names_features <- bind_rows(mu_lst, .id = "color_group") %>% 
    rowwise() %>%
    mutate(sum = sum(abs(V1)+abs(V2)), 
           pos_c2 = if_else(V2>0, "pos", "neg")) %>%
    group_by(color_group, pos_c2) %>% arrange(-sum, .by_group = TRUE) %>% ungroup() %>% 
    mutate(rnum = row_number()) %>%
    group_by(name) %>% slice_head() %>% ungroup() %>%
    arrange(color_group, rnum) %>%
    dplyr::select(name, color_group)
  
  # Values for g --> x association
  valueGtoX <- c(lapply(lucidus_fit$res_Beta$Beta, 
                        function(x)(x[-1])) %>%
                   unlist(), 
                 rep(0, dimG*n_layers))
  
  # For Cluster 2 (which needs effect estimates): 
  valueGtoX_c1 <- do.call(rbind, lucidus_fit$res_Beta$Beta)[,-1] %>%
    as_tibble() %>%
    dplyr::mutate(layer = str_c("(Layer ", row_number(), ")"),
                  cluster = "Cluster 2") 
  
  # For cluster 1 (ref. cluster, effect est = 0):
  valueGtoX_c2 <- valueGtoX_c1 %>%
    mutate(across(where(is.numeric), ~0), 
           cluster = "Cluster 1")
  
  # combine, pivot longer, and create source and target columns
  GtoX <- bind_rows(valueGtoX_c1, valueGtoX_c2) %>%
    mutate(target = str_c(cluster, layer, sep = " ")) %>%
    pivot_longer(cols = setdiff(colnames(valueGtoX_c1), 
                                c("layer", "cluster")), 
                 names_to = "source", values_to = "value") %>%
    mutate(color_group = as.factor(value > 0), 
           value = abs(value)) %>%
    dplyr::select(source, target, value, color_group) %>%
    as.data.frame()
  
  valueXtoZ <- c(lapply(lucidus_fit$res_Mu, 
                        function(x)x[, 1]) %>% 
                   unlist(), 
                 lapply(lucidus_fit$res_Mu, 
                        function(x)x[, 2]) %>% 
                   unlist())
  
  valueXtoY <- c(rep(0, n_layers), 
                 # rep(lucidus_fit$res_Delta$Delta$mu[1] / n_layers, n_layers),
                 lucidus_fit$res_Gamma$Gamma$mu[-1])
  
  # n features in each layer
  XtoZ <- data.frame(source = c(rep("Cluster 1 (Layer 1)", n_features[1]),
                                rep("Cluster 1 (Layer 2)", n_features[2]),
                                # rep("Cluster 1 (Layer 3)", n_features[3]),
                                # rep("Cluster 1 (Layer 4)", n_features[4]),
                                rep("Cluster 2 (Layer 1)", n_features[1]),
                                rep("Cluster 2 (Layer 2)", n_features[2])
                                # rep("Cluster 2 (Layer 3)", n_features[3]) 
                                # rep("Cluster 2 (Layer 4)", n_features[4])
  ), 
  target = rep(c(lapply(lucidus_fit$res_Mu,
                        rownames) %>% unlist()),
               K[1]), 
  value = abs(valueXtoZ), 
  color_group = as.factor(valueXtoZ > 0))
  
  # To change the outcome from left to right hand side, flip source and target
  XtoY <- data.frame(target = rep("Outcome", 2*n_layers), 
                     source = c("Cluster 1 (Layer 1)", 
                                "Cluster 1 (Layer 2)",
                                # "Cluster 1 (Layer 3)",
                                # "Cluster 1 (Layer 4)",
                                "Cluster 2 (Layer 1)",
                                "Cluster 2 (Layer 2)"
                                # "Cluster 2 (Layer 3)" 
                                # "Cluster 2 (Layer 4)"
                     ), 
                     value = abs(valueXtoY), 
                     color_group = as.factor(valueXtoY > 0))
  
  # create Sankey diagram
  # Create Links ----
  links <- rbind(GtoX, XtoZ, XtoY) %>%
    mutate(
      # Group: one of exposure, clusters, or outcomes 
      # (doesn't include Z.order by desired order)
      source_group = case_when(
        str_detect(source, "Cluster") ~ "2_Cluster", 
        # source == "Outcome" ~ "3_outcome", # removed when moving outcome to right
        TRUE ~ "1_exposure"), 
      # Source Omics Layer: lc1-lc4 (for omics layers), outcome, or other 
      source_layer = case_when(
        str_detect(source, "Layer 1") ~ "lc1", 
        str_detect(source, "Layer 2") ~ "lc2", 
        # str_detect(source, "Layer 3") ~ "lc3", 
        # str_detect(source, "Layer 4") ~ "lc4",
        # source == "Outcome" ~ str_sub(target, start = -3, end = -2),  # removed when moving outcome to right
        TRUE ~ "exposure"), 
      # Source group_ for color (one of: exposure, : lc1-lc4 (for omics layers), outcome, or other 
      color_group_node = if_else(source == "Outcome", 
                                 "Outcome", 
                                 source_layer)) %>%
    group_by(source_group) %>%
    arrange(source_layer, .by_group = TRUE) %>%
    ungroup() %>%
    dplyr::select(source, target, value, color_group, color_group_node)
  
  
  # Create Nodes ----
  nodes <- links %>%
    dplyr::select(source, color_group_node) %>%
    mutate(rownum = row_number()) %>%
    rename(name = source, 
           color_group = color_group_node) %>%
    # Add outcome (only if outcome is on right side)
    bind_rows(data.frame(name = "Outcome", color_group = "Outcome")) %>%
    group_by(name) %>%
    slice_head() %>%
    ungroup() %>%
    arrange(rownum) %>%
    dplyr::select(-rownum) %>%
    # Add feature names
    bind_rows(names_features) %>% 
    mutate(id = row_number()-1) %>%
    left_join(sankey_colors, by = c( "color_group"= "domain"))
  
  # Join links and nodes for color names -----
  links <- links %>%
    left_join(nodes %>% 
                dplyr::select(id, name), 
              by = c("source" = "name")) %>%
    rename(source_id = id) %>% 
    dplyr::select(source_id, everything()) %>%
    left_join(nodes %>% 
                dplyr::select(id, name), 
              by = c("target" = "name")) %>%
    rename(target_id = id) %>% 
    dplyr::select(source_id,source, target_id,everything()) 
  
  
  # Manually change colors ----
  links <- links  %>%
    mutate(
      link_color = case_when(
        # Ref link color
        value == 0 ~     "#f3f6f4", 
        # Outcome
        str_detect(target, "Outcome") &  color_group == TRUE  ~  "red",
        str_detect(target, "Outcome") &  color_group == FALSE  ~  "#d9d2e9",
        # Methylation 
        str_detect(source, "Layer 2") &  color_group == TRUE  ~  "#bf9000",
        str_detect(source, "Layer 2") &  color_group == FALSE ~  "#ffd966",
        # Transcriptome
        str_detect(source, "Layer 1") &  color_group == TRUE  ~  "#38761d",
        str_detect(source, "Layer 1") &  color_group == FALSE ~  "#b6d7a8",
        # mirna
        # str_detect(source, "Layer 3") &  color_group == TRUE  ~  "#a64d79",
        # str_detect(source, "Layer 3") &  color_group == FALSE ~  "#ead1dc",
        
        links$color_group == FALSE ~ "#d9d2e9", # Negative association
        links$color_group == TRUE ~  "red"))
  
  ## change node names
  nodes <- nodes %>%
    mutate(name = case_when(name == "value" ~ "<b>Hg</b>",
                            name == "Cluster 1 (Layer 1)" ~ "<b>Microbiome\nProfile 0</b>",
                            name == "Cluster 2 (Layer 1)" ~ "<b>Microbiome\nProfile 1</b>",
                            name == "Cluster 1 (Layer 2)" ~ "<b>Metabolome\nProfile 0</b>",
                            name == "Cluster 2 (Layer 2)" ~ "<b>Metabolome\nProfile 1</b>",
                            TRUE ~ name))
  nodes <- nodes %>% 
    mutate(x = case_when(color_group=="exposure" ~ 0,
                         str_detect(name, "Metabolome") | str_detect(name,"Microbiome")|
                           color_group=="layer1" | color_group=="layer2" ~ 1/3, 
                         str_detect(name, "Outcome") ~ 2/3))
  
  (fig <- plot_ly(
    type = "sankey",
    orientation = "h",
    domain = list(
      x =  c(0,0.8),
      y =  c(0,1)),
    # arrangement = "snap",
    node = list(
      label = nodes$name,
      color = nodes$range,
      pad = 15,
      thickness = 20,
      line = list(
        color = "black",
        width = 0.5
      ),
      x = nodes$x
    ),
    
    link = list(
      source = links$source_id,
      target = links$target_id,
      value =  links$value+.00000000000000000000001,
      # label = links$source,
      color = links$link_color
    )
  )
  )
  
  fig <- fig %>% layout(
    font = list(
      size = text_size
    ),
    xaxis = list(showgrid = F, zeroline = F),
    yaxis = list(showgrid = F, zeroline = F)
  )
  
  fig
}

#C. Late Integration ----
#' Plot Sankey Diagram for LUCID in late integration
#' # Set Color Palettes 
col_pal <- RColorBrewer::brewer.pal(n = 8, name = "Dark2")

# Set Sankey Colors ----
# Color pallet for sankey diagrams
sankey_colors <- matrix(c("exposure", col_pal[6],
                          "lc1",      col_pal[1],
                          "lc2",      col_pal[2],
                          "lc3",      col_pal[3],
                          "lc4",      col_pal[4],
                          "layer1",   col_pal[1],
                          "layer2",   col_pal[2],
                          "layer3",   col_pal[3],
                          "layer4",   col_pal[4],
                          "Outcome",  col_pal[8],
                          "TRUE",     "#6372e0", # Blue
                          "FALSE",    "#d1d4ff", # Light grey
                          "pos_clus_to_out", "red",
                          "neg_clus_to_out", "#e4e5f2"),
                        byrow = TRUE, nrow = 14)

col_pal <- RColorBrewer::brewer.pal(n = 8, name = "Dark2")
color_pal_sankey <- matrix(c("exposure", col_pal[6],
                             "lc"      , "#b3d8ff",
                             "Microbiome","#1B9E77",
                             "CpG"     , "#1B9E77",
                             "Metabolome","#7570B3",
                             "outcome" , "#D95F02"),
                           ncol = 2, byrow = TRUE) %>%
  as_tibble(.name_repair = "unique") %>%
  janitor::clean_names() %>%
  dplyr::rename(group = x1, color = x2)
# x <- fit2
# Get sankey dataframe
get_sankey_df <- function(x,
                          G_color = "dimgray", 
                          X_color = "#eb8c30",
                          Z_color = "#2fa4da", 
                          Y_color = "#afa58e", 
                          pos_link_color = "#67928b", 
                          neg_link_color = "#d1e5eb", 
                          fontsize = 7) {
  K <- x$K
  var.names <- x$var.names
  pars <- x$pars
  dimG <- length(var.names$Gnames)
  dimZ <- length(var.names$Znames)
  valueGtoX <- as.vector(t(x$res_Beta[, -1]))
  valueXtoZ <- as.vector(t(x$res_Mu))
  valueXtoY <- as.vector(x$res_Gamma$beta)[1:K]
  
  # GtoX
  GtoX <- data.frame(
    source = rep(x$var.names$Gnames, K), 
    target = paste0("Latent Cluster", 
                    as.vector(sapply(1:K, function(x) rep(x, dimG)))), 
    value = abs(valueGtoX), 
    group = as.factor(valueGtoX > 0))
  
  # XtoZ
  XtoZ <- data.frame(
    source = paste0("Latent Cluster", 
                    as.vector(sapply(1:K, 
                                     function(x) rep(x, dimZ)))), 
    target = rep(var.names$Znames, 
                 K), value = abs(valueXtoZ),
    group = as.factor(valueXtoZ > 
                        0))
  # XtoY
  XtoY <- data.frame(source = paste0("Latent Cluster", 1:K), 
                     target = rep(var.names$Ynames, K), value = abs(valueXtoY), 
                     group = as.factor(valueXtoY > 0))
  
  links <- rbind(GtoX, XtoZ, XtoY)
  
  nodes <- data.frame(
    name = unique(c(as.character(links$source), 
                    as.character(links$target))), 
    group = as.factor(c(rep("exposure", 
                            dimG), rep("lc", K), rep("biomarker", dimZ), "outcome")))
  
  ## the following two lines were used to exclude covars from the plot 
  links <- links %>% filter(!grepl("cohort", source) & 
                              !grepl("age", source) & 
                              !grepl("fish", source) &
                              !grepl("sex", source))
  nodes <- nodes %>% filter(!grepl("cohort", name) &
                              !grepl("age", name) & 
                              !grepl("fish", name) &
                              !grepl("sex", name))  
  
  
  links$IDsource <- match(links$source, nodes$name) - 1
  links$IDtarget <- match(links$target, nodes$name) - 1
  
  color_scale <- data.frame(
    domain = c("exposure", "lc", "biomarker", 
               "outcome", "TRUE", "FALSE"), 
    range = c(G_color, X_color, 
              Z_color, Y_color, pos_link_color, neg_link_color))
  
  sankey_df = list(links = links, 
                   nodes = nodes)
  
  # p <- sankeyNetwork(
  #   Links = sankey_df$links, 
  #   Nodes = sankey_df$nodes, 
  #   Source = "IDsource", 
  #   Target = "IDtarget",
  #   Value = "value", 
  #   NodeID = "name", 
  #   colourScale = JS(sprintf("d3.scaleOrdinal()\n .domain(%s)\n .range(%s)\n ", 
  #                            jsonlite::toJSON(color_scale$domain), 
  #                            jsonlite::toJSON(color_scale$range))), 
  #   LinkGroup = "group", 
  #   NodeGroup = "group", 
  #   sinksRight = FALSE, 
  #   fontSize = fontsize)
  # p
  return(sankey_df)
}


# sankey_in_serial Function ----
# lucid_fit1 <- fit1
# lucid_fit2 <- fit2
sankey_in_serial <- function(lucid_fit1, lucid_fit2, color_pal_sankey, text_size = 15) {
  
  # 1. Get sankey dataframes ----
  sankey_dat1 <- get_sankey_df(lucid_fit1)
  sankey_dat2 <- get_sankey_df(lucid_fit2)
  # sankey_dat3 <- get_sankey_df(lucid_fit3)
  
  n_omics_1 <- length(lucid_fit1$var.names$Znames)
  n_omics_2 <- length(lucid_fit2$var.names$Znames)
  # n_omics_3 <- length(lucid_fit3$var.names$Znames)
  
  # combine link data
  lnks1_microbiome <- sankey_dat1[["links"]] %>% mutate(analysis = "1_microbiome")
  lnks2_metabolome  <- sankey_dat2[["links"]] %>% mutate(analysis = "2_metabolome")
  # lnks3_transcription    <- sankey_dat3[["links"]] %>% mutate(analysis = "3_transcript")
  links <- bind_rows(lnks1_microbiome, lnks2_metabolome)
  
  # combine node data
  nodes1_microbiome <- sankey_dat1[["nodes"]] %>% mutate(analysis = "1_microbiome")
  nodes2_metabolome  <- sankey_dat2[["nodes"]] %>% mutate(analysis = "2_metabolome")
  # nodes3_transcription    <- sankey_dat3[["nodes"]] %>% mutate(analysis = "3_transcript")
  nodes <- bind_rows(nodes1_microbiome, nodes2_metabolome)
  
  
  # 2. Modify analysis 1 ----
  # For analysis 1, latent clusters need to be renamed to names from analysis 2:
  ## 2.1 Get new and original latent cluster names (from the next analysis) ----
  names_clusters_1 <- data.frame(
    name_og = c("Latent Cluster1", "Latent Cluster2"), 
    name_new = c("<b>Microbiome\nProfile 0</b>", "<b>Microbiome\nProfile 1</b>"))
  
  ## 2.2 Change link names ----
  # Change link names and 
  lnks1_microbiome_new <- sankey_dat1[["links"]] %>%
    mutate(
      analysis = "1_microbiome",
      source = case_when(
        source == names_clusters_1$name_og[1] ~ names_clusters_1$name_new[1],
        source == names_clusters_1$name_og[2] ~ names_clusters_1$name_new[2],
        TRUE ~ source),
      target = case_when(
        target == names_clusters_1$name_og[1] ~ names_clusters_1$name_new[1],
        target == names_clusters_1$name_og[2] ~ names_clusters_1$name_new[2],
        TRUE ~ target)) %>%
    filter(target != "outcome")
  
  ## 2.3 Change node names ----
  # first, change latent cluster names to analysis specific cluster names
  nodes1_microbiome_new <- sankey_dat1[["nodes"]] %>%
    mutate(
      name = case_when(
        name == names_clusters_1$name_og[1] ~ names_clusters_1$name_new[1],
        name == names_clusters_1$name_og[2] ~ names_clusters_1$name_new[2],
        TRUE ~ name), 
      group = if_else(group == "biomarker", "CpG", as.character(group))) %>%
    filter(group != "outcome")
  
  
  # Visualize
  # sankeyNetwork(
  #   Links = lnks1_microbiome_new,
  #   Nodes = nodes1_microbiome_new,
  #   Source = "IDsource", Target = "IDtarget",
  #   Value = "value", NodeID = "name", LinkGroup = "group", NodeGroup = "group",
  #   sinksRight = FALSE)
  
  
  # 3. Modify analysis 2 ----
  # For analysis 2, latent clusters need to be renamed to names from analysis 3:
  ## 3.1 Get new and og latent cluster names ----
  names_clusters_2 <- data.frame(
    name_og = c("Latent Cluster1", "Latent Cluster2"), 
    name_new = c("<b>Metabolome\nProfile 0</b>", "<b>Metabolome\nProfile 1</b>"))
  
  ## 3.2 Change cluster names ----
  lnks2_metabolome_new <- sankey_dat2[["links"]] %>% 
    mutate(
      analysis = "2_metabolome", 
      source = case_when(
        source == names_clusters_2$name_og[1] ~ names_clusters_2$name_new[1], 
        source == names_clusters_2$name_og[2] ~ names_clusters_2$name_new[2], 
        TRUE ~ source), 
      target = case_when(
        target == names_clusters_2$name_og[1] ~ names_clusters_2$name_new[1], 
        target == names_clusters_2$name_og[2] ~ names_clusters_2$name_new[2], 
        TRUE ~ target)) %>%
    filter(target != "outcome")
  # filter(target != "outcome")
  
  ## 3.3 Change node names ----
  nodes2_metabolome_new <- sankey_dat2[["nodes"]] %>% 
    mutate(
      name = case_when(
        name == names_clusters_2$name_og[1] ~ names_clusters_2$name_new[1], 
        name == names_clusters_2$name_og[2] ~ names_clusters_2$name_new[2], 
        TRUE ~ name), 
      group = case_when(group == "exposure" ~ "lc", 
                        group == "biomarker" ~ "Metabolome",
                        TRUE ~ as.character(group))) %>%
    filter(name != "outcome")
  
  # Visualize
  # sankeyNetwork(
  #   Links = lnks2_transcript_new, 
  #   Nodes = nodes2_transcript_new,
  #   Source = "IDsource", Target = "IDtarget",
  #   Value = "value", NodeID = "name", 
  #   LinkGroup = "group", NodeGroup = "group",
  #   sinksRight = FALSE)
  ##
  
  # # 4. Modify analysis 3 ----
  # # For analysis 2, latent clusters need to be renamed to names from analysis 3:
  # ## 4.1 Get new and og latent cluster names ----
  # names_clusters_3 <- tibble(
  #   name_og = c("Latent Cluster1", "Latent Cluster2"),
  #   name_new = c("<b>Transcriptome\nProfile 0</b>", "<b>Transcriptome\nProfile 1</b>")) 
  # 
  # 
  # ## 4.2 Change cluster names ----
  # lnks3_transcript_new <- sankey_dat3[["links"]] %>% 
  #   mutate(
  #     analysis = "3_transcript", 
  #     source = case_when(
  #       source == names_clusters_3$name_og[1] ~ names_clusters_3$name_new[1], 
  #       source == names_clusters_3$name_og[2] ~ names_clusters_3$name_new[2], 
  #       TRUE ~ source), 
  #     target = case_when(
  #       target == names_clusters_3$name_og[1] ~ names_clusters_3$name_new[1], 
  #       target == names_clusters_3$name_og[2] ~ names_clusters_3$name_new[2], 
  #       TRUE ~ target))
  # 
  # ## 4.3 Change node names ----
  # nodes3_transcript_new <- sankey_dat3[["nodes"]] %>% 
  #   mutate(
  #     name = case_when(
  #       name == names_clusters_3$name_og[1] ~ names_clusters_3$name_new[1], 
  #       name == names_clusters_3$name_og[2] ~ names_clusters_3$name_new[2], 
  #       TRUE ~ name), 
  #     group = case_when(group == "exposure" ~ "lc", 
  #                       group == "biomarker" ~ "TC",
  #                       TRUE ~ as.character(group)))
  # 
  # Test/Visualize
  # sankeyNetwork(
  #   Links = lnks3_protein_new, 
  #   Nodes = nodes3_protein_new,
  #   Source = "IDsource", Target = "IDtarget",
  #   Value = "value", NodeID = "name", LinkGroup = "group", NodeGroup = "group",
  #   sinksRight = FALSE)
  
  
  
  # 5. Combine analysis 1-3 ----
  
  ## 5.1 Final Links ----
  links_all_1 <- bind_rows(lnks1_microbiome_new, 
                           lnks2_metabolome_new) %>% 
    # lnks3_transcript_new) %>%
    dplyr::select(-IDsource, -IDtarget)
  
  
  ### 5.1.1 Arrange by magnitude ----
  omics_priority <- links_all_1 %>% 
    filter(str_detect(source, "Profile 0"), 
           str_detect(target, "Profile 0", negate = TRUE), 
           str_detect(target, "Profile 1", negate = TRUE), 
           str_detect(target, "eGFR_log_scaled", negate = TRUE)) %>%
    group_by(source) %>%
    arrange(desc(group), desc(value), .by_group = TRUE) %>%
    mutate(omics_order = row_number()) %>%
    ungroup() %>%
    dplyr::select(target, omics_order)
  
  
  
  links_all <- links_all_1 %>%
    left_join(omics_priority) %>%
    mutate(
      # arrange_me = if_else(is.na(omics_order), 
      #                           "dont_arrange", 
      #                           "arrange"), 
      row_num = row_number(), 
      # row_num_order_comb = if_else(is.na(omics_order), 
      #                              row_num, 
      #                              omics_order), 
      row_num_to_add = if_else(is.na(omics_order), 
                               as.numeric(row_num), 
                               NA_real_) %>%
        zoo::na.locf(),
      order = if_else(is.na(omics_order), 
                      row_num_to_add, 
                      row_num_to_add+omics_order)
    ) %>%
    arrange(order)
  
  
  ### 5.1.2 Get new source and target IDs ----
  # First, combine all layers, get unique identifier
  node_ids <- tibble(name = unique(c(unique(links_all$source), 
                                     unique(links_all$target)))) %>%
    mutate(ID = row_number()-1)
  
  # Then combine with original data 
  links_new <- links_all %>%
    left_join(node_ids, by = c("source" = "name")) %>%
    dplyr::rename(IDsource = ID) %>%
    left_join(node_ids, by = c("target" = "name")) %>%
    dplyr::rename(IDtarget = ID)
  
  
  ## 5.2 Final Nodes ----
  nodes_new <- node_ids %>%
    dplyr::select(name) %>%
    left_join(bind_rows(nodes1_microbiome_new, 
                        nodes2_metabolome_new))
  # nodes3_transcript_new))
  # remove duplicates 
  nodes_new_nodup <- nodes_new[!base::duplicated(nodes_new),] %>%
    base::as.data.frame()
  
  
  # 6. Plotly Version ----
  # library(plotly)
  
  # Add color scheme to nodes
  nodes_new_plotly <- nodes_new_nodup %>%
    left_join(color_pal_sankey) %>%
    mutate(
      x = case_when(
        group == "exposure" ~ 0,
        str_detect(group, "CpG") |
          str_detect(group, "Microbiome") ~ 1/5,
        str_detect(group, "Metabolome") ~ 2/5,
        str_detect(group, "outcome") ~ 3/5,
      ))
  
  
  ## 6.2 Get links for Plotly, set color ----
  links_new <- links_new  %>%
    mutate(
      link_color = case_when(
        # Ref link color
        value == 0 ~     "#f3f6f4", 
        # microbiome 
        str_detect(target, "eGFR_log_scaled") &  group == TRUE  ~  "red",
        str_detect(target, "eGFR_log_scaled") &  group == FALSE  ~  "#e4e5f2",
        
        # str_detect(source, "Transcriptome") &  group == TRUE  ~  "#bf9000",
        # str_detect(source, "Transcriptome") &  group == FALSE ~  "#ffd966",
        # Transcriptome
        str_detect(source, "Microbiome") &  group == TRUE  ~  "#38761d",
        str_detect(source, "Microbiome") &  group == FALSE ~  "#b6d7a8",
        # proteome
        str_detect(source, "Metabolome") &  group == TRUE  ~  "#a64d79",
        str_detect(source, "Metabolome") &  group == FALSE ~  "#ead1dc",
        
        links_new$group == FALSE ~ "#d9d2e9", # Negative association
        links_new$group == TRUE ~  "red")) # Positive association
  
  plotly_link <- list(
    source = links_new$IDsource,
    target = links_new$IDtarget,
    value = links_new$value+.00000000000000000000001, 
    color = links_new$link_color)  
  
  
  # Get list of nodes for Plotly
  plotly_node <- list(
    label = nodes_new_plotly$name, 
    color = nodes_new_plotly$color,
    pad = 15,
    thickness = 20,
    line = list(color = "black",width = 0.5), 
    x = nodes_new_plotly$x, 
    y = c(0.01, 
          0.1, 0.3, # microbiome clusters
          .45, .55, # Transcriptome clusters
          .80, .95, # Proteome clusters
          seq(from = .01, to = 1, by = 0.035)[1:n_omics_1], # Cpgs (10 total)
          seq(from = 0.35, to = 1, by = 0.025)[1:n_omics_2], # metabolome (8 total)
          # seq(from = 0.75, to = 1, by = 0.03)[1:n_omics_3], # Transcript (10 total)
          .95
    ))
  
  
  ## 6.3 Plot Figure ----
  fig <- plot_ly(
    type = "sankey",
    domain = list(
      x =  c(0,1),
      y =  c(0,1)),
    orientation = "h",
    node = plotly_node,
    link = plotly_link)
  
  (fig <- fig %>% layout(
    # title = "Basic Sankey Diagram",
    font = list(
      size = text_size
    ))
  )
  
  return(fig)
}
# col_pal <- RColorBrewer::brewer.pal(n = 8, name = "Dark2")
# color_pal_sankey <- matrix(c("exposure", "red",
#                              "lc"      , "#b3d8ff",
#                              "Microbiome","#66A61E",
#                              "CpG"     , "#66A61E",
#                              "Metabolome","#E7298A",
#                              "outcome" , "yellow"),
#                            ncol = 2, byrow = TRUE) %>%
#   as_tibble(.name_repair = "unique") %>%
#   janitor::clean_names() %>%
#   dplyr::rename(group = x1, color = x2)



# p3<- sankey_in_serial(fit1, 
#                       fit2, 
#                       # fit3, 
#                       color_pal_sankey,
#                       text_size = 24)
