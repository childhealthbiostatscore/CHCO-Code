# ============================================================
# STEP 1: BUILD phylogenic_tree_names_df FROM Z-MATRICES ONLY
# ============================================================

phylogenic_tree_names_df <- data.frame(species_names_full = rownames(Z.s.g)) %>%
  mutate(
    genus_name_full  = sapply(species_names_full, function(s) {
      g <- which(Z.s.g[s, ] == 1)
      if (length(g) == 0) NA else colnames(Z.s.g)[g]
    }),
    family_name_full = sapply(genus_name_full, function(g) {
      if (is.na(g)) return(NA)
      f <- which(Z.g.f[g, ] == 1)
      if (length(f) == 0) NA else colnames(Z.g.f)[f]
    }),
    order_name_full  = sapply(family_name_full, function(f) {
      if (is.na(f)) return(NA)
      o <- which(Z.f.o[f, ] == 1)
      if (length(o) == 0) NA else colnames(Z.f.o)[o]
    }),
    class_name_full  = sapply(order_name_full, function(o) {
      if (is.na(o)) return(NA)
      cl <- which(Z.o.c[o, ] == 1)
      if (length(cl) == 0) NA else colnames(Z.o.c)[cl]
    }),
    phylum_name_full = sapply(class_name_full, function(cl) {
      if (is.na(cl)) return(NA)
      p <- which(Z.c.p[cl, ] == 1)
      if (length(p) == 0) NA else colnames(Z.c.p)[p]
    })
  )

# ============================================================
# STEP 2: BUILD NEWICK TREE
# ============================================================

build_tree <- function(tree_names_df) {
  species_names <- tree_names_df$species_names_full
  genus_names   <- unique(tree_names_df$genus_name_full)
  family_names  <- unique(tree_names_df$family_name_full)
  order_names   <- unique(tree_names_df$order_name_full)
  class_names   <- unique(tree_names_df$class_name_full)
  phylum_names  <- unique(tree_names_df$phylum_name_full)
  
  l1 <- list()
  for (g in genus_names) {
    children <- species_names[tree_names_df$genus_name_full == g]
    l1[[g]]  <- paste0("(", paste(children, collapse = ","), ")", g)
  }
  
  l2 <- list()
  for (f in family_names) {
    children <- unique(tree_names_df$genus_name_full[tree_names_df$family_name_full == f])
    l2[[f]]  <- paste0("(", paste(unlist(l1[children]), collapse = ","), ")", f)
  }
  
  l3 <- list()
  for (o in order_names) {
    children <- unique(tree_names_df$family_name_full[tree_names_df$order_name_full == o])
    l3[[o]]  <- paste0("(", paste(unlist(l2[children]), collapse = ","), ")", o)
  }
  
  l4 <- list()
  for (cl in class_names) {
    children <- unique(tree_names_df$order_name_full[tree_names_df$class_name_full == cl])
    l4[[cl]] <- paste0("(", paste(unlist(l3[children]), collapse = ","), ")", cl)
  }
  
  l5 <- list()
  for (p in phylum_names) {
    children <- unique(tree_names_df$class_name_full[tree_names_df$phylum_name_full == p])
    l5[[p]]  <- paste0("(", paste(unlist(l4[children]), collapse = ","), ")", p)
  }
  
  newick <- paste0("(", paste(unlist(l5), collapse = ","), ")root;")
  ape::read.tree(text = newick)
}

# ============================================================
# STEP 3: EXTRACT EXPECTED VALUES FROM FORMATTED DATA
# expectedLogOdds is already correctly calculated in format.fxn
# Just pull unique value per taxa_full for a given exposure type
# ============================================================

get_expected_data <- function(formatted.data, exposure_type = "individual") {
  
  # Build complete taxa list from Z-matrices at every level
  all_taxa_df <- data.frame(label = c(
    rownames(Z.s.g),  # all species
    colnames(Z.s.g),  # all genera
    colnames(Z.g.f),  # all families
    colnames(Z.f.o),  # all orders
    colnames(Z.o.c),  # all classes
    colnames(Z.c.p)   # all phyla
  ))
  
  # Get expected values from formatted data
  if (exposure_type == "individual") {
    exp_data <- formatted.data %>%
      filter(exposure != "Mixture") %>%
      dplyr::select(taxa_full, expectedLogOdds) %>%
      group_by(taxa_full) %>%
      summarize(ExpVal = mean(expectedLogOdds), .groups = "drop") %>%
      rename(label = taxa_full)
    
  } else {
    exp_data <- formatted.data %>%
      filter(exposure == "Mixture") %>%
      dplyr::select(taxa_full, expectedLogOdds) %>%
      group_by(taxa_full) %>%
      summarize(ExpVal = mean(expectedLogOdds), .groups = "drop") %>%
      rename(label = taxa_full)
  }
  
  # Join - all taxa not in formatted.data get ExpVal = 0 (gray dot)
  all_taxa_df %>%
    left_join(exp_data, by = "label") %>%
    mutate(ExpVal = replace_na(ExpVal, 0))
}

build_dendrograms <- function(formatted.data) {
  
  scen_num   <- unique(formatted.data$Scenario)
  sim.par    <- sim.par.all[scen_num, ]
  setOR      <- sim.par$OR.exposure
  P.e.causal <- sim.par$P.e.causal
  scen_label <- paste0("Scenario ", scen_num)
  
  par      <- sim.par$P.s.scenario
  causal.s <- if (par == 1) NULL else paste0("species", simulation.parameters[[par]])
  
  # Always use FULL tree with all species/taxa from Z-matrices
  tree_i <- build_tree(phylogenic_tree_names_df)
  
  exp_indiv <- get_expected_data(formatted.data, exposure_type = "individual")
  exp_mix   <- get_expected_data(formatted.data, exposure_type = "mixture")
  
  # For null scenarios set tiny max so scale doesn't break
  max_val_indiv <- if (is.null(causal.s) || P.e.causal == 0) 0.001 else max(exp_indiv$ExpVal, na.rm = TRUE)
  max_val_mix   <- if (is.null(causal.s) || P.e.causal == 0) 0.001 else max(exp_mix$ExpVal,   na.rm = TRUE)
  
  list(
    individual = plot_dendrogram(
      tree          = tree_i,
      expected_data = exp_indiv,
      min_val       = 0,
      max_val       = max_val_indiv,
      label         = paste0(scen_label, " - Individual")
    ),
    mixture = plot_dendrogram(
      tree          = tree_i,
      expected_data = exp_mix,
      min_val       = 0,
      max_val       = max_val_mix,
      label         = paste0(scen_label, " - Mixture")
    )
  )
}

plot_dendrogram <- function(tree, expected_data, min_val, max_val,
                            open_angle = 225, label = NULL) {
  p <- ggtree(tree, layout = "fan", open.angle = open_angle) %<+% expected_data +
    geom_point(
      aes(
        color = ifelse(ExpVal > 0, "Associated", "Not Associated"),
        size  = ifelse(ExpVal > 0, ExpVal, 0.001)
      ),
      alpha = 0.9
    ) +
    scale_color_manual(
      values   = c("Associated" = "firebrick4", "Not Associated" = "gray60"),
      na.value = "gray60"
    ) +
    scale_size_continuous(
      range  = c(1, 8),
      limits = c(0, max(max_val, 0.001))
    ) +
    theme(legend.position = "none")
  
  if (!is.null(label)) {
    p <- p + labs(title = label) +
      theme(plot.title = element_text(family = "Times", size = 12, hjust = 0.5))
  }
  return(p)
}

# Rebuild all dendrograms
all_dendros <- all_formatted_results %>%
  group_split(Scenario) %>%
  lapply(build_dendrograms)

# View
wrap_plots(lapply(all_dendros, function(x) x$individual), ncol = 3)
wrap_plots(lapply(all_dendros, function(x) x$mixture),    ncol = 3)

# ============================================================
# STEP 5: BUILD ALL DENDROGRAMS FROM FORMATTED DATA
# Call this after running format.fxn on each scenario
# ============================================================
# formatted.data <- all_formatted_results %>% 
#   filter(Scenario==3)
build_dendrograms <- function(formatted.data) {
  
  scen_num   <- unique(formatted.data$Scenario)
  sim.par    <- sim.par.all[scen_num, ]
  setOR      <- sim.par$OR.exposure
  P.e.causal <- sim.par$P.e.causal
  scen_label <- paste0("Scenario ", scen_num)
  
  # Get causal species for this scenario to filter tree
  par      <- sim.par$P.s.causal.all
  causal.s <- if (par == 1) NULL else paste0("species", simulation.parameters[[par]])
  
  # Filter tree to causal subtree
  if (is.null(causal.s)) {
    causal_tree_names_df <- phylogenic_tree_names_df
  } else {
    causal_tree_names_df <- phylogenic_tree_names_df %>%
      filter(species_names_full %in% causal.s)
  }
  
  tree_i <- build_tree(causal_tree_names_df)
  
  # Extract expected values directly from formatted.data
  exp_indiv <- get_expected_data(formatted.data, exposure_type = "individual")
  exp_mix   <- get_expected_data(formatted.data, exposure_type = "mixture")
  
  # Max values for scale
  max_val_indiv <- log(setOR)
  max_val_mix   <- P.e.causal * log(setOR)
  
  # Null scenario - no causal taxa
  if (is.null(causal.s) || P.e.causal == 0) {
    max_val_indiv <- 0.001
    max_val_mix   <- 0.001
  }
  
  list(
    individual = plot_dendrogram(
      tree          = tree_i,
      expected_data = exp_indiv,
      min_val       = 0,
      max_val       = max_val_indiv,
      label         = paste0(scen_label, " - Individual")
    ),
    mixture = plot_dendrogram(
      tree          = tree_i,
      expected_data = exp_mix,
      min_val       = 0,
      max_val       = max_val_mix,
      label         = paste0(scen_label, " - Mixture")
    )
  )
}

# ============================================================
# USAGE: after running format.fxn on each scenario
# ============================================================

all_dendros <- all_formatted_results %>%
  group_split(Scenario) %>%
  lapply(build_dendrograms)

# Access results
all_dendros[[3]]$individual
all_dendros[[3]]$mixture

# View all in a grid
library(patchwork)
wrap_plots(lapply(all_dendros, function(x) x$mixture),    ncol = 3)
wrap_plots(lapply(all_dendros, function(x) x$individual), ncol = 3)