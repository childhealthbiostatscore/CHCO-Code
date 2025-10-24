######## Brain Biomarkers and Proteomics Analysis 



############ Brain Biomarkers Analysis 

library(scran)
library(future)
library(future.apply)
library(tidyverse)
library(colorspace)
library(patchwork)
library(ggdendro)
library(cowplot)
library(ggpubr)
library(rstatix)
library(arsenal)
library(Biobase)
library(msigdbr)
library(kableExtra)
library(knitr)
library(REDCapR)
library(data.table)
library(emmeans)
library(NMF)
library(pheatmap)
library(UpSetR)
library(enrichR)
library(WriteXLS)
library(SAVER)
library(readxl)
library(limma)
library(edgeR)
library(BiocGenerics)
library(GSEABase)
library(slingshot)
library(SingleCellExperiment)
library(MAST)
library(muscat)
library(scater)
library(Seurat)
library(jsonlite)
library(dplyr)
library(glmmTMB)
library(reshape2)
library(broom.mixed)
library(nebula)
#library(table1)
library(clusterProfiler)
library('org.Hs.eg.db')
library(dplyr)
library(ggplot2)
library(tidyr)
library(gridExtra)
library(readxl)
library(stringr)
library(httr)

qx_var <- c("ab40_avg_conc","ab42_avg_conc","tau_avg_conc",
            "nfl_avg_conc","gfap_avg_conc","ptau_181_avg_conc","ptau_217_avg_conc")







# Function to get protein names from UniProt IDs
get_protein_names <- function(uniprot_ids) {
  
  # Remove NAs and duplicates
  unique_ids <- unique(uniprot_ids[!is.na(uniprot_ids)])
  
  if(length(unique_ids) == 0) {
    return(data.frame(Entry = character(), 
                      Gene.Names = character(), 
                      Protein.names = character()))
  }
  
  # Split into batches of 100 (API limit)
  batch_size <- 100
  batches <- split(unique_ids, ceiling(seq_along(unique_ids) / batch_size))
  
  all_results <- data.frame()
  
  cat("Retrieving protein names from UniProt...\n")
  
  for (i in seq_along(batches)) {
    cat(paste0("Processing batch ", i, " of ", length(batches), "...\n"))
    
    batch <- batches[[i]]
    
    # Query UniProt API
    url <- "https://rest.uniprot.org/uniprotkb/search"
    query <- paste0("accession:(", paste(batch, collapse = " OR "), ")")
    
    response <- GET(
      url,
      query = list(
        query = query,
        format = "tsv",
        fields = "accession,gene_names,protein_name"
      )
    )
    
    if (status_code(response) == 200) {
      # Parse response
      content <- content(response, "text", encoding = "UTF-8")
      if(nchar(content) > 0) {
        batch_results <- read.delim(text = content, sep = "\t", stringsAsFactors = FALSE)
        all_results <- rbind(all_results, batch_results)
      }
    } else {
      warning(paste("Failed to retrieve batch", i, ". Status code:", status_code(response)))
    }
    
    # Be nice to the API
    Sys.sleep(0.5)
  }
  
  cat(paste0("Retrieved information for ", nrow(all_results), " proteins\n"))
  
  return(all_results)
}

# Read harmonized data
harmonized_data <- read.csv("C:/Users/netio/OneDrive - UW/Laura Pyle's files - Biostatistics Core Shared Drive/Data Harmonization/Data Clean/soma_harmonized_dataset.csv", na = '')

dat <- harmonized_data %>% 
  dplyr::select(-dob) %>% 
  arrange(date_of_screen) %>% 
  dplyr::summarise(
    across(where(negate(is.numeric)), ~ ifelse(all(is.na(.x)), NA_character_, last(na.omit(.x)))),
    across(where(is.numeric), ~ ifelse(all(is.na(.x)), NA_real_, mean(na.omit(.x), na.rm=T))),
    .by = c(record_id, visit)
  )

dat <- dat %>% filter(group %in% c('Lean Control', 'Type 2 Diabetes'))

# Read data dictionary
data_dictionary <- readxl::read_xlsx('/Users/netio/Downloads/data_dictionary_master.xlsx')

form_names <- unique(data_dictionary$form_name)
proteo <- form_names[str_which(form_names, pattern = 'proteom')]
metab <- form_names[str_which(form_names, pattern = 'metab')]

variables_class <- c('proteomics', 'az_urine_metabolites', 'metabolomics_aq')

data_dictionary_small <- data_dictionary %>% 
  filter(form_name %in% variables_class)

# Extract UniProt IDs from proteomics variable labels
cat("\n=== Extracting UniProt IDs from proteomics data ===\n")
proteomics_vars <- data_dictionary %>%
  filter(form_name == 'proteomics')

# Extract UniProt IDs (pattern matches standard UniProt accession format)
# This regex matches patterns like P12345, Q9Y123, O95238, etc.
uniprot_ids <- proteomics_vars$label %>%
  str_extract("([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})")

# Get protein names from UniProt
protein_mapping <- get_protein_names(uniprot_ids)

# Create enhanced data dictionary with protein information
data_dictionary_enhanced <- data_dictionary %>%
  mutate(uniprot_id = str_extract(label, "([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})")) %>%
  left_join(
    protein_mapping %>% 
      rename(uniprot_id = Entry,
             gene_name = Gene.Names,
             protein_name = Protein.names),
    by = "uniprot_id"
  ) %>%
  mutate(
    # Create clean gene name (take first gene if multiple)
    gene_name_clean = ifelse(!is.na(gene_name), 
                             str_trim(str_split(gene_name, " ", simplify = TRUE)[,1]),
                             NA),
    # Create enhanced label with protein info
    label_enhanced = case_when(
      !is.na(gene_name_clean) & !is.na(protein_name) ~ paste0(gene_name_clean, " - ", protein_name),
      !is.na(gene_name_clean) ~ gene_name_clean,
      TRUE ~ label
    )
  )

cat(paste0("\nSuccessfully mapped ", sum(!is.na(data_dictionary_enhanced$gene_name)), 
           " proteins out of ", nrow(proteomics_vars), " proteomics variables\n\n"))

# Replace original data dictionary with enhanced version
data_dictionary <- data_dictionary_enhanced

data_dictionary_small <- data_dictionary %>% 
  filter(form_name %in% variables_class)

data_dictionary_small$variable_name <- str_replace_all(data_dictionary_small$variable_name, pattern = '_', replacement = '.')

# Find which columns actually exist
existing_cols <- intersect(data_dictionary_small$variable_name, names(dat))

# Select only those
dat_omics <- dat %>% 
  dplyr::select(record_id, group, sex, age, all_of(qx_var), all_of(existing_cols))

# Check how many were found vs missing
cat("Found:", length(existing_cols), "out of", length(data_dictionary_small$variable_name), "\n")
cat("Missing:", length(data_dictionary_small$variable_name) - length(existing_cols), "\n")





# ============================================================================
# Brain Biomarker vs Proteomics Analysis Pipeline
# ============================================================================

library(tidyverse)
library(broom)

# Set output directory
output_dir <- "/Users/netio/Documents/UofW/Projects/Brain_fMRI_Analysis/omics/"

# Create directory if it doesn't exist
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Define brain biomarkers
biomarkers <- c("ab40_avg_conc", "ab42_avg_conc", "tau_avg_conc",
                "nfl_avg_conc", "gfap_avg_conc", "ptau_181_avg_conc", "ptau_217_avg_conc")

# Get proteomics columns (exclude metadata columns)
proteomics_cols <- setdiff(names(dat_omics), 
                           c("record_id", "group", "sex", "age", biomarkers))

cat(paste0("\nAnalyzing ", length(biomarkers), " biomarkers against ", 
           length(proteomics_cols), " proteins\n"))
cat(paste0("Total comparisons: ", length(biomarkers) * length(proteomics_cols), "\n\n"))

# ============================================================================
# Model 1: biomarker ~ protein + age + sex
# ============================================================================

cat("Running Model 1: biomarker ~ protein + age + sex\n")

results_model1 <- tibble()

for (biomarker in biomarkers) {
  cat(paste0("  Processing ", biomarker, "...\n"))
  
  for (protein in proteomics_cols) {
    # Prepare data - filter for complete cases
    analysis_data <- dat_omics %>%
      select(all_of(c(biomarker, protein, "age", "sex"))) %>%
      drop_na()
    
    # Skip if insufficient data
    if (nrow(analysis_data) < 10) next
    
    # Fit model
    tryCatch({
      formula_str <- paste0(biomarker, " ~ ", protein, " + age + sex")
      model <- lm(as.formula(formula_str), data = analysis_data)
      
      # Extract protein coefficient
      protein_coef <- tidy(model) %>%
        filter(term == protein) %>%
        mutate(
          biomarker = biomarker,
          protein = protein,
          n_obs = nrow(analysis_data),
          model = "Model 1"
        ) %>%
        select(model, biomarker, protein, estimate, std.error, 
               statistic, p.value, n_obs)
      
      results_model1 <- bind_rows(results_model1, protein_coef)
      
    }, error = function(e) {
      warning(paste("Error with", biomarker, "~", protein, ":", e$message))
    })
  }
}

# Add adjusted p-values
results_model1 <- results_model1 %>%
  mutate(
    p.adj.fdr = p.adjust(p.value, method = "fdr"),
    p.adj.bonferroni = p.adjust(p.value, method = "bonferroni"),
    significant_fdr = p.adj.fdr < 0.05,
    significant_bonf = p.adj.bonferroni < 0.05
  ) %>%
  arrange(p.value)

# ============================================================================
# Model 2: biomarker ~ protein * diabetes_status + age + sex
# ============================================================================

cat("\nRunning Model 2: biomarker ~ protein * diabetes_status + age + sex\n")

# Create binary diabetes variable
dat_omics <- dat_omics %>%
  mutate(diabetes = ifelse(group == "Type 2 Diabetes", 1, 0))

results_model2_main <- tibble()
results_model2_interaction <- tibble()

for (biomarker in biomarkers) {
  cat(paste0("  Processing ", biomarker, "...\n"))
  
  for (protein in proteomics_cols) {
    # Prepare data
    analysis_data <- dat_omics %>%
      select(all_of(c(biomarker, protein, "diabetes", "age", "sex"))) %>%
      drop_na()
    
    if (nrow(analysis_data) < 10) next
    
    tryCatch({
      formula_str <- paste0(biomarker, " ~ ", protein, " * diabetes + age + sex")
      model <- lm(as.formula(formula_str), data = analysis_data)
      
      # Extract coefficients
      model_coefs <- tidy(model)
      
      # Main protein effect
      protein_main <- model_coefs %>%
        filter(term == protein) %>%
        mutate(
          biomarker = biomarker,
          protein = protein,
          n_obs = nrow(analysis_data),
          model = "Model 2 - Main Effect"
        ) %>%
        select(model, biomarker, protein, estimate, std.error, 
               statistic, p.value, n_obs)
      
      # Interaction effect
      interaction_term <- paste0(protein, ":diabetes")
      protein_interaction <- model_coefs %>%
        filter(term == interaction_term) %>%
        mutate(
          biomarker = biomarker,
          protein = protein,
          n_obs = nrow(analysis_data),
          model = "Model 2 - Interaction"
        ) %>%
        select(model, biomarker, protein, estimate, std.error, 
               statistic, p.value, n_obs)
      
      results_model2_main <- bind_rows(results_model2_main, protein_main)
      results_model2_interaction <- bind_rows(results_model2_interaction, protein_interaction)
      
    }, error = function(e) {
      warning(paste("Error with", biomarker, "~", protein, ":", e$message))
    })
  }
}

# Add adjusted p-values for Model 2
results_model2_main <- results_model2_main %>%
  mutate(
    p.adj.fdr = p.adjust(p.value, method = "fdr"),
    p.adj.bonferroni = p.adjust(p.value, method = "bonferroni"),
    significant_fdr = p.adj.fdr < 0.05,
    significant_bonf = p.adj.bonferroni < 0.05
  ) %>%
  arrange(p.value)

results_model2_interaction <- results_model2_interaction %>%
  mutate(
    p.adj.fdr = p.adjust(p.value, method = "fdr"),
    p.adj.bonferroni = p.adjust(p.value, method = "bonferroni"),
    significant_fdr = p.adj.fdr < 0.05,
    significant_bonf = p.adj.bonferroni < 0.05
  ) %>%
  arrange(p.value)

# ============================================================================
# Save results to CSV files
# ============================================================================

cat("\nSaving results to CSV files...\n")

write.csv(results_model1, 
          "biomarker_protein_model1_results.csv", 
          row.names = FALSE)

write.csv(results_model2_main, 
          "biomarker_protein_model2_main_effects.csv", 
          row.names = FALSE)

write.csv(results_model2_interaction, 
          "biomarker_protein_model2_interactions.csv", 
          row.names = FALSE)

# Combined results
all_results <- bind_rows(
  results_model1,
  results_model2_main,
  results_model2_interaction
)

write.csv(all_results, 
          "biomarker_protein_all_results.csv", 
          row.names = FALSE)

# Summary statistics
cat("\n=== RESULTS SUMMARY ===\n")
cat("\nModel 1 Results:\n")
cat(paste0("  Total comparisons: ", nrow(results_model1), "\n"))
cat(paste0("  Significant (FDR < 0.05): ", sum(results_model1$significant_fdr), "\n"))
cat(paste0("  Significant (Bonferroni < 0.05): ", sum(results_model1$significant_bonf), "\n"))

cat("\nModel 2 Main Effects:\n")
cat(paste0("  Total comparisons: ", nrow(results_model2_main), "\n"))
cat(paste0("  Significant (FDR < 0.05): ", sum(results_model2_main$significant_fdr), "\n"))
cat(paste0("  Significant (Bonferroni < 0.05): ", sum(results_model2_main$significant_bonf), "\n"))

cat("\nModel 2 Interactions:\n")
cat(paste0("  Total comparisons: ", nrow(results_model2_interaction), "\n"))
cat(paste0("  Significant (FDR < 0.05): ", sum(results_model2_interaction$significant_fdr), "\n"))
cat(paste0("  Significant (Bonferroni < 0.05): ", sum(results_model2_interaction$significant_bonf), "\n"))

# ============================================================================
# Generate plots for significant results
# ============================================================================

cat("\n=== Generating plots for significant results ===\n")

# Create output directory for plots
dir.create("biomarker_protein_plots", showWarnings = FALSE)

# Function to create scatter plot with regression line
create_scatter_plot <- function(biomarker, protein, model_type, p_value, p_adj) {
  
  # Get protein information for title
  protein_info <- data_dictionary %>%
    filter(variable_name == str_replace_all(protein, "\\.", "_")) %>%
    pull(label_enhanced)
  
  protein_label <- if(length(protein_info) > 0) protein_info[1] else protein
  
  plot_data <- dat_omics %>%
    select(all_of(c(biomarker, protein, "group", "age", "sex"))) %>%
    drop_na()
  
  if (nrow(plot_data) < 10) return(NULL)
  
  # Create base plot
  p <- ggplot(plot_data, aes(x = .data[[protein]], y = .data[[biomarker]])) +
    geom_point(aes(color = group), alpha = 0.6, size = 2.5) +
    geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "dashed") +
    labs(
      title = paste0(biomarker, " vs ", protein_label),
      subtitle = paste0(model_type, " | p = ", format.pval(p_value, digits = 3), 
                        " | FDR = ", format.pval(p_adj, digits = 3)),
      x = protein_label,
      y = biomarker,
      color = "Group"
    ) +
    theme_classic() +
    theme(
      plot.title = element_text(size = 12, face = "bold"),
      plot.subtitle = element_text(size = 10),
      legend.position = "bottom"
    )
  
  return(p)
}

# Plot significant Model 1 results (FDR < 0.05)
model1_sig <- results_model1 %>% filter(significant_fdr)

if (nrow(model1_sig) > 0) {
  cat(paste0("\nCreating ", nrow(model1_sig), " plots for Model 1...\n"))
  
  for (i in 1:nrow(model1_sig)) {
    row <- model1_sig[i,]
    
    p <- create_scatter_plot(
      row$biomarker, 
      row$protein, 
      "Model 1",
      row$p.value,
      row$p.adj.fdr
    )
    
    if (!is.null(p)) {
      filename <- paste0("biomarker_protein_plots/model1_", 
                         row$biomarker, "_", row$protein, ".png")
      ggsave(filename, p, width = 8, height = 6, dpi = 300)
    }
  }
}

# Plot significant Model 2 main effects (FDR < 0.05)
model2_main_sig <- results_model2_main %>% filter(significant_fdr)

if (nrow(model2_main_sig) > 0) {
  cat(paste0("\nCreating ", nrow(model2_main_sig), " plots for Model 2 main effects...\n"))
  
  for (i in 1:nrow(model2_main_sig)) {
    row <- model2_main_sig[i,]
    
    # Add diabetes stratification
    plot_data <- dat_omics %>%
      select(all_of(c(row$biomarker, row$protein, "group", "age", "sex"))) %>%
      drop_na()
    
    if (nrow(plot_data) < 10) next
    
    protein_info <- data_dictionary %>%
      filter(variable_name == str_replace_all(row$protein, "\\.", "_")) %>%
      pull(label_enhanced)
    
    protein_label <- if(length(protein_info) > 0) protein_info[1] else row$protein
    
    p <- ggplot(plot_data, aes(x = .data[[row$protein]], y = .data[[row$biomarker]], 
                               color = group)) +
      geom_point(alpha = 0.6, size = 2.5) +
      geom_smooth(method = "lm", se = TRUE) +
      labs(
        title = paste0(row$biomarker, " vs ", protein_label),
        subtitle = paste0("Model 2 Main Effect | p = ", format.pval(row$p.value, digits = 3), 
                          " | FDR = ", format.pval(row$p.adj.fdr, digits = 3)),
        x = protein_label,
        y = row$biomarker,
        color = "Group"
      ) +
      theme_classic() +
      theme(
        plot.title = element_text(size = 12, face = "bold"),
        plot.subtitle = element_text(size = 10),
        legend.position = "bottom"
      )
    
    filename <- paste0("biomarker_protein_plots/model2_main_", 
                       row$biomarker, "_", row$protein, ".png")
    ggsave(filename, p, width = 8, height = 6, dpi = 300)
  }
}

# Plot significant Model 2 interactions (FDR < 0.05)
model2_int_sig <- results_model2_interaction %>% filter(significant_fdr)

if (nrow(model2_int_sig) > 0) {
  cat(paste0("\nCreating ", nrow(model2_int_sig), " plots for Model 2 interactions...\n"))
  
  for (i in 1:nrow(model2_int_sig)) {
    row <- model2_int_sig[i,]
    
    plot_data <- dat_omics %>%
      select(all_of(c(row$biomarker, row$protein, "group", "age", "sex"))) %>%
      drop_na()
    
    if (nrow(plot_data) < 10) next
    
    protein_info <- data_dictionary %>%
      filter(variable_name == str_replace_all(row$protein, "\\.", "_")) %>%
      pull(label_enhanced)
    
    protein_label <- if(length(protein_info) > 0) protein_info[1] else row$protein
    
    p <- ggplot(plot_data, aes(x = .data[[row$protein]], y = .data[[row$biomarker]], 
                               color = group)) +
      geom_point(alpha = 0.6, size = 2.5) +
      geom_smooth(method = "lm", se = TRUE) +
      labs(
        title = paste0(row$biomarker, " vs ", protein_label, " (Interaction)"),
        subtitle = paste0("Model 2 Interaction | p = ", format.pval(row$p.value, digits = 3), 
                          " | FDR = ", format.pval(row$p.adj.fdr, digits = 3)),
        x = protein_label,
        y = row$biomarker,
        color = "Group"
      ) +
      theme_classic() +
      theme(
        plot.title = element_text(size = 12, face = "bold"),
        plot.subtitle = element_text(size = 10),
        legend.position = "bottom"
      )
    
    filename <- paste0("biomarker_protein_plots/model2_interaction_", 
                       row$biomarker, "_", row$protein, ".png")
    ggsave(filename, p, width = 8, height = 6, dpi = 300)
  }
}

cat("\n=== Analysis Complete! ===\n")
cat("\nOutput files created:\n")
cat("  - biomarker_protein_model1_results.csv\n")
cat("  - biomarker_protein_model2_main_effects.csv\n")
cat("  - biomarker_protein_model2_interactions.csv\n")
cat("  - biomarker_protein_all_results.csv\n")
cat("  - biomarker_protein_plots/ (folder with all significant result plots)\n")

























########################################### LC vs. T1D; LC vs. PKD 


library(scran)
library(future)
library(future.apply)
library(tidyverse)
library(colorspace)
library(patchwork)
library(ggdendro)
library(cowplot)
library(ggpubr)
library(rstatix)
library(arsenal)
library(Biobase)
library(msigdbr)
library(kableExtra)
library(knitr)
library(REDCapR)
library(data.table)
library(emmeans)
library(NMF)
library(pheatmap)
library(UpSetR)
library(enrichR)
library(WriteXLS)
library(SAVER)
library(readxl)
library(limma)
library(edgeR)
library(BiocGenerics)
library(GSEABase)
library(slingshot)
library(SingleCellExperiment)
library(MAST)
library(muscat)
library(scater)
library(Seurat)
library(jsonlite)
library(dplyr)
library(glmmTMB)
library(reshape2)
library(broom.mixed)
library(nebula)
#library(table1)
library(clusterProfiler)
library('org.Hs.eg.db')
library(dplyr)
library(ggplot2)
library(tidyr)
library(gridExtra)
library(readxl)
library(stringr)
library(httr)

qx_var <- c("ab40_avg_conc","ab42_avg_conc","tau_avg_conc",
            "nfl_avg_conc","gfap_avg_conc","ptau_181_avg_conc","ptau_217_avg_conc")







# Function to get protein names from UniProt IDs
get_protein_names <- function(uniprot_ids) {
  
  # Remove NAs and duplicates
  unique_ids <- unique(uniprot_ids[!is.na(uniprot_ids)])
  
  if(length(unique_ids) == 0) {
    return(data.frame(Entry = character(), 
                      Gene.Names = character(), 
                      Protein.names = character()))
  }
  
  # Split into batches of 100 (API limit)
  batch_size <- 100
  batches <- split(unique_ids, ceiling(seq_along(unique_ids) / batch_size))
  
  all_results <- data.frame()
  
  cat("Retrieving protein names from UniProt...\n")
  
  for (i in seq_along(batches)) {
    cat(paste0("Processing batch ", i, " of ", length(batches), "...\n"))
    
    batch <- batches[[i]]
    
    # Query UniProt API
    url <- "https://rest.uniprot.org/uniprotkb/search"
    query <- paste0("accession:(", paste(batch, collapse = " OR "), ")")
    
    response <- GET(
      url,
      query = list(
        query = query,
        format = "tsv",
        fields = "accession,gene_names,protein_name"
      )
    )
    
    if (status_code(response) == 200) {
      # Parse response
      content <- content(response, "text", encoding = "UTF-8")
      if(nchar(content) > 0) {
        batch_results <- read.delim(text = content, sep = "\t", stringsAsFactors = FALSE)
        all_results <- rbind(all_results, batch_results)
      }
    } else {
      warning(paste("Failed to retrieve batch", i, ". Status code:", status_code(response)))
    }
    
    # Be nice to the API
    Sys.sleep(0.5)
  }
  
  cat(paste0("Retrieved information for ", nrow(all_results), " proteins\n"))
  
  return(all_results)
}

# Read harmonized data
harmonized_data <- read.csv("C:/Users/netio/OneDrive - UW/Laura Pyle's files - Biostatistics Core Shared Drive/Data Harmonization/Data Clean/soma_harmonized_dataset.csv", na = '')

dat <- harmonized_data %>% 
  dplyr::select(-dob) %>% 
  arrange(date_of_screen) %>% 
  dplyr::summarise(
    across(where(negate(is.numeric)), ~ ifelse(all(is.na(.x)), NA_character_, last(na.omit(.x)))),
    across(where(is.numeric), ~ ifelse(all(is.na(.x)), NA_real_, mean(na.omit(.x), na.rm=T))),
    .by = c(record_id, visit)
  )

#dat <- dat %>% filter(group %in% c('Lean Control', 'Type 2 Diabetes'))

# Read data dictionary
data_dictionary <- readxl::read_xlsx('/Users/netio/Downloads/data_dictionary_master.xlsx')

form_names <- unique(data_dictionary$form_name)
proteo <- form_names[str_which(form_names, pattern = 'proteom')]
metab <- form_names[str_which(form_names, pattern = 'metab')]

variables_class <- c('proteomics', 'az_urine_metabolites', 'metabolomics_aq')

data_dictionary_small <- data_dictionary %>% 
  filter(form_name %in% variables_class)

# Extract UniProt IDs from proteomics variable labels
cat("\n=== Extracting UniProt IDs from proteomics data ===\n")
proteomics_vars <- data_dictionary %>%
  filter(form_name == 'proteomics')

# Extract UniProt IDs (pattern matches standard UniProt accession format)
# This regex matches patterns like P12345, Q9Y123, O95238, etc.
uniprot_ids <- proteomics_vars$label %>%
  str_extract("([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})")

# Get protein names from UniProt
protein_mapping <- get_protein_names(uniprot_ids)

# Create enhanced data dictionary with protein information
data_dictionary_enhanced <- data_dictionary %>%
  mutate(uniprot_id = str_extract(label, "([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})")) %>%
  left_join(
    protein_mapping %>% 
      rename(uniprot_id = Entry,
             gene_name = Gene.Names,
             protein_name = Protein.names),
    by = "uniprot_id"
  ) %>%
  mutate(
    # Create clean gene name (take first gene if multiple)
    gene_name_clean = ifelse(!is.na(gene_name), 
                             str_trim(str_split(gene_name, " ", simplify = TRUE)[,1]),
                             NA),
    # Create enhanced label with protein info
    label_enhanced = case_when(
      !is.na(gene_name_clean) & !is.na(protein_name) ~ paste0(gene_name_clean, " - ", protein_name),
      !is.na(gene_name_clean) ~ gene_name_clean,
      TRUE ~ label
    )
  )

cat(paste0("\nSuccessfully mapped ", sum(!is.na(data_dictionary_enhanced$gene_name)), 
           " proteins out of ", nrow(proteomics_vars), " proteomics variables\n\n"))

# Replace original data dictionary with enhanced version
data_dictionary <- data_dictionary_enhanced

data_dictionary_small <- data_dictionary %>% 
  filter(form_name %in% variables_class)

data_dictionary_small$variable_name <- str_replace_all(data_dictionary_small$variable_name, pattern = '_', replacement = '.')

# Find which columns actually exist
existing_cols <- intersect(data_dictionary_small$variable_name, names(dat))

# Select only those
dat_omics <- dat %>% 
  dplyr::select(record_id, group, sex, age, all_of(qx_var), all_of(existing_cols))

# Check how many were found vs missing
cat("Found:", length(existing_cols), "out of", length(data_dictionary_small$variable_name), "\n")
cat("Missing:", length(data_dictionary_small$variable_name) - length(existing_cols), "\n")



library(gt)
library(gtsummary)


dat_small <- dat %>% filter(group %in% c('Lean Control', 'Type 1 Diabetes', 'Type 2 Diabetes', 'PKD')) %>%
  filter(!is.na(ab40_avg_conc))


desc_table1_fixed <- dat_small %>%
  select(age, sex, group, race_ethnicity, bmi, hba1c, study) %>%
  tbl_summary(
    by = group,
    type = list(
      age ~ "continuous",
      bmi ~ "continuous", 
      sex ~ 'categorical',
      hba1c ~ "continuous",
      race_ethnicity ~ "categorical",
      study ~ "categorical"
      
    ),
    statistic = list(
      all_continuous() ~ "{mean} ({sd})",
      all_categorical() ~ "{n} ({p}%)"
    ),
    digits = list(
      age ~ 1,
      bmi ~ 1,
      hba1c ~ 2,
      all_categorical() ~ c(0, 1)
    ),
    label = list(
      age ~ "Age, years",
      race_ethnicity ~ "Race/Ethnicity",
      sex ~ 'Sex', 
      bmi ~ "BMI, kg/m²",
      hba1c ~ "HbA1c, %",
      study ~ "Study"
    ),
    missing_text = "Missing"
  ) %>%
  add_p(test = list(
    all_continuous() ~ "t.test"
    # Skip categorical p-values if they cause issues
  )) %>%
  add_overall(col_label = "**Overall**\nN = {N}") %>%
  modify_header(label ~ "**Characteristic**") %>%
  modify_spanning_header(all_stat_cols() ~ "**Group**") %>%
  modify_footnote(all_stat_cols() ~ "Mean (SD) for continuous variables; n (%) for categorical variables")

# Save version with epic
desc_table1_fixed %>%
  as_gt() %>%
  tab_options(
    table.font.size = 11,
    heading.title.font.size = 14,
    column_labels.font.size = 12
  ) %>%
  gtsave(paste0("/Users/netio/Documents/UofW/Projects/Brain_fMRI_Analysis/omics/Omics_demographics.png"), 
         vwidth = 1200, vheight = 800)





############## Comparisons for PKD

# ============================================================================
# Brain Biomarker vs Proteomics Analysis: LC vs PKD
# ============================================================================

library(tidyverse)
library(broom)

# Set output directory
output_dir <- "/Users/netio/Documents/UofW/Projects/Brain_fMRI_Analysis/omics/"

# Create directory if it doesn't exist
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Define brain biomarkers
biomarkers <- c("ab40_avg_conc", "ab42_avg_conc", "tau_avg_conc",
                "nfl_avg_conc", "gfap_avg_conc", "ptau_181_avg_conc", "ptau_217_avg_conc")

# Filter data for LC vs PKD comparison
dat_omics_pkd <- dat_omics %>% 
  filter(group %in% c('Lean Control', 'PKD'))

cat("\n=== LC vs PKD Analysis ===\n")
cat(paste0("Lean Control: ", sum(dat_omics_pkd$group == 'Lean Control'), " subjects\n"))
cat(paste0("PKD: ", sum(dat_omics_pkd$group == 'PKD'), " subjects\n"))

# Get proteomics columns
proteomics_cols <- setdiff(names(dat_omics_pkd), 
                           c("record_id", "group", "sex", "age", "diabetes", biomarkers))

cat(paste0("\nAnalyzing ", length(biomarkers), " biomarkers against ", 
           length(proteomics_cols), " proteins\n"))
cat(paste0("Total comparisons: ", length(biomarkers) * length(proteomics_cols), "\n\n"))

# ============================================================================
# Model 1: biomarker ~ protein + age + sex
# ============================================================================

cat("Running Model 1: biomarker ~ protein + age + sex\n")

results_model1 <- tibble()

for (biomarker in biomarkers) {
  cat(paste0("  Processing ", biomarker, "...\n"))
  
  for (protein in proteomics_cols) {
    # Prepare data - filter for complete cases
    analysis_data <- dat_omics_pkd %>%
      select(all_of(c(biomarker, protein, "age", "sex"))) %>%
      drop_na()
    
    # Skip if insufficient data
    if (nrow(analysis_data) < 10) next
    
    # Fit model
    tryCatch({
      formula_str <- paste0(biomarker, " ~ ", protein, " + age + sex")
      model <- lm(as.formula(formula_str), data = analysis_data)
      
      # Extract protein coefficient
      protein_coef <- tidy(model) %>%
        filter(term == protein) %>%
        mutate(
          biomarker = biomarker,
          protein = protein,
          n_obs = nrow(analysis_data),
          model = "Model 1"
        ) %>%
        select(model, biomarker, protein, estimate, std.error, 
               statistic, p.value, n_obs)
      
      results_model1 <- bind_rows(results_model1, protein_coef)
      
    }, error = function(e) {
      warning(paste("Error with", biomarker, "~", protein, ":", e$message))
    })
  }
}

# Add adjusted p-values
results_model1 <- results_model1 %>%
  mutate(
    p.adj.fdr = p.adjust(p.value, method = "fdr"),
    p.adj.bonferroni = p.adjust(p.value, method = "bonferroni"),
    significant_fdr = p.adj.fdr < 0.05,
    significant_bonf = p.adj.bonferroni < 0.05
  ) %>%
  arrange(p.value)

# ============================================================================
# Model 2: biomarker ~ protein * pkd_status + age + sex
# ============================================================================

cat("\nRunning Model 2: biomarker ~ protein * pkd_status + age + sex\n")

# Create binary PKD variable
dat_omics_pkd <- dat_omics_pkd %>%
  mutate(pkd = ifelse(group == "PKD", 1, 0))

results_model2_main <- tibble()
results_model2_interaction <- tibble()

for (biomarker in biomarkers) {
  cat(paste0("  Processing ", biomarker, "...\n"))
  
  for (protein in proteomics_cols) {
    # Prepare data
    analysis_data <- dat_omics_pkd %>%
      select(all_of(c(biomarker, protein, "pkd", "age", "sex"))) %>%
      drop_na()
    
    if (nrow(analysis_data) < 10) next
    
    tryCatch({
      formula_str <- paste0(biomarker, " ~ ", protein, " * pkd + age + sex")
      model <- lm(as.formula(formula_str), data = analysis_data)
      
      # Extract coefficients
      model_coefs <- tidy(model)
      
      # Main protein effect
      protein_main <- model_coefs %>%
        filter(term == protein) %>%
        mutate(
          biomarker = biomarker,
          protein = protein,
          n_obs = nrow(analysis_data),
          model = "Model 2 - Main Effect"
        ) %>%
        select(model, biomarker, protein, estimate, std.error, 
               statistic, p.value, n_obs)
      
      # Interaction effect
      interaction_term <- paste0(protein, ":pkd")
      protein_interaction <- model_coefs %>%
        filter(term == interaction_term) %>%
        mutate(
          biomarker = biomarker,
          protein = protein,
          n_obs = nrow(analysis_data),
          model = "Model 2 - Interaction"
        ) %>%
        select(model, biomarker, protein, estimate, std.error, 
               statistic, p.value, n_obs)
      
      results_model2_main <- bind_rows(results_model2_main, protein_main)
      results_model2_interaction <- bind_rows(results_model2_interaction, protein_interaction)
      
    }, error = function(e) {
      warning(paste("Error with", biomarker, "~", protein, ":", e$message))
    })
  }
}

# Add adjusted p-values for Model 2
results_model2_main <- results_model2_main %>%
  mutate(
    p.adj.fdr = p.adjust(p.value, method = "fdr"),
    p.adj.bonferroni = p.adjust(p.value, method = "bonferroni"),
    significant_fdr = p.adj.fdr < 0.05,
    significant_bonf = p.adj.bonferroni < 0.05
  ) %>%
  arrange(p.value)

results_model2_interaction <- results_model2_interaction %>%
  mutate(
    p.adj.fdr = p.adjust(p.value, method = "fdr"),
    p.adj.bonferroni = p.adjust(p.value, method = "bonferroni"),
    significant_fdr = p.adj.fdr < 0.05,
    significant_bonf = p.adj.bonferroni < 0.05
  ) %>%
  arrange(p.value)

# ============================================================================
# Save results to CSV files
# ============================================================================

cat("\nSaving results to CSV files...\n")

write.csv(results_model1, 
          file.path(output_dir, "LC_vs_PKD_model1_results.csv"), 
          row.names = FALSE)

write.csv(results_model2_main, 
          file.path(output_dir, "LC_vs_PKD_model2_main_effects.csv"), 
          row.names = FALSE)

write.csv(results_model2_interaction, 
          file.path(output_dir, "LC_vs_PKD_model2_interactions.csv"), 
          row.names = FALSE)

# Combined results
all_results <- bind_rows(
  results_model1,
  results_model2_main,
  results_model2_interaction
)

write.csv(all_results, 
          file.path(output_dir, "LC_vs_PKD_all_results.csv"), 
          row.names = FALSE)

# Summary statistics
cat("\n=== RESULTS SUMMARY (LC vs PKD) ===\n")
cat("\nModel 1 Results:\n")
cat(paste0("  Total comparisons: ", nrow(results_model1), "\n"))
cat(paste0("  Significant (FDR < 0.05): ", sum(results_model1$significant_fdr), "\n"))
cat(paste0("  Significant (Bonferroni < 0.05): ", sum(results_model1$significant_bonf), "\n"))

cat("\nModel 2 Main Effects:\n")
cat(paste0("  Total comparisons: ", nrow(results_model2_main), "\n"))
cat(paste0("  Significant (FDR < 0.05): ", sum(results_model2_main$significant_fdr), "\n"))
cat(paste0("  Significant (Bonferroni < 0.05): ", sum(results_model2_main$significant_bonf), "\n"))

cat("\nModel 2 Interactions:\n")
cat(paste0("  Total comparisons: ", nrow(results_model2_interaction), "\n"))
cat(paste0("  Significant (FDR < 0.05): ", sum(results_model2_interaction$significant_fdr), "\n"))
cat(paste0("  Significant (Bonferroni < 0.05): ", sum(results_model2_interaction$significant_bonf), "\n"))

# ============================================================================
# Generate plots for significant results
# ============================================================================

cat("\n=== Generating plots for significant results ===\n")

# Create output directory for plots
plots_dir <- file.path(output_dir, "LC_vs_PKD_plots")
dir.create(plots_dir, showWarnings = FALSE, recursive = TRUE)

# Function to create scatter plot with regression line
create_scatter_plot <- function(biomarker, protein, model_type, p_value, p_adj) {
  
  # Get protein information for title
  protein_info <- data_dictionary %>%
    filter(variable_name == str_replace_all(protein, "\\.", "_")) %>%
    pull(label_enhanced)
  
  protein_label <- if(length(protein_info) > 0) protein_info[1] else protein
  
  plot_data <- dat_omics_pkd %>%
    select(all_of(c(biomarker, protein, "group", "age", "sex"))) %>%
    drop_na()
  
  if (nrow(plot_data) < 10) return(NULL)
  
  # Create base plot
  p <- ggplot(plot_data, aes(x = .data[[protein]], y = .data[[biomarker]])) +
    geom_point(aes(color = group), alpha = 0.6, size = 2.5) +
    geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "dashed") +
    labs(
      title = paste0(biomarker, " vs ", protein_label, " (LC vs PKD)"),
      subtitle = paste0(model_type, " | p = ", format.pval(p_value, digits = 3), 
                        " | FDR = ", format.pval(p_adj, digits = 3)),
      x = protein_label,
      y = biomarker,
      color = "Group"
    ) +
    theme_classic() +
    theme(
      plot.title = element_text(size = 12, face = "bold"),
      plot.subtitle = element_text(size = 10),
      legend.position = "bottom"
    )
  
  return(p)
}

# Plot significant Model 1 results (FDR < 0.05)
model1_sig <- results_model1 %>% filter(significant_fdr)

if (nrow(model1_sig) > 0) {
  cat(paste0("\nCreating ", nrow(model1_sig), " plots for Model 1...\n"))
  
  for (i in 1:nrow(model1_sig)) {
    row <- model1_sig[i,]
    
    p <- create_scatter_plot(
      row$biomarker, 
      row$protein, 
      "Model 1",
      row$p.value,
      row$p.adj.fdr
    )
    
    if (!is.null(p)) {
      filename <- file.path(plots_dir, paste0("model1_", 
                                              row$biomarker, "_", row$protein, ".png"))
      ggsave(filename, p, width = 8, height = 6, dpi = 300)
    }
  }
}

# Plot significant Model 2 main effects (FDR < 0.05)
model2_main_sig <- results_model2_main %>% filter(significant_fdr)

if (nrow(model2_main_sig) > 0) {
  cat(paste0("\nCreating ", nrow(model2_main_sig), " plots for Model 2 main effects...\n"))
  
  for (i in 1:nrow(model2_main_sig)) {
    row <- model2_main_sig[i,]
    
    # Add PKD stratification
    plot_data <- dat_omics_pkd %>%
      select(all_of(c(row$biomarker, row$protein, "group", "age", "sex"))) %>%
      drop_na()
    
    if (nrow(plot_data) < 10) next
    
    protein_info <- data_dictionary %>%
      filter(variable_name == str_replace_all(row$protein, "\\.", "_")) %>%
      pull(label_enhanced)
    
    protein_label <- if(length(protein_info) > 0) protein_info[1] else row$protein
    
    p <- ggplot(plot_data, aes(x = .data[[row$protein]], y = .data[[row$biomarker]], 
                               color = group)) +
      geom_point(alpha = 0.6, size = 2.5) +
      geom_smooth(method = "lm", se = TRUE) +
      labs(
        title = paste0(row$biomarker, " vs ", protein_label, " (LC vs PKD)"),
        subtitle = paste0("Model 2 Main Effect | p = ", format.pval(row$p.value, digits = 3), 
                          " | FDR = ", format.pval(row$p.adj.fdr, digits = 3)),
        x = protein_label,
        y = row$biomarker,
        color = "Group"
      ) +
      theme_classic() +
      theme(
        plot.title = element_text(size = 12, face = "bold"),
        plot.subtitle = element_text(size = 10),
        legend.position = "bottom"
      )
    
    filename <- file.path(plots_dir, paste0("model2_main_", 
                                            row$biomarker, "_", row$protein, ".png"))
    ggsave(filename, p, width = 8, height = 6, dpi = 300)
  }
}

# Plot significant Model 2 interactions (FDR < 0.05)
model2_int_sig <- results_model2_interaction %>% filter(significant_fdr)

if (nrow(model2_int_sig) > 0) {
  cat(paste0("\nCreating ", nrow(model2_int_sig), " plots for Model 2 interactions...\n"))
  
  for (i in 1:nrow(model2_int_sig)) {
    row <- model2_int_sig[i,]
    
    plot_data <- dat_omics_pkd %>%
      select(all_of(c(row$biomarker, row$protein, "group", "age", "sex"))) %>%
      drop_na()
    
    if (nrow(plot_data) < 10) next
    
    protein_info <- data_dictionary %>%
      filter(variable_name == str_replace_all(row$protein, "\\.", "_")) %>%
      pull(label_enhanced)
    
    protein_label <- if(length(protein_info) > 0) protein_info[1] else row$protein
    
    p <- ggplot(plot_data, aes(x = .data[[row$protein]], y = .data[[row$biomarker]], 
                               color = group)) +
      geom_point(alpha = 0.6, size = 2.5) +
      geom_smooth(method = "lm", se = TRUE) +
      labs(
        title = paste0(row$biomarker, " vs ", protein_label, " (LC vs PKD - Interaction)"),
        subtitle = paste0("Model 2 Interaction | p = ", format.pval(row$p.value, digits = 3), 
                          " | FDR = ", format.pval(row$p.adj.fdr, digits = 3)),
        x = protein_label,
        y = row$biomarker,
        color = "Group"
      ) +
      theme_classic() +
      theme(
        plot.title = element_text(size = 12, face = "bold"),
        plot.subtitle = element_text(size = 10),
        legend.position = "bottom"
      )
    
    filename <- file.path(plots_dir, paste0("model2_interaction_", 
                                            row$biomarker, "_", row$protein, ".png"))
    ggsave(filename, p, width = 8, height = 6, dpi = 300)
  }
}

cat("\n=== LC vs PKD Analysis Complete! ===\n")
cat("\nOutput files created in:", output_dir, "\n")
cat("  - LC_vs_PKD_model1_results.csv\n")
cat("  - LC_vs_PKD_model2_main_effects.csv\n")
cat("  - LC_vs_PKD_model2_interactions.csv\n")
cat("  - LC_vs_PKD_all_results.csv\n")
cat("  - LC_vs_PKD_plots/ (folder with all significant result plots)\n")



















########### Comparisons for T1D 



# ============================================================================
# Brain Biomarker vs Proteomics Analysis: LC vs T1D
# ============================================================================

library(tidyverse)
library(broom)

# Set output directory
output_dir <- "/Users/netio/Documents/UofW/Projects/Brain_fMRI_Analysis/omics/"

# Create directory if it doesn't exist
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Define brain biomarkers
biomarkers <- c("ab40_avg_conc", "ab42_avg_conc", "tau_avg_conc",
                "nfl_avg_conc", "gfap_avg_conc", "ptau_181_avg_conc", "ptau_217_avg_conc")

# Filter data for LC vs T1D comparison
dat_omics_t1d <- dat_omics %>% 
  filter(group %in% c('Lean Control', 'Type 1 Diabetes'))

cat("\n=== LC vs T1D Analysis ===\n")
cat(paste0("Lean Control: ", sum(dat_omics_t1d$group == 'Lean Control'), " subjects\n"))
cat(paste0("Type 1 Diabetes: ", sum(dat_omics_t1d$group == 'Type 1 Diabetes'), " subjects\n"))

# Get proteomics columns
proteomics_cols <- setdiff(names(dat_omics_t1d), 
                           c("record_id", "group", "sex", "age", "diabetes", biomarkers))

cat(paste0("\nAnalyzing ", length(biomarkers), " biomarkers against ", 
           length(proteomics_cols), " proteins\n"))
cat(paste0("Total comparisons: ", length(biomarkers) * length(proteomics_cols), "\n\n"))

# ============================================================================
# Model 1: biomarker ~ protein + age + sex
# ============================================================================

cat("Running Model 1: biomarker ~ protein + age + sex\n")

results_model1 <- tibble()

for (biomarker in biomarkers) {
  cat(paste0("  Processing ", biomarker, "...\n"))
  
  for (protein in proteomics_cols) {
    # Prepare data - filter for complete cases
    analysis_data <- dat_omics_t1d %>%
      select(all_of(c(biomarker, protein, "age", "sex"))) %>%
      drop_na()
    
    # Skip if insufficient data
    if (nrow(analysis_data) < 10) next
    
    # Fit model
    tryCatch({
      formula_str <- paste0(biomarker, " ~ ", protein, " + age + sex")
      model <- lm(as.formula(formula_str), data = analysis_data)
      
      # Extract protein coefficient
      protein_coef <- tidy(model) %>%
        filter(term == protein) %>%
        mutate(
          biomarker = biomarker,
          protein = protein,
          n_obs = nrow(analysis_data),
          model = "Model 1"
        ) %>%
        select(model, biomarker, protein, estimate, std.error, 
               statistic, p.value, n_obs)
      
      results_model1 <- bind_rows(results_model1, protein_coef)
      
    }, error = function(e) {
      warning(paste("Error with", biomarker, "~", protein, ":", e$message))
    })
  }
}

# Add adjusted p-values
results_model1 <- results_model1 %>%
  mutate(
    p.adj.fdr = p.adjust(p.value, method = "fdr"),
    p.adj.bonferroni = p.adjust(p.value, method = "bonferroni"),
    significant_fdr = p.adj.fdr < 0.05,
    significant_bonf = p.adj.bonferroni < 0.05
  ) %>%
  arrange(p.value)

# ============================================================================
# Model 2: biomarker ~ protein * t1d_status + age + sex
# ============================================================================

cat("\nRunning Model 2: biomarker ~ protein * t1d_status + age + sex\n")

# Create binary T1D variable
dat_omics_t1d <- dat_omics_t1d %>%
  mutate(t1d = ifelse(group == "Type 1 Diabetes", 1, 0))

results_model2_main <- tibble()
results_model2_interaction <- tibble()

for (biomarker in biomarkers) {
  cat(paste0("  Processing ", biomarker, "...\n"))
  
  for (protein in proteomics_cols) {
    # Prepare data
    analysis_data <- dat_omics_t1d %>%
      select(all_of(c(biomarker, protein, "t1d", "age", "sex"))) %>%
      drop_na()
    
    if (nrow(analysis_data) < 10) next
    
    tryCatch({
      formula_str <- paste0(biomarker, " ~ ", protein, " * t1d + age + sex")
      model <- lm(as.formula(formula_str), data = analysis_data)
      
      # Extract coefficients
      model_coefs <- tidy(model)
      
      # Main protein effect
      protein_main <- model_coefs %>%
        filter(term == protein) %>%
        mutate(
          biomarker = biomarker,
          protein = protein,
          n_obs = nrow(analysis_data),
          model = "Model 2 - Main Effect"
        ) %>%
        select(model, biomarker, protein, estimate, std.error, 
               statistic, p.value, n_obs)
      
      # Interaction effect
      interaction_term <- paste0(protein, ":t1d")
      protein_interaction <- model_coefs %>%
        filter(term == interaction_term) %>%
        mutate(
          biomarker = biomarker,
          protein = protein,
          n_obs = nrow(analysis_data),
          model = "Model 2 - Interaction"
        ) %>%
        select(model, biomarker, protein, estimate, std.error, 
               statistic, p.value, n_obs)
      
      results_model2_main <- bind_rows(results_model2_main, protein_main)
      results_model2_interaction <- bind_rows(results_model2_interaction, protein_interaction)
      
    }, error = function(e) {
      warning(paste("Error with", biomarker, "~", protein, ":", e$message))
    })
  }
}

# Add adjusted p-values for Model 2
results_model2_main <- results_model2_main %>%
  mutate(
    p.adj.fdr = p.adjust(p.value, method = "fdr"),
    p.adj.bonferroni = p.adjust(p.value, method = "bonferroni"),
    significant_fdr = p.adj.fdr < 0.05,
    significant_bonf = p.adj.bonferroni < 0.05
  ) %>%
  arrange(p.value)

results_model2_interaction <- results_model2_interaction %>%
  mutate(
    p.adj.fdr = p.adjust(p.value, method = "fdr"),
    p.adj.bonferroni = p.adjust(p.value, method = "bonferroni"),
    significant_fdr = p.adj.fdr < 0.05,
    significant_bonf = p.adj.bonferroni < 0.05
  ) %>%
  arrange(p.value)

# ============================================================================
# Save results to CSV files
# ============================================================================

cat("\nSaving results to CSV files...\n")

write.csv(results_model1, 
          file.path(output_dir, "LC_vs_T1D_model1_results.csv"), 
          row.names = FALSE)

write.csv(results_model2_main, 
          file.path(output_dir, "LC_vs_T1D_model2_main_effects.csv"), 
          row.names = FALSE)

write.csv(results_model2_interaction, 
          file.path(output_dir, "LC_vs_T1D_model2_interactions.csv"), 
          row.names = FALSE)

# Combined results
all_results <- bind_rows(
  results_model1,
  results_model2_main,
  results_model2_interaction
)

write.csv(all_results, 
          file.path(output_dir, "LC_vs_T1D_all_results.csv"), 
          row.names = FALSE)

# Summary statistics
cat("\n=== RESULTS SUMMARY (LC vs T1D) ===\n")
cat("\nModel 1 Results:\n")
cat(paste0("  Total comparisons: ", nrow(results_model1), "\n"))
cat(paste0("  Significant (FDR < 0.05): ", sum(results_model1$significant_fdr), "\n"))
cat(paste0("  Significant (Bonferroni < 0.05): ", sum(results_model1$significant_bonf), "\n"))

cat("\nModel 2 Main Effects:\n")
cat(paste0("  Total comparisons: ", nrow(results_model2_main), "\n"))
cat(paste0("  Significant (FDR < 0.05): ", sum(results_model2_main$significant_fdr), "\n"))
cat(paste0("  Significant (Bonferroni < 0.05): ", sum(results_model2_main$significant_bonf), "\n"))

cat("\nModel 2 Interactions:\n")
cat(paste0("  Total comparisons: ", nrow(results_model2_interaction), "\n"))
cat(paste0("  Significant (FDR < 0.05): ", sum(results_model2_interaction$significant_fdr), "\n"))
cat(paste0("  Significant (Bonferroni < 0.05): ", sum(results_model2_interaction$significant_bonf), "\n"))

# ============================================================================
# Generate plots for significant results
# ============================================================================

cat("\n=== Generating plots for significant results ===\n")

# Create output directory for plots
plots_dir <- file.path(output_dir, "LC_vs_T1D_plots")
dir.create(plots_dir, showWarnings = FALSE, recursive = TRUE)

# Function to create scatter plot with regression line
create_scatter_plot <- function(biomarker, protein, model_type, p_value, p_adj) {
  
  # Get protein information for title
  protein_info <- data_dictionary %>%
    filter(variable_name == str_replace_all(protein, "\\.", "_")) %>%
    pull(label_enhanced)
  
  protein_label <- if(length(protein_info) > 0) protein_info[1] else protein
  
  plot_data <- dat_omics_t1d %>%
    select(all_of(c(biomarker, protein, "group", "age", "sex"))) %>%
    drop_na()
  
  if (nrow(plot_data) < 10) return(NULL)
  
  # Create base plot
  p <- ggplot(plot_data, aes(x = .data[[protein]], y = .data[[biomarker]])) +
    geom_point(aes(color = group), alpha = 0.6, size = 2.5) +
    geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "dashed") +
    labs(
      title = paste0(biomarker, " vs ", protein_label, " (LC vs T1D)"),
      subtitle = paste0(model_type, " | p = ", format.pval(p_value, digits = 3), 
                        " | FDR = ", format.pval(p_adj, digits = 3)),
      x = protein_label,
      y = biomarker,
      color = "Group"
    ) +
    theme_classic() +
    theme(
      plot.title = element_text(size = 12, face = "bold"),
      plot.subtitle = element_text(size = 10),
      legend.position = "bottom"
    )
  
  return(p)
}

# Plot significant Model 1 results (FDR < 0.05)
model1_sig <- results_model1 %>% filter(significant_fdr)

if (nrow(model1_sig) > 0) {
  cat(paste0("\nCreating ", nrow(model1_sig), " plots for Model 1...\n"))
  
  for (i in 1:nrow(model1_sig)) {
    row <- model1_sig[i,]
    
    p <- create_scatter_plot(
      row$biomarker, 
      row$protein, 
      "Model 1",
      row$p.value,
      row$p.adj.fdr
    )
    
    if (!is.null(p)) {
      filename <- file.path(plots_dir, paste0("model1_", 
                                              row$biomarker, "_", row$protein, ".png"))
      ggsave(filename, p, width = 8, height = 6, dpi = 300)
    }
  }
}

# Plot significant Model 2 main effects (FDR < 0.05)
model2_main_sig <- results_model2_main %>% filter(significant_fdr)

if (nrow(model2_main_sig) > 0) {
  cat(paste0("\nCreating ", nrow(model2_main_sig), " plots for Model 2 main effects...\n"))
  
  for (i in 1:nrow(model2_main_sig)) {
    row <- model2_main_sig[i,]
    
    # Add T1D stratification
    plot_data <- dat_omics_t1d %>%
      select(all_of(c(row$biomarker, row$protein, "group", "age", "sex"))) %>%
      drop_na()
    
    if (nrow(plot_data) < 10) next
    
    protein_info <- data_dictionary %>%
      filter(variable_name == str_replace_all(row$protein, "\\.", "_")) %>%
      pull(label_enhanced)
    
    protein_label <- if(length(protein_info) > 0) protein_info[1] else row$protein
    
    p <- ggplot(plot_data, aes(x = .data[[row$protein]], y = .data[[row$biomarker]], 
                               color = group)) +
      geom_point(alpha = 0.6, size = 2.5) +
      geom_smooth(method = "lm", se = TRUE) +
      labs(
        title = paste0(row$biomarker, " vs ", protein_label, " (LC vs T1D)"),
        subtitle = paste0("Model 2 Main Effect | p = ", format.pval(row$p.value, digits = 3), 
                          " | FDR = ", format.pval(row$p.adj.fdr, digits = 3)),
        x = protein_label,
        y = row$biomarker,
        color = "Group"
      ) +
      theme_classic() +
      theme(
        plot.title = element_text(size = 12, face = "bold"),
        plot.subtitle = element_text(size = 10),
        legend.position = "bottom"
      )
    
    filename <- file.path(plots_dir, paste0("model2_main_", 
                                            row$biomarker, "_", row$protein, ".png"))
    ggsave(filename, p, width = 8, height = 6, dpi = 300)
  }
}

# Plot significant Model 2 interactions (FDR < 0.05)
model2_int_sig <- results_model2_interaction %>% filter(significant_fdr)

if (nrow(model2_int_sig) > 0) {
  cat(paste0("\nCreating ", nrow(model2_int_sig), " plots for Model 2 interactions...\n"))
  
  for (i in 1:nrow(model2_int_sig)) {
    row <- model2_int_sig[i,]
    
    plot_data <- dat_omics_t1d %>%
      select(all_of(c(row$biomarker, row$protein, "group", "age", "sex"))) %>%
      drop_na()
    
    if (nrow(plot_data) < 10) next
    
    protein_info <- data_dictionary %>%
      filter(variable_name == str_replace_all(row$protein, "\\.", "_")) %>%
      pull(label_enhanced)
    
    protein_label <- if(length(protein_info) > 0) protein_info[1] else row$protein
    
    p <- ggplot(plot_data, aes(x = .data[[row$protein]], y = .data[[row$biomarker]], 
                               color = group)) +
      geom_point(alpha = 0.6, size = 2.5) +
      geom_smooth(method = "lm", se = TRUE) +
      labs(
        title = paste0(row$biomarker, " vs ", protein_label, " (LC vs T1D - Interaction)"),
        subtitle = paste0("Model 2 Interaction | p = ", format.pval(row$p.value, digits = 3), 
                          " | FDR = ", format.pval(row$p.adj.fdr, digits = 3)),
        x = protein_label,
        y = row$biomarker,
        color = "Group"
      ) +
      theme_classic() +
      theme(
        plot.title = element_text(size = 12, face = "bold"),
        plot.subtitle = element_text(size = 10),
        legend.position = "bottom"
      )
    
    filename <- file.path(plots_dir, paste0("model2_interaction_", 
                                            row$biomarker, "_", row$protein, ".png"))
    ggsave(filename, p, width = 8, height = 6, dpi = 300)
  }
}

cat("\n=== LC vs T1D Analysis Complete! ===\n")
cat("\nOutput files created in:", output_dir, "\n")
cat("  - LC_vs_T1D_model1_results.csv\n")
cat("  - LC_vs_T1D_model2_main_effects.csv\n")
cat("  - LC_vs_T1D_model2_interactions.csv\n")
cat("  - LC_vs_T1D_all_results.csv\n")
cat("  - LC_vs_T1D_plots/ (folder with all significant result plots)\n")


























