#Create name reference dataframe
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

# Create taxonomic structure reference dataframe
create_taxa_structure <- function() {
  
  # Start with species (base level)
  taxa_structure <- data.frame(
    species = rownames(Z.s.g),
    stringsAsFactors = FALSE
  )
  
  # Add genus
  taxa_structure$genus <- sapply(taxa_structure$species, function(s) {
    colnames(Z.s.g)[which(Z.s.g[s, ] == 1)]
  })
  
  # Add family
  taxa_structure$family <- sapply(taxa_structure$genus, function(g) {
    colnames(Z.g.f)[which(Z.g.f[g, ] == 1)]
  })
  
  # Add order
  taxa_structure$order <- sapply(taxa_structure$family, function(f) {
    colnames(Z.f.o)[which(Z.f.o[f, ] == 1)]
  })
  
  # Add class
  taxa_structure$class <- sapply(taxa_structure$order, function(o) {
    colnames(Z.o.c)[which(Z.o.c[o, ] == 1)]
  })
  
  # Add phylum
  taxa_structure$phylum <- sapply(taxa_structure$class, function(c) {
    colnames(Z.c.p)[which(Z.c.p[c, ] == 1)]
  })
  
  # Add full names column
  taxa_structure$full.names <- paste(
    taxa_structure$species,
    taxa_structure$genus,
    taxa_structure$family,
    taxa_structure$order,
    taxa_structure$class,
    taxa_structure$phylum,
    sep = "_"
  )
  
  return(taxa_structure)
}

# Create the dataframe
taxa_structure <- create_taxa_structure()

# Save it
# saveRDS(taxa_structure, fs::path(dir.dat,"taxa_structure.RDS"))
