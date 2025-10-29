library(aws.s3)
library(jsonlite)
library(biomaRt)

user <- Sys.info()[["user"]]

if (user == "choiyej") { # local version
  root_path <- "/Users/choiyej/Library/CloudStorage/OneDrive-SharedLibraries-UW/Laura Pyle - Bjornstad/Biostatistics Core Shared Drive"
  git_path <- "/Users/choiyej/GitHub/CHCO-Code/Petter Bjornstad"
  keys <- fromJSON("/Users/choiyej/Library/CloudStorage/OneDrive-TheUniversityofColoradoDenver/Bjornstad Pyle Lab/keys.json")
} else if (user == "rameshsh") { # hyak version
  root_path <- ""
  git_path <- "/mmfs1/gscratch/togo/rameshsh/CHCO-Code/Petter Bjornstad"
} else if (user == "yejichoi") { # hyak version
  root_path <- "/mmfs1/gscratch/togo/yejichoi/"
  git_path <- "/mmfs1/gscratch/togo/yejichoi/CHCO-Code/Petter Bjornstad"
  keys <- fromJSON("/mmfs1/home/yejichoi/keys.json")
} else if (user == "pylell") {
  root_path <- "/Users/pylell/Library/CloudStorage/OneDrive-SharedLibraries-UW/Bjornstad/Biostatistics Core Shared Drive"
  git_path <- "/Users/pylell/Documents/GitHub/CHCO-Code/Petter Bjornstad"
  keys <- fromJSON("/mmfs1/home/pylell/keys.json")
} else {
  stop("Unknown user: please specify root path for this user.")
}

source(file.path(git_path, "Renal HEIRitage/RH_RH2_IMPROVE_functions.R"))

# pull in nebula results
folders <- c("DKD_vs_nonDKD_100" = "dkd100", 
             "DKD_100_vs_hc" = "dkd100_hc", 
             "nonDKD_100_vs_hc" = "nondkd100_hc", 
             "DKD_vs_nonDKD_30" = "dkd30", 
             "DKD_30_vs_hc" = "dkd30_hc", 
             "nonDKD_30_vs_hc" = "nondkd30_hc", 
             "GLP_N_vs_HC" = "glpn_hc", 
             "GLP_Y_vs_GLP_N" = "glpy_glpn")

for (folder in names(folders)) {
  for (cell in names(celltype_groups)){
    s3readRDS(object = paste0("Renal HERITAGE/Results/nebula/", folder, "/", cell, "_rh_rh2_imp_nebula_kpmp_", folders[folder], "_processed.rds"),
              bucket = "scrna", region = "")
  }
}
