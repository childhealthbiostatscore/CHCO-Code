suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(here))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(reticulate    ))
suppressPackageStartupMessages(library(tibble    ))
suppressPackageStartupMessages(library(tidyr    ))
suppressPackageStartupMessages(library(caret    ))
suppressPackageStartupMessages(library(pROC    ))
suppressPackageStartupMessages(library(broom))
suppressPackageStartupMessages(library(colorspace))

conda_env_path = "/mmfs1/gscratch/togo/matteo/software/conda/r_deps"
use_condaenv(conda_env_path)

Sys.setenv(
  #CPATH           = paste(paste(conda_env_path, "include", sep = "/"), Sys.getenv("CPATH"          ), sep = ":"),
  #LIBRARY_PATH    = paste(paste(conda_env_path, "lib"    , sep = "/"), Sys.getenv("LIBRARY_PATH"   ), sep = ":"),
  #LD_LIBRARY_PATH = paste(paste(conda_env_path, "lib"    , sep = "/"), Sys.getenv("LD_LIBRARY_PATH"), sep = ":"),
  PATH            = paste(paste(conda_env_path, "bin"    , sep = "/"), Sys.getenv("PATH"           ), sep = ":")
)

#my_lib_path = "/mmfs1/gscratch/togo/matteo/software/R"
#.libPaths(my_lib_path)

# tell the compiler where to find headers/libraries
#Sys.setenv(
#  #PYTHONHOME = "/gscratch/scrubbed/mdanto/envs/r_deps/bin/python", # OLD, maybe useless
#  CPATH           = paste(paste(conda_env_path, "include", sep = "/"), Sys.getenv("CPATH"          ), sep = ":"),
#  LIBRARY_PATH    = paste(paste(conda_env_path, "lib"    , sep = "/"), Sys.getenv("LIBRARY_PATH"   ), sep = ":"),
#  LD_LIBRARY_PATH = paste(paste(conda_env_path, "lib"    , sep = "/"), Sys.getenv("LD_LIBRARY_PATH"), sep = ":"),
#  PATH            = paste(paste(conda_env_path, "bin"    , sep = "/"), Sys.getenv("PATH"           ), sep = ":")
#)

################################################################################
# Colors

analysis2color = data.frame(
  analysis = c("metabolomics_plasma.standard_covariates", "metabolomics_urine.standard_covariates", "proteomics.standard_covariates"), 
  color    = qualitative_hcl(n = 3, palette = "Cold"))
