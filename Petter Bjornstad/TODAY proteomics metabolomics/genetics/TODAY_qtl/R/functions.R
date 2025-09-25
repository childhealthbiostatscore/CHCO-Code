history(0)

suppressPackageStartupMessages(library(reticulate    ))

my_lib_path = "/mmfs1/gscratch/togo/matteo/software/R"
.libPaths(my_lib_path)

conda_env_path = "/mmfs1/gscratch/togo/matteo/software/conda/r_deps"
# tell the compiler where to find headers/libraries
Sys.setenv(
  #PYTHONHOME = "/gscratch/scrubbed/mdanto/envs/r_deps/bin/python", # OLD, maybe useless
  CPATH           = paste(paste(conda_env_path, "include", sep = "/"), Sys.getenv("CPATH"          ), sep = ":"),
  LIBRARY_PATH    = paste(paste(conda_env_path, "lib"    , sep = "/"), Sys.getenv("LIBRARY_PATH"   ), sep = ":"),
  LD_LIBRARY_PATH = paste(paste(conda_env_path, "lib"    , sep = "/"), Sys.getenv("LD_LIBRARY_PATH"), sep = ":"),
  PATH            = paste(paste(conda_env_path, "bin"    , sep = "/"), Sys.getenv("PATH"           ), sep = ":")
)

#Sys.setenv(
#  PYTHONHOME = "/mmfs1/gscratch/scrubbed/mdanto/envs/bioinfo",
#  PATH = paste("/mmfs1/gscratch/scrubbed/mdanto/envs/bioinfo/bin", Sys.getenv("PATH"), sep = ":")
#)


suppressPackageStartupMessages(library(data.table    ))
suppressPackageStartupMessages(library(dplyr         ))
suppressPackageStartupMessages(library(tidyr         ))
suppressPackageStartupMessages(library(stringr       ))
suppressPackageStartupMessages(library(here          ))
suppressPackageStartupMessages(library(tibble        ))
suppressPackageStartupMessages(library(openxlsx      ))
suppressPackageStartupMessages(library(RNOmni        ))
suppressPackageStartupMessages(library(umap          ))
suppressPackageStartupMessages(library(purrr         ))
suppressPackageStartupMessages(library(qqman         ))
