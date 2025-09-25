source("R/functions.R")
source("R/idat_to_vcf_wrapper.r")

idat_to_vcf(here("data"), here("output/variants.vcf.gz"))
