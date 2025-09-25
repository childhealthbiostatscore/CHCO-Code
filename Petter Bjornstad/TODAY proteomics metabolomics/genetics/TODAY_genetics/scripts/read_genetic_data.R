source("R/connect_to_bucket.R")

## Read file in the TODAY TAR files

temp_file1 = tempfile(fileext = ".tar")
temp_file2 = "data/TODAY_WHOLE_EXOMES.tar"

#s3$download_file(bucket, "TODAY_GWAS_DATA.tar"   , temp_file1)
#s3$download_file(bucket, "TODAY_WHOLE_EXOMES.tar", temp_file2)

#writeLines(untar(temp_file1, list = TRUE), con = "data/TODAY_GWAS_DATA.list_files.txt"   , sep = "\n")
#writeLines(untar(temp_file2, list = TRUE), con = "data/TODAY_WHOLE_EXOMES.list_files.txt", sep = "\n")


# Unpack IDAT files
#dir.create("data/TODAY_GWAS_DATA", showWarnings = FALSE)
#command = paste("tar", "-xf", temp_file1, "-C", "data/TODAY_GWAS_DATA")
#system(command)

# Unpack CRAM files
dir.create("data/TODAY_WHOLE_EXOMES", showWarnings = FALSE)
command = paste("tar", "-xf", temp_file2, "-C", "data/TODAY_WHOLE_EXOMES")
system(command)


wxs_files = readLines("data/TODAY_WHOLE_EXOMES.list_files.txt")
wxs_files = here(wxs_files[grepl("cram$", wxs_files) == TRUE])

writeLines(wxs_files, con = "data/TODAY_WHOLE_EXOMES.cram_files.txt", sep = "\n")

