
infiles  = list.files("/mmfs1/home/mdanto/projects/TODAY/data/TODAY_GWAS_DATA", full.names = TRUE)
outfiles = sub("/mmfs1/home/mdanto/projects/TODAY/data/TODAY_GWAS_DATA", "/mmfs1/gscratch/togo/matteo/projects/TODAY_array/data/array", infiles)

to_convert = data.frame(infile = infiles, outfile = outfiles)

invisible(lapply(1:nrow(to_convert), function(ii)
{
  file.copy(from = to_convert[ii, "infile"], 
            to   = to_convert[ii, "outfile"])
}))
