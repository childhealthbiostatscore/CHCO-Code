setwd('G:/My Drive/lab files/Bryan Bergman/cohort 11-8-24 CAT VAT SAT')


file_list = as.data.frame(list.files('./raw data/20240830_LH00407_0078_B22NV2JLT3/'))
colnames(file_list) = 'file_name'
file_list
file_list = as.data.frame(file_list[!grepl('20240830_LH00407_0078_B22NV2JLT3_L123_Bergman_demux.csv', file_list$file_name),])

colnames(file_list) = 'file_name'
file_list
file_list$sample = gsub('_R1_001.fastq.gz', '', file_list$file_name)
file_list$sample = gsub('_R2_001.fastq.gz', '', file_list$sample)

sample_list = as.data.frame(unique(file_list$sample))

colnames(sample_list) = 'sample'
head(sample_list)
sample_list$sub_command = paste0('kallisto quant -i /dfs4/som/mseldin/genome_files/kallisto_human_index/Homo_sapiens.GRCh38.cdna.all.release-94_k31.idx -o /dfs4/som/mseldin/bergman_crosstalk/set2/kallisto_outs/', sample_list$sample, ' /dfs4/som/mseldin/bergman_crosstalk/set2/20240830_LH00407_0078_B22NV2JLT3/',  sample_list$sample, '_R1_001.fastq.gz ', ' /dfs4/som/mseldin/bergman_crosstalk/set2/20240830_LH00407_0078_B22NV2JLT3/', sample_list$sample, '_R2_001.fastq.gz')
sample_list$sub_command[1]
write.csv(sample_list$sub_command, file = 'kallisto_run.csv', row.names = F, quote = F)

Mus_musculus.GRCm39.cdna.idx
###############################
#assemble into counts matrix
setwd('G:/My Drive/lab files/Bryan Bergman/pancreatic islets analysis/kallisto tpms')
sample_list_k = list.files('./kallisto_outs/')

get_tpm_from_files = function(sample_name){
  tpm_file = read.delim(paste0('./kallisto_outs/', sample_name, '/abundance.tsv'))
  tpm_file$sample_ID = paste0(sample_name)
  df_name = paste0(sample_name, '_tpms')
  return(data.frame(assign(  df_name, tpm_file )))
}
sample_list_k
#"1"  "10" "12" "14" "15" "16" "17" "18" "2"  "21" "23" "24" "25" "26" "28" "29" "4"  "5"  "6"  "7"  "9"

d1 = get_tpm_from_files("1")
d10 = get_tpm_from_files("10")
d12 = get_tpm_from_files("12")
d14 = get_tpm_from_files("14")

#d15 = get_tpm_from_files("15")
d16 = get_tpm_from_files("16")
d17 = get_tpm_from_files("17")
d18 = get_tpm_from_files("18")

d2 = get_tpm_from_files("2")
d21 = get_tpm_from_files("21")
#d23= get_tpm_from_files("23")
#d24 = get_tpm_from_files("24")

d25 = get_tpm_from_files("25")
d26 = get_tpm_from_files("26")
#d28= get_tpm_from_files("28")
d29 = get_tpm_from_files("29")

d4 = get_tpm_from_files("4")
d5 = get_tpm_from_files("5")
d6= get_tpm_from_files("6")
d7 = get_tpm_from_files("7")

d9 = get_tpm_from_files("9")

full_melted_tpms = as.data.frame(rbind(d1,  d10, d12, d14,  d16, d17, d18, d2,  d21,  d25, d26, d29, d4,  d5,  d6,  d7,  d9))
write.csv(full_melted_tpms, file = 'melted_tpms unfiltered.csv', row.names = F)

head(full_melted_tpms)
tpm_matrix = dcast(full_melted_tpms, target_id ~ sample_ID, value.var = 'tpm', fun.aggregate = mean)       
write.csv(tpm_matrix, file = 'tpm matrix by sample.csv', row.names = F)