# Annotate VCF files with gene information (snpEff)
- annotate_vcf_genes.R

# Compare QTL signals with outcome/comorbidity signals
- compare_qtl_outcome.R

# Analyze QTLs vs outcomes
- analyze_qtl_outcome.R




# Combine QTL studies from TODAY with publicly available

# Download QTL data:
## eQTLGen
- https://eqtlgen.org/sc/datasets/1m-scbloodnl-eqtls.html
- cis-eQTLs per cell type and condition (genome wide): https://molgenis26.gcc.rug.nl/downloads/1m-scbloodnl/eqtls_20201106_genome_wide.tar.gz
- eQTLs phase 1: https://eqtlgen.org/cis-eqtls.html: https://molgenis26.gcc.rug.nl/downloads/eqtlgen/cis-eqtl/2019-12-11-cis-eQTLsFDR-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz

## eQTL catalogue
- https://www.ebi.ac.uk/eqtl/Data_access/
- https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/
- metadata: https://github.com/eQTL-Catalogue/eQTL-Catalogue-resources/blob/master/data_tables/dataset_metadata.tsv

## Download additional data
- Gencode V.39 (hg38)
