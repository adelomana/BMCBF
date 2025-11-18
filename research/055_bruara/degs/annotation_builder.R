rm(list = ls())

#
# -1. install libraries
# 
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("biomaRt")

#
# 0. load libraries
#
library(biomaRt)


#
# 1. generate gene to transcript mapping
#

# Ensembl 108 Oct 2022 https://oct2022.archive.ensembl.org     108
# last time I ran this I had to call fist host = 'https://www.ensembl.org', then it worked
# 2025.11.17, it seems it works now as it is

mart = biomaRt::useMart(biomart="ENSEMBL_MART_ENSEMBL", 
                        dataset="drerio_gene_ensembl",
                        host = 'https://oct2022.archive.ensembl.org', 
                        verbose = TRUE)

# attributes = listAttributes(mart)
# hgnc_symbol gives less than external_gene_name
working_attributes = c('ensembl_gene_id', 
                       'external_gene_name',
                       'gene_biotype',
                       'description')
annotation = biomaRt::getBM(attributes=working_attributes, 
                            mart=mart,
                            verbose=TRUE)
dim(annotation) # Version 108 gives 37241 entries.
View(annotation)

store = paste('/Users/adrian/software/kallisto/zebrafish_index_standard', '/annotation.tsv', sep='')
write.table(annotation, file=store, quote=FALSE, sep='\t', col.names=NA)
