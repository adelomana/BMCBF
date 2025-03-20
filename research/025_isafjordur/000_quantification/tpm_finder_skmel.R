rm(list = ls())

#
# -1. install libraries
# 
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("DESeq2")
# BiocManager::install("tximport")

#
# 0. load libraries
#
library(DESeq2)
library(tximport)
library(biomaRt)
library(stringr)

#
# 1. user-defined variables
#
setwd("~/scratch/")
kallisto_dir = "/Users/adrian/research/bmcbf/025_isafjordur/results/000_quantification"
results_dir = '/Users/adrian/research/bmcbf/025_isafjordur/results/000_quantification'
a_cases = c('SkMel28-MITFKO_ev_skmel28_rep1', 'SkMel28-MITFKO_ev_skmel28_rep2', 'SkMel28-MITFKO_ev_skmel28_rep3', 'SkMel28-MITFKO_ev_skmel28_rep4', 'SkMel28-MITFKO_mitf_x6_rep1', 'SkMel28-MITFKO_mitf_x6_rep3', 'SkMel28-MITFKO_mitf_x6_rep4')

#
# 1. generate gene to transcript mapping
#
mart = biomaRt::useMart(biomart="ENSEMBL_MART_ENSEMBL", 
                        dataset="hsapiens_gene_ensembl",
                        host = 'https://oct2022.archive.ensembl.org', 
                        # Ensembl 108 Oct 2022 https://oct2022.archive.ensembl.org     108
                        # last time I ran this I had to call fist host = 'https://www.ensembl.org', then it worked
                        verbose = TRUE)
# attributes = listAttributes(mart)
# hgnc_symbol gives less than external_gene_name
working_attributes = c('ensembl_transcript_id', 
                      'ensembl_gene_id', 
                      'external_gene_name',
                      'gene_biotype',
                      'description')
t2g = biomaRt::getBM(attributes=working_attributes, 
                     mart=mart,
                     verbose=TRUE)
dim(t2g) # Version 108 gives 274081 entries. Version 113 gives 412034
View(t2g)

#
# 2. define metadata
#
dirnames = list.dirs(kallisto_dir, full.names=TRUE, recursive=FALSE)
dirnames = dirnames[grep('SkMel28', dirnames)]
paths = file.path(dirnames, 'kallisto_output_c/abundance.h5')
for (i in 1:length(paths))
{
  print(i)
  print(paths[i])
  for (j in 1:length(a_cases))
  {
    
    if (grepl(a_cases[j], paths[i]) == TRUE) {
      paths[i] = sub('output_c', 'output_a', paths[i])
    }
  }
} 
labels = sapply(strsplit(paths, split='/',fixed=TRUE), function(x) (x[9]))
labels = str_remove(labels, '_processed')

metadata = data.frame(labels)
metadata$path = paths
View(metadata)

#
# 3. read files
#
txi = tximport(metadata$path, type="kallisto", tx2gene=t2g, ignoreTxVersion=TRUE)

#
# 4. find abundance
#
tpm = txi$abundance
colnames(tpm) = metadata$labels
dim(tpm)
View(tpm)

#
# 5. store
#
store = paste(results_dir, '/DESeq2_TPM_values.skmel.tsv', sep='')
write.table(tpm, file=store, quote=FALSE, sep='\t', col.names=NA)
