rm(list = ls())

#
# -1. install libraries
# 
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("DESeq2")

#
# 0. load libraries
#
library(DESeq2)
library(tximport)
library(crayon) 
library(ggplot2)
library(stringr)
library(ramify) # this is for clip, for the volcano
library(dplyr)

#
# 0. user-defined variables
#
setwd("~/scratch/")
kallisto_dir = "/Users/adrian/research/bmcbf/061_drangajokull/profiles"
results_dir = '/Users/adrian/research/bmcbf/061_drangajokull/degs'

#
# 1. generate gene to transcript mapping
#
t2g_file = '/Users/adrian/software/kallisto/human_index_standard/t2g.txt'
t2g = read.csv(t2g_file, sep='\t', header=FALSE)
View(t2g)
dim(t2g)

annotation_file = '/Users/adrian/software/kallisto/human_index_standard/annotation.tsv'
full_annotation = read.csv(annotation_file, sep='\t')
annotation <- full_annotation %>% distinct(ensembl_gene_id, .keep_all = TRUE)

#
# 2. contrast thresholds
#
count_threshold = 20
effect_size_threshold = log2(1.5) 
tpm_threshold = 2

#
# 3. define full metadata
#
dirnames = list.dirs(kallisto_dir, full.names=TRUE, recursive=FALSE)
dirnames = dirnames[grep('_processed', dirnames)]
labels = sapply(strsplit(dirnames, split='/',fixed=TRUE), function(x) (x[8]))
labels = str_remove(labels, '_processed')
metadata = data.frame(labels)

paths = file.path(dirnames, paste0("kallisto_output_", labels, "_unstranded"), 'abundance.h5')
metadata$path = paths

metadata$genotype = c(rep('C1', 6), rep('C4', 6), rep('C0', 6), rep('P', 6))
metadata$treatment = rep(c(rep('si', 3), rep('NT', 3)), 4)
View(metadata)

seta_indexes = 1:3
setb_indexes = 4:6

#
# 4. iterate contrasts
#
contrasts = list(
  c('C_1_8_siC', 'C_Ctrl_siC', 'genotype'),
  c('C_4_2_siC', 'C_Ctrl_siC', 'genotype'),
  c('P_siC', 'C_Ctrl_siC', 'genotype'),
  
  c('C_4_2_siC', 'C_1_8_siC', 'genotype'),
  
  c('C_Ctrl_siA', 'C_Ctrl_siC', 'treatment'),
  c('C_1_8_siA', 'C_1_8_siC', 'treatment'),
  c('C_4_2_siA', 'C_4_2_siC', 'treatment'),
  c('P_siA', 'P_siC', 'treatment')
)

for (contrast in contrasts) {
  sample_flag = contrast[1]
  control_flag = contrast[2]
  design_name    <- contrast[3] 
  print(c('working with', sample_flag, control_flag, design_name))
  
  # 3.1. define working metadata
  rules = grepl(sample_flag, metadata$labels) | grepl(control_flag, metadata$labels)
  working_metadata <- metadata[rules, ]
  print('working metadata')
  print(working_metadata)
  
  # 3.2. read quantification files
  txi = tximport(working_metadata$path, type="kallisto", tx2gene=t2g)
  design_formula <- as.formula(paste("~", design_name))
  dds = DESeqDataSetFromTximport(txi, colData=working_metadata, design=design_formula) 
  map <- c(genotype = "C0", treatment = "NT")
  reference <- map[design_name]
  if (control_flag == "C_1_8_siC" && design_name == "genotype") {
    reference <- "C1"
  }
  print(c('and the reference is', reference))
  dds[[ design_name ]] <- relevel(dds[[ design_name ]], ref = reference)
  
  # keep features with at least 20 counts median difference
  cat(blue(paste('size before counts filtering:', dim(dds)[1], sep=' ')), fill=TRUE)
  a = counts(dds)[ , seta_indexes]
  b = counts(dds)[ , setb_indexes]
  c = rowMedians(a) - rowMedians(b)
  keep = abs(c) >= count_threshold
  dds = dds[keep, ]
  cat(blue(paste('size after counts filtering:', dim(dds)[1], sep=' ')), fill=TRUE)
  
  # keep features with at least a max median expression of the TPM threshold
  cat(blue(paste('size before counts filtering:', dim(dds)[1], sep=' ')), fill=TRUE)
  subset = txi$abundance[names(dds), ]
  a = rowMedians(subset[ , seta_indexes])
  b = rowMedians(subset[ , setb_indexes])
  c = pmax(a, b)
  keep = c >= tpm_threshold
  dds = dds[keep, ]
  cat(blue(paste('size after counts filtering:', dim(dds)[1], sep=' ')), fill=TRUE)
  
  dds = DESeq(dds, test="LRT", reduced=~1)
  
  res = results(dds, parallel=TRUE, alpha=0.05) # it does not seem to affect  https://www.biostars.org/p/209118/ 
  filtred_results = res[which(res$padj < 0.05 & abs(res$log2FoldChange) > effect_size_threshold), ]
  sorted_filtred_results = filtred_results[order(filtred_results[["padj"]]),]
  anti_results = res[which(res$padj > 0.05 | abs(res$log2FoldChange) < effect_size_threshold), ]
  cat(blue(paste('contrast rank:', dim(filtred_results)[1], sep=' ')), fill=TRUE)
  write.table(sorted_filtred_results, file=paste(results_dir, '/', paste(sample_flag, control_flag, sep='_'), '.raw.tsv', sep=''), quote=FALSE, sep='\t')
  
  ensembl_results_wo = sapply(strsplit(rownames(sorted_filtred_results), split='.',fixed=TRUE), function(x) (x[1]))
  rownames(sorted_filtred_results) = ensembl_results_wo
  length(ensembl_results_wo)
  sub = annotation[annotation$ensembl_gene_id %in% ensembl_results_wo, ]
  dim(sub)
  sub = sub[, c(3, 4, 5, 6)]
  sub$description2 = sapply(strsplit(sub$description, split='[Source',fixed=TRUE), function(x) (x[1]))
  rownames(sub) <- sub$ensembl_gene_id
  
  sorted_filtred_results$ensembl_id = sub[rownames(sorted_filtred_results), 'ensembl_gene_id']
  sorted_filtred_results$gene_name = sub[rownames(sorted_filtred_results), 'external_gene_name']
  sorted_filtred_results$biotype = sub[rownames(sorted_filtred_results), 'gene_biotype']
  sorted_filtred_results$description2 = sub[rownames(sorted_filtred_results), 'description2']
  
  # this is very dangerous, but lets go
  wo = sapply(strsplit(rownames(subset), split='.',fixed=TRUE), function(x) (x[1]))
  rownames(subset) = wo
  sorted_filtred_results$medianTPM_WT = rowMedians(subset[rownames(sorted_filtred_results), seta_indexes]) # SEN
  sorted_filtred_results$medianTPM_HET = rowMedians(subset[rownames(sorted_filtred_results), setb_indexes]) # RES
  
  write.table(sorted_filtred_results, file=paste(results_dir, '/', paste(sample_flag, control_flag, sep='_'), '.for.tsv', sep=''), quote=FALSE, sep='\t')
  write.table(anti_results, file=paste(results_dir, '/', paste(sample_flag, control_flag, sep='_'), '.anti.tsv', sep=''), quote=FALSE, sep='\t')
  
}