rm(list = ls())

#
# -1. install libraries
# 
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install()
##BiocManager::install("devtools")    # only if devtools not yet installed
#BiocManager::install("pachterlab/sleuth") # for some reason this line needs to be run twice

library(devtools)
library(sleuth)
library(ggplot2)
library(dplyr)
library(stringr)
library(data.table)
library(matrixStats) # required for rowMedians

#
# 0. user-defined variables
#
setwd("/Users/adrian/scratch/")
kallisto_dir = "/Users/adrian/research/bmcbf/054_selfoss/results/profiles"
results_dir = '/Users/adrian/research/bmcbf/054_selfoss/results/degs_sleuth'

# thresholds
count_threshold = 20
effect_size_threshold = log2(2)
tpm_threshold = 2

#
# 1. generate gene to transcript mapping and annotation
#
t2g_file = '/Users/adrian/software/kallisto/human_index_standard/t2g.txt'
t2g = read.csv(t2g_file, sep='\t', header=FALSE)
names(t2g)[names(t2g) == "V1"] <- "target_id"
names(t2g)[names(t2g) == "V2"] <- "ens_gene"
View(t2g)
dim(t2g)

annotation_file = '/Users/adrian/software/kallisto/human_index_standard/annotation.tsv'
full_annotation = read.csv(annotation_file, sep='\t')
annotation <- full_annotation %>% distinct(ensembl_gene_id, .keep_all = TRUE)

#
# 2. define metadata
#
dirnames = list.dirs(kallisto_dir, full.names=TRUE, recursive=FALSE)
dirnames = dirnames[grep('_processed', dirnames)]
paths = file.path(dirnames, 'kallisto_output_c/abundance.h5')
labels = sapply(strsplit(paths, split='/',fixed=TRUE), function(x) (x[9]))
sample = str_remove(labels, '_processed')
print(sample)

metadata = data.frame(sample)
metadata$path = paths

genotypes = c(rep('RES', 3), rep('SEN', 3))
metadata$genotype = genotypes

dim(metadata)
View(metadata)

seta_indexes = 1:3
setb_indexes = 4:6

# make sure to relevel for the appropriate reference

#
# 3. preliminary data filter
#

# prepare object
so = sleuth_prep(metadata,
                 target_mapping = t2g,
                 aggregation_column = 'ens_gene',
                 #transform_fun_counts = function(x) (log2(x+0.5)), # surprisingly this has an effect on significance: 1,906 DEGs without, 1,879 with it. Oh, boy.
                 read_bootstrap_tpm = TRUE)

nrow(so$sample_to_covariates)
length(so$target_mapping$target_id)

dim(so$obs_norm)
dim(so$obs_raw)
dim(so$target_mapping) # these are the targets?

# Convert sleuth object to TPM and count matrix. We are working with transcripts
tpm_matrix = sleuth_to_matrix(so, which_df = "obs_norm", 'tpm')
dim(tpm_matrix)
count_matrix = sleuth_to_matrix(so, which_df='obs_norm', 'est_counts') 
dim(count_matrix)

# filter transcripts that have less than 20 counts in difference
a = count_matrix[ , seta_indexes]
b = count_matrix[ , setb_indexes]
c = rowMedians(a) - rowMedians(b)
keep = abs(c) >= count_threshold
sum(keep)
so$target_mapping <- so$target_mapping[so$target_mapping$target_id %in% rownames(count_matrix[keep,]), ]

nrow(so$sample_to_covariates)
length(so$target_mapping$target_id)

# keep features with at least a max median expression of the TPM threshold
a = rowMedians(tpm_matrix[ , seta_indexes])
b = rowMedians(tpm_matrix[ , setb_indexes])
c = pmax(a, b)
keep = c >= tpm_threshold
sum(keep)
so$target_mapping <- so$target_mapping[so$target_mapping$target_id %in% rownames(count_matrix[keep,]), ]

nrow(so$sample_to_covariates)
length(so$target_mapping$target_id)

#
# 4. contrast 
#

so = sleuth_fit(so, ~genotype, 'full')
so = sleuth_fit(so, ~1, 'reduced')
so = sleuth_lrt(so, 'reduced', 'full')

# using LRT instead of Wald because authors mentioned that it gives lots of false positives
# do not use gene mode in prep, use Lancaster aggregation method for transcripts into genes. 
# see https://pachterlab.github.io/sleuth/docs/sleuth_results.html
sleuth_table = sleuth_results(so, 'reduced:full', 'lrt', show_all = FALSE, pval_aggregate = TRUE) 

# filter simple duplicates
sleuth_significant = dplyr::filter(sleuth_table, qval < 0.05)
dim(sleuth_significant)
anti = dplyr::filter(sleuth_table, qval > 0.05)

# filter table as expected, see https://pachterlab.github.io/sleuth/docs/sleuth_results.html
filtered_df <- sleuth_significant[!duplicated(sleuth_significant$target_id), ] # filtering repetitives
dim(filtered_df)

# without any filter: 1,906 significant DEGs
# with filter of 20 estimated counts difference: 2,784 DEGs
# with filter of counts and TPMs: 2,018 DEGs

#
# 5. get TPMs and log2FC from DESeq2. Oh, well.
# from sleuth_to_matrix: Note this currently does not support returning raw values for gene-level counts or TPMs.
#
#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("DESeq2")
#BiocManager::install("tximport")
library(DESeq2)
library(tximport)

txi = tximport(metadata$path, type="kallisto", tx2gene=t2g)
dds = DESeqDataSetFromTximport(txi, colData=metadata, design=~genotype) 
dds$genotype = relevel(dds$genotype, ref="SEN")
dds = DESeq(dds, test="LRT", reduced=~1)
res = results(dds, parallel=TRUE)

filtered_df$log2FC = res[filtered_df$target_id, "log2FoldChange"]
dim(filtered_df)
final = filtered_df[abs(filtered_df$log2FC) >= 1, ]
dim(final)

subset = txi$abundance[final$target_id, ]
a = rowMedians(subset[ , seta_indexes])
b = rowMedians(subset[ , setb_indexes])
final$TPMaDESeq = a
final$TPMb = b

#
# 6. add annotation
#
ensembl_results_wo = sapply(strsplit(final$target_id, split='.',fixed=TRUE), function(x) (x[1]))
length(ensembl_results_wo)
sub = annotation[annotation$ensembl_gene_id %in% ensembl_results_wo, ]
dim(sub)
sub = sub[, c(3, 4, 5, 6)]
sub$description2 = sapply(strsplit(sub$description, split='[Source',fixed=TRUE), function(x) (x[1]))
rownames(sub) <- sub$ensembl_gene_id

final$ensembl_id = sub[ensembl_results_wo, 'ensembl_gene_id']
final$gene_name = sub[ensembl_results_wo, 'external_gene_name']
final$biotype = sub[ensembl_results_wo, 'gene_biotype']
final$description2 = sub[ensembl_results_wo, 'description2']

plot_pca(so, color_by = 'genotype') 

write.table(final, 
            file = paste(results_dir, '/effect_genotype.tsv', sep=''), 
            sep = '\t',
            quote = FALSE)



