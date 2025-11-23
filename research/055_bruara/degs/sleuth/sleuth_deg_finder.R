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
setwd("~/scratch/")
kallisto_dir = "/Users/adrian/research/bmcbf/055_bruara/results/transcriptomics/quantification"
results_dir = '/Users/adrian/research/bmcbf/055_bruara/results/transcriptomics/degs_sleuth'

# thresholds
count_threshold = 20
effect_size_threshold = log2(2)
tpm_threshold = 2

#
# 1. generate gene to transcript mapping and annotation
#
t2g_file = '/Users/adrian/software/kallisto/zebrafish_index_standard/t2g.txt'
t2g = read.csv(t2g_file, sep='\t', header=FALSE)
names(t2g)[names(t2g) == "V1"] <- "target_id"
names(t2g)[names(t2g) == "V2"] <- "ens_gene"
View(t2g)
dim(t2g)

annotation_file = '/Users/adrian/software/kallisto/zebrafish_index_standard/annotation.tsv'
full_annotation = read.csv(annotation_file, sep='\t')
annotation <- full_annotation %>% distinct(ensembl_gene_id, .keep_all = TRUE)
View(annotation)

#
# 2. define metadata
#
dirnames = list.dirs(kallisto_dir, full.names=TRUE, recursive=FALSE)
dirnames = dirnames[grep('_processed', dirnames)]
paths = file.path(dirnames, 'kallisto_output_a/abundance.h5')
labels = sapply(strsplit(paths, split='/',fixed=TRUE), function(x) (x[10]))
sample = str_remove(labels, '_processed')
print(sample)

metadata = data.frame(sample)
metadata$path = paths

genotypes = c(rep('DMSO', 3), rep('R', 3), rep('S', 3))
metadata$genotype = genotypes

#metadata = metadata[1:6, ]
metadata = metadata[c(1,2,3,7,8,9),]

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
                 #transform_fun_counts = function(x) (log2(x+0.5)), # surprisingly this has an effect on significance. Oh, boy.
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

# R gives zero DEGs
# S gives zero DEGs


