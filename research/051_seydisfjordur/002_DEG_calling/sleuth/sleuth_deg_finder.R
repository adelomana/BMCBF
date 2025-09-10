rm(list = ls())

#
# -1. install libraries
# 
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install()
#BiocManager::install("devtools")    # only if devtools not yet installed
#BiocManager::install("pachterlab/sleuth")

library(devtools)
library(sleuth)
library(ggplot2)
library(dplyr)
library(stringr)
library(data.table)

#
# 0. user-defined variables
#
setwd("/Users/adrian/scratch/")
kallisto_dir = "/Users/adrian/research/bmcbf/051_seydisfjordur/results/profiles"
results_dir = '/Users/adrian/research/bmcbf/051_seydisfjordur/results/degs_sleuth'

#
# 1. generate gene to transcript mapping
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

genotypes = c(rep('SEN', 3), rep('RES', 3))
metadata$genotype = genotypes

dim(metadata)
View(metadata)

# make sure to relevel for the appropriate reference

#
# 3. contrasts
#

# using LRT instead of Wald because authors mentioned that it gives lots of false positives

# prepare contrast
so = sleuth_prep(metadata,
                 target_mapping = t2g,
                 aggregation_column = 'ens_gene',
                 transform_fun_counts = function(x) (log2(x+0.5)),
                 read_bootstrap_tpm = TRUE)


#############\\
#filter_low_expression <- function(row) {
#  mean(row > 5) >= 0.5  # expressed (>5 estimated counts) in at least 50% of samples
#}

#so <- sleuth_prep(s2c, ~condition, target_mapping = t2g, aggregation_column = "gene", 
#                  extra_bootstrap_summary = TRUE,
#                  filter_fun = filter_low_expression)
###########

# contrast 
so = sleuth_fit(so, ~genotype, 'full')
so = sleuth_fit(so, ~1, 'reduced')
so = sleuth_lrt(so, 'reduced', 'full')

# do not use gene mode in prep, use Lancaster aggregation method for transcripts into genes. 
# see https://pachterlab.github.io/sleuth/docs/sleuth_results.html
sleuth_table = sleuth_results(so, 'reduced:full', 'lrt', show_all = FALSE, pval_aggregate = TRUE) 


#! implement filters on ammount, etc
#! get log2fc from the wald test
#! add annotation, save final table
#! check the relevel


# filter
sleuth_significant = dplyr::filter(sleuth_table, qval < 0.05)
dim(sleuth_significant)
anti = dplyr::filter(sleuth_table, qval > 0.05)



# filter table as expected, see https://pachterlab.github.io/sleuth/docs/sleuth_results.html
filtered_df <- sleuth_significant[!duplicated(sleuth_significant$target_id), ] # filtering repetitives
dim(filtered_df)

plot_pca(so, color_by = 'genotype') 

write.table(filtered_df, 
            file = paste(results_dir, '/effect_genotype.tsv', sep=''), 
            sep = '\t',
            quote = FALSE)
write.table(anti, 
            file = paste(results_dir, '/effect_genotype.anti.tsv', sep=''), 
            sep = '\t',
            quote = FALSE)


