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



#
# 0. user-defined variables
#
setwd("~/scratch/")
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

seta_indexes = 1:3
setb_indexes = 4:6

#
# 3. contrasts
#

# old, need to check
# using LRT instead of Wald because Wald gives three times more and authors mentioned that it gives lots of false positives
# Using pval_aggregate = TRUE bc otherwise no DEGs. Get FC from est counts file


# prepare contrast
so = sleuth_prep(metadata,
                 target_mapping = t2g,
                 aggregation_column = 'ens_gene',
                 read_bootstrap_tpm = TRUE)
# contrast 
so = sleuth_fit(so, ~genotype, 'full')
so = sleuth_fit(so, ~1, 'reduced')
so = sleuth_lrt(so, 'reduced', 'full')
sleuth_table = sleuth_results(so, 
                              'reduced:full', 
                              'lrt',
                              show_all = FALSE,
                              pval_aggregate = TRUE)
sleuth_significant = dplyr::filter(sleuth_table, qval < 0.05)
anti = dplyr::filter(sleuth_table, qval > 0.05)
dim(sleuth_significant)
plot_pca(so, color_by = 'time') + ggtitle('effect time for 2D')
ggsave(file.path(results_dir, 'effect_time_2D.png'))
write.table(sleuth_significant, 
            file = paste(results_dir, '/effect_time_2D.tsv', sep=''), 
            sep = '\t',
            quote = FALSE)
write.table(anti, 
            file = paste(results_dir, '/effect_time_2D.anti.tsv', sep=''), 
            sep = '\t',
            quote = FALSE)

#
# 3.2. contrast effect of time for 3D
#
rule = metadata$culture == '3D'
s2c = metadata[rule, ]
dim(s2c)
View(s2c)
# prepare contrast
so = sleuth_prep(s2c,
                 target_mapping = t2g,
                 aggregation_column = 'ens_gene',
                 read_bootstrap_tpm = TRUE)
# contrast 
so = sleuth_fit(so, ~time, 'full')
so = sleuth_fit(so, ~1, 'reduced')
so = sleuth_lrt(so, 'reduced', 'full')
sleuth_table = sleuth_results(so, 
                              'reduced:full', 
                              'lrt',
                              show_all = FALSE,
                              pval_aggregate = TRUE)
sleuth_significant = dplyr::filter(sleuth_table, qval < 0.05)
anti = dplyr::filter(sleuth_table, qval > 0.05)
dim(sleuth_significant)
plot_pca(so, color_by = 'time') + ggtitle('effect time for 3D')
ggsave(file.path(results_dir, 'effect_time_3D.png'))
write.table(sleuth_significant, 
            file = paste(results_dir, '/effect_time_3D.tsv', sep=''), 
            sep = '\t',
            quote = FALSE)
write.table(anti, 
            file = paste(results_dir, '/effect_time_3D.anti.tsv', sep=''), 
            sep = '\t',
            quote = FALSE)

#
# 3.3. contrast effect of culture at 2 days
#
rule = metadata$time == 'two'
s2c = metadata[rule, ]
dim(s2c)
View(s2c)
# prepare contrast
so = sleuth_prep(s2c,
                 target_mapping = t2g,
                 aggregation_column = 'ens_gene',
                 read_bootstrap_tpm = TRUE)
# contrast 
so = sleuth_fit(so, ~culture, 'full')
so = sleuth_fit(so, ~1, 'reduced')
so = sleuth_lrt(so, 'reduced', 'full')
sleuth_table = sleuth_results(so, 
                              'reduced:full', 
                              'lrt',
                              show_all = FALSE,
                              pval_aggregate = TRUE)
sleuth_significant = dplyr::filter(sleuth_table, qval < 0.05)
anti = dplyr::filter(sleuth_table, qval > 0.05)
dim(sleuth_significant)
plot_pca(so, color_by = 'culture') + ggtitle('effect culture at 2 days')
ggsave(file.path(results_dir, 'effect_culture_2days.png'))
write.table(sleuth_significant, 
            file = paste(results_dir, '/effect_culture_day2.tsv', sep=''), 
            sep = '\t',
            quote = FALSE)
write.table(anti, 
            file = paste(results_dir, '/effect_culture_day2.anti.tsv', sep=''), 
            sep = '\t',
            quote = FALSE)

#
# 3.4. contrast effect of culture at 14 days
#
rule = metadata$time == 'fourteen'
s2c = metadata[rule, ]
dim(s2c)
View(s2c)
# prepare contrast
so = sleuth_prep(s2c,
                 target_mapping = t2g,
                 aggregation_column = 'ens_gene',
                 read_bootstrap_tpm = TRUE)
# contrast 
so = sleuth_fit(so, ~culture, 'full')
so = sleuth_fit(so, ~1, 'reduced')
so = sleuth_lrt(so, 'reduced', 'full')
sleuth_table = sleuth_results(so, 
                              'reduced:full', 
                              'lrt',
                              show_all = FALSE,
                              pval_aggregate = TRUE)
sleuth_significant = dplyr::filter(sleuth_table, qval < 0.05)
anti = dplyr::filter(sleuth_table, qval > 0.05)
dim(sleuth_significant)
plot_pca(so, color_by = 'culture') + ggtitle('effect culture at 14 days')
ggsave(file.path(results_dir, 'effect_culture_day14.png'))
write.table(sleuth_significant, 
            file = paste(results_dir, '/effect_culture_day14.tsv', sep=''), 
            sep = '\t',
            quote = FALSE)
write.table(anti, 
            file = paste(results_dir, '/effect_culture_day14.anti.tsv', sep=''), 
            sep = '\t',
            quote = FALSE)

#
# 3.5 interaction
#
s2c = metadata
dim(s2c)
View(s2c)
# prepare contrast
so = sleuth_prep(s2c,
                 target_mapping = t2g,
                 aggregation_column = 'ens_gene',
                 read_bootstrap_tpm = TRUE)
# contrast 
so = sleuth_fit(so, ~time+culture+time:culture, 'full')
so = sleuth_fit(so, ~time+culture, 'reduced')
so = sleuth_lrt(so, 'reduced', 'full')
sleuth_table = sleuth_results(so, 
                              'reduced:full', 
                              'lrt',
                              show_all = FALSE,
                              pval_aggregate = TRUE)
sleuth_significant = dplyr::filter(sleuth_table, qval <= 0.05)
dim(sleuth_significant)
plot_pca(so, color_by = 'culture') + ggtitle('interaction')
ggsave(file.path(results_dir, 'interaction.png'))
write.table(sleuth_significant, 
            file = paste(results_dir, '/interaction.tsv', sep=''), 
            sep = '\t',
            quote = FALSE)




