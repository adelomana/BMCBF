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
kallisto_dir = "/Users/adrian/research/bmcbf/051_seydisfjordur/results/profiles"
results_dir = '/Users/adrian/research/bmcbf/051_seydisfjordur/results/degs_deseq'

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
# 2. define metadata
#
dirnames = list.dirs(kallisto_dir, full.names=TRUE, recursive=FALSE)
dirnames = dirnames[grep('_processed', dirnames)]
paths = file.path(dirnames, 'kallisto_output_c/abundance.h5')
labels = sapply(strsplit(paths, split='/',fixed=TRUE), function(x) (x[9]))
labels = str_remove(labels, '_processed')
print(labels)

metadata = data.frame(labels)
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
count_threshold = 20
effect_size_threshold = log2(2)
tpm_threshold = 2

#
# 3.1. contrast 
#
txi = tximport(metadata$path, type="kallisto", tx2gene=t2g)
dds = DESeqDataSetFromTximport(txi, colData=metadata, design=~genotype) 
dds$genotype = relevel(dds$genotype, ref="SEN")

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
cat(blue(paste('RES vs SEN:', dim(filtred_results)[1], sep=' ')), fill=TRUE)
write.table(sorted_filtred_results, file=paste(results_dir, '/effect_RES_vs_SEN.raw.tsv', sep=''), quote=FALSE, sep='\t')


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

write.table(sorted_filtred_results, file=paste(results_dir, '/effect_RES_vs_SEN.for.tsv', sep=''), quote=FALSE, sep='\t')
write.table(anti_results, file=paste(results_dir, '/effect_RES_vs_SEN.anti.tsv', sep=''), quote=FALSE, sep='\t')

plotPCA(rlog(dds), intgroup=c('genotype')) + ggtitle('effect RES vs SEN')

#               
# volcano
#
plotting_x = sorted_filtred_results$log2FoldChange
y = sorted_filtred_results$padj
epsilon = min(y[y !=0])
plotting_y = -log10(y + epsilon) 
wo = sapply(strsplit(rownames(txi$abundance), split='.',fixed=TRUE), function(x) (x[1]))
rownames(txi$abundance) = wo
z = log10(rowMedians(txi$abundance[rownames(sorted_filtred_results), ]) + 1)
df = data.frame(plotting_x=clip(plotting_x, .min=-6, .max=6), plotting_y=clip(plotting_y, .min=0, .max=20), plotting_z=clip(z, .min=0, .max=3))
reds = df[df$plotting_x > 0, ]
blues = df[df$plotting_x < 0, ]

plotting_x = anti_results$log2FoldChange
plotting_y = -log10(anti_results$padj)
blacks = data.frame(plotting_x=clip(plotting_x, .min=-6, .max=6), plotting_y=clip(plotting_y, .min=0, .max=20))

ggplot() + 
  geom_point(data=reds, aes(x=plotting_x, y=plotting_y, color=plotting_z), , size=3, shape=19, alpha=2/3, stroke=0) + 
  geom_point(data=blues, aes(plotting_x, plotting_y, color=plotting_z), size=3, shape=19, alpha=2/3, stroke=0) +
  geom_point(data=blacks, aes(plotting_x, plotting_y), color = "black", size=1, shape=19, alpha=0.2, stroke=0) +
  labs(x=expression('Expression [log'[2]~'FC]'), y=expression('Significance [log'[10]~'adjusted P]'), color=expression('Expression [log'[10]~'TPM]')) + 
  theme_linedraw() +
  geom_segment(aes(x=-1, xend=-1, y=-log10(0.05), yend=20), linetype=2) +
  geom_segment(aes(x=1, xend=1, y=-log10(0.05), yend=20), linetype=2) +
  geom_segment(aes(x=-6, xend=-1, y=-log10(0.05), yend=-log10(0.05)), linetype=2) +
  geom_segment(aes(x=1, xend=6, y=-log10(0.05), yend=-log10(0.05)), linetype=2) +
  xlim(-6, 6) +
  scale_color_viridis_c(option = "cividis") 
#ggsave(paste('effect_T0570_vs_T84', '.png', sep=''))