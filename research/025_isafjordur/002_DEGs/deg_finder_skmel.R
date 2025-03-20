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
library(biomaRt)
library(BiocParallel)
library(crayon) 
library(ggplot2)
library(stringr)
library(ramify) # this is for clip, for the volcano

#
# 0. user-defined variables
#
setwd("~/scratch/")
kallisto_dir = "/Users/adrian/research/bmcbf/025_isafjordur/results/000_quantification"
results_dir = '/Users/adrian/research/bmcbf/025_isafjordur/results/002_DEGs'
a_cases = c('SkMel28-MITFKO_ev_skmel28_rep1', 'SkMel28-MITFKO_ev_skmel28_rep2', 'SkMel28-MITFKO_ev_skmel28_rep3', 'SkMel28-MITFKO_ev_skmel28_rep4', 'SkMel28-MITFKO_mitf_x6_rep1', 'SkMel28-MITFKO_mitf_x6_rep3', 'SkMel28-MITFKO_mitf_x6_rep4')

#
# 1. generate gene to transcript mapping
#
mart = biomaRt::useMart(biomart="ENSEMBL_MART_ENSEMBL", 
                        dataset="hsapiens_gene_ensembl",
                        host = 'https://oct2022.archive.ensembl.org', # because of 108
                        verbose = TRUE)
# attributes = listAttributes(mart)
working_attributes = c('ensembl_transcript_id', 
                       'ensembl_gene_id', 
                       'external_gene_name',
                       'gene_biotype',
                       'description')
t2g = biomaRt::getBM(attributes=working_attributes, 
                     mart=mart,
                     verbose=TRUE)
dim(t2g)

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

replicates = c(c('A', 'B', 'C', 'D'), c('A', 'B', 'C'))
genotypes = c(rep('wt', 4), rep('ko', 3))

metadata = data.frame(labels)
metadata$path = paths
metadata$replicate = replicates
metadata$genotype = genotypes
View(metadata)

#
# 3. contrasts
#
threshold = 20
effect_size_threshold = log2(2)
tpm_threshold = 1

#
# 3.1. contrast 
#
txi = tximport(metadata$path, type="kallisto", tx2gene=t2g, ignoreTxVersion=TRUE)
dds = DESeqDataSetFromTximport(txi, colData=metadata, design=~genotype) 
dds$time = relevel(dds$genotype, ref="wt")

# keep features with at least 20 counts median difference
cat(blue(paste('size before counts filtering:', dim(dds)[1], sep=' ')), fill=TRUE)
a = counts(dds)[ , 1:4]
b = counts(dds)[ , 5:7]
c = rowMedians(a) - rowMedians(b)
keep = abs(c) >= threshold
dds = dds[keep, ]
cat(blue(paste('size after counts filtering:', dim(dds)[1], sep=' ')), fill=TRUE)

# keep features with at least a max median expression of 1 TPM.
cat(blue(paste('size before counts filtering:', dim(dds)[1], sep=' ')), fill=TRUE)
subset = txi$abundance[names(dds), ]
a = rowMedians(subset[ , 1:4])
b = rowMedians(subset[ , 5:7])
c = pmax(a, b)
keep = c >= tpm_threshold
dds = dds[keep, ]
cat(blue(paste('size after counts filtering:', dim(dds)[1], sep=' ')), fill=TRUE)

dds = DESeq(dds, test="LRT", reduced=~1)

res = results(dds, parallel=TRUE, alpha=0.05) # it does not seem to affect  https://www.biostars.org/p/209118/ 
filtred_results = res[which(res$padj < 0.05 & abs(res$log2FoldChange) > effect_size_threshold), ]
sorted_filtred_results = filtred_results[order(filtred_results[["padj"]]),]
anti_results = res[which(res$padj > 0.05 | abs(res$log2FoldChange) < effect_size_threshold), ]
cat(blue(paste('contrast wt vs control:', dim(filtred_results)[1], sep=' ')), fill=TRUE)
write.table(sorted_filtred_results, file=paste(results_dir, '/effect_ko_vs_wt.skmel.for.tsv', sep=''), quote=FALSE, sep='\t')
write.table(anti_results, file=paste(results_dir, '/effect_ko_vs_wt.skmel.anti.tsv', sep=''), quote=FALSE, sep='\t')

plotPCA(rlog(dds), intgroup=c('genotype')) + ggtitle('effect ko vs wt')

#               
# volcano
#
plotting_x = sorted_filtred_results$log2FoldChange
y = sorted_filtred_results$padj
epsilon = min(y[y !=0])
plotting_y = -log10(y + epsilon) 
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
  labs(x=expression('Expression [log'[2]~'FC]'), y=expression('Significance [log'[10]~'adjusted P]'), color=expression('Expression [log'[10]~'TPM]'), title='KO vs WT | skmel') + 
  theme_linedraw() +
  geom_segment(aes(x=-1, xend=-1, y=-log10(0.05), yend=20), linetype=2) +
  geom_segment(aes(x=1, xend=1, y=-log10(0.05), yend=20), linetype=2) +
  geom_segment(aes(x=-6, xend=-1, y=-log10(0.05), yend=-log10(0.05)), linetype=2) +
  geom_segment(aes(x=1, xend=6, y=-log10(0.05), yend=-log10(0.05)), linetype=2) +
  xlim(-6, 6) +
  scale_color_viridis_c(option = "cividis") 
#ggsave(paste(label, '.png', sep=''))