rm(list = ls())

#
# -1. install libraries
# 
# if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
# BiocManager::install("vsn")
# BiocManager::install("tximport")

#
# 0. load libraries
#
library(ggplot2)
library(DESeq2)
library(stringr)
library(dplyr)
library(vsn)
library(crayon) 

#
# 1. user-defined variables
#
setwd("~/scratch/")
counts_file = '/Users/adrian/research/bmcbf/031_reykjavik/data/raw_counts_LVRV.tsv'
metadata_file = '/Users/adrian/research/bmcbf/031_reykjavik/metadata/E-MTAB-13553.sdrf.txt'
results_dir = '/Users/adrian/research/bmcbf/031_reykjavik/results/000_differential'

effect_size_threshold = log2(2)
face_value_threshold = 25

#
# 2. read counts
#
cts = as.matrix(read.csv(counts_file, sep="\t"))

#
# 3. define metadata
#
metadata = read.csv(metadata_file, sep="\t")
rownames(metadata) <- metadata$Source.Name
metadata <- metadata[,c("Characteristics.organism.part.","Characteristics.disease.")]
metadata$axis <- factor(metadata$Characteristics.organism.part.)
metadata$health <- factor(metadata$Characteristics.disease.)
sorted_indexes = match(colnames(cts), rownames(metadata))
metadata = metadata[sorted_indexes, ]
metadata %>% mutate(across(where(is.character), str_remove_all, pattern = fixed(" ")))
View(metadata)
levels(metadata$axis)
levels(metadata$health)

metadata = metadata[metadata$Characteristics.disease. == 'heart failure with reduced ejection fraction', ]
dim(metadata)

#
# 2. subset counts
#
dim(cts)
cts = cts[ ,colnames(cts) %in% rownames(metadata)]
dim(cts)
cts = as.matrix(cts)
View(cts)

#
# 4. create working object
#
dds = DESeqDataSetFromMatrix(countData=cts, colData=metadata, design=~axis)
dds$axis = relevel(dds$axis, ref="heart left ventricle")

#
# 5. filter out genes with a median difference smaller than a threshold of reads
#
cat(blue(paste('size before counts filtering:', dim(dds)[1], sep=' ')), fill=TRUE)
others = metadata$Characteristics.organism.part. == 'heart left ventricle'
a = colnames(cts)[others]
reference = metadata$Characteristics.organism.part. == 'heart right ventricle'
b = colnames(cts)[reference]
suba = cts[ , colnames(cts) %in% a]
subb = cts[ , colnames(cts) %in% b]
c(dim(suba), dim(subb))
c = rowMedians(suba) - rowMedians(subb)
keep = abs(c) >= face_value_threshold
dds = dds[keep, ]
cat(blue(paste('size after counts filtering:', dim(dds)[1], sep=' ')), fill=TRUE)

#
# 6. perform test
#
dds = DESeq(dds, test="LRT", reduced=~1)

#
# 7. select DEGs
#
res = results(dds, parallel=TRUE, alpha=0.05) # it does not seem to affect  https://www.biostars.org/p/209118/ 
filtred_results = res[which(res$padj < 0.05 & abs(res$log2FoldChange) > effect_size_threshold), ]
sorted_filtred_results = filtred_results[order(filtred_results[["padj"]]),]
anti_results = res[which(res$padj > 0.05 | abs(res$log2FoldChange) < effect_size_threshold), ]
cat(blue(paste('contrast right vs left:', dim(filtred_results)[1], sep=' ')), fill=TRUE)
plotPCA(rlog(dds), intgroup=c('axis')) + ggtitle('effect right vs left')

#
# 8. write table with annotation
#
t2g = read.csv('/Users/adrian/software/kallisto/human_index_standard/annotation.tsv', sep='\t')
dim(t2g)

sub = t2g[t2g$ensembl_gene_id %in% rownames(sorted_filtred_results), ]
sub = sub[, c(3, 4, 5, 6)]
subu <- sub[!duplicated(sub), ]
rownames(subu) <- subu$ensembl_gene_id

sorted_filtred_results$ensembl_id = subu[rownames(sorted_filtred_results), 'ensembl_gene_id']
sorted_filtred_results$gene_name = subu[rownames(sorted_filtred_results), 'external_gene_name']
sorted_filtred_results$biotype = subu[rownames(sorted_filtred_results), 'gene_biotype']
sorted_filtred_results$description = subu[rownames(sorted_filtred_results), 'description']

write.table(sorted_filtred_results, file=paste(results_dir, '/effect_right_vs_left_reduced.for.tsv', sep=''), quote=FALSE, sep='\t')
write.table(anti_results, file=paste(results_dir, '/effect_right_vs_left_reduced.anti.tsv', sep=''), quote=FALSE, sep='\t')
