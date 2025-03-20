#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("org.Dm.eg.db")

#
# 0. load libraries
#
library(crayon)
library(clusterProfiler)
library(enrichplot)
library(tictoc)
library(viridis)
library(ggplot2)

#
setwd('/Users/adrian/research/bmcbf/016.saudarkrokur/results/deseq2/h1_sets/')

#
# 2. read files and generate lists of genes
#
filename = '146.tsv' # probably what we rescue
df = read.csv(filename, sep='\t', header=TRUE)
list_146 = df$entrezgene_id[!is.na(df$entrezgene_id)]
list_146

filename = '142.tsv' 
df = read.csv(filename, sep='\t', header=TRUE)
list_142 = df$entrezgene_id[!is.na(df$entrezgene_id)]
list_142

filename = '183.tsv' 
df = read.csv(filename, sep='\t', header=TRUE)
list_183 = df$entrezgene_id[!is.na(df$entrezgene_id)]
list_183

filename = '423.tsv' 
df = read.csv(filename, sep='\t', header=TRUE)
list_423 = df$entrezgene_id[!is.na(df$entrezgene_id)]
list_423

filename = '91.tsv' 
df = read.csv(filename, sep='\t', header=TRUE)
list_91 = df$entrezgene_id[!is.na(df$entrezgene_id)]
list_91

filename = '387.tsv' # probably what we don't rescue
df = read.csv(filename, sep='\t', header=TRUE)
list_387 = df$entrezgene_id[!is.na(df$entrezgene_id)]
list_387

filename = '55.tsv' 
df = read.csv(filename, sep='\t', header=TRUE)
list_55 = df$entrezgene_id[!is.na(df$entrezgene_id)]
list_55

geneLists = list('rescue'=list_146, 'no rescue'=list_387, '142'=list_142, '183'=list_183, 'particular of WT_KO'=list_423, '91'=list_91, '55'=list_55)

#
# 3. run the analysis on different Ontologies
#
# this step takes surprisingly long time. It took xx in an M1 chip
#ck = compareCluster(geneLists, fun="enrichGO", pvalueCutoff=0.05, OrgDb='org.Dm.eg.db')
ck = compareCluster(geneLists, fun="enrichPathway", pvalueCutoff=0.05, organism='fly')

p1 = dotplot(ck, size='count', showCategory=10, font.size=12) + scale_size_area(max_size=9)
print(p1)

# and I have a preference for cividis, but this is just personal preference
my_log_breaks = seq(from=round(log10(0.05)), to=round(log10(min(ck@compareClusterResult$p.adjust))), by=-4)
my_breaks = 10**my_log_breaks
p5 = p1 +  scale_fill_viridis(direction=-1, trans="log", breaks=my_breaks, option='cividis')
print(p5)

# importantly, store your fuctional enrichment in a form of table which will be a supplementary file of your paper
storage_file = 'clusterProfiler_enrichments.h1.tsv'
write.table(ck@compareClusterResult, storage_file, quote=FALSE, sep='\t')

#ggsave('/Users/adrian/scratch/h1.enrichment.svg')
#dev.off()
