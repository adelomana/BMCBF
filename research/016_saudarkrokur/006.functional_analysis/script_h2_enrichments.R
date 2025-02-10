if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("org.Dm.eg.db")

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
setwd('/Users/adrian/research/016.saudarkrokur/results/deseq2/h2_sets/')

#
# 2. read files and generate lists of genes
#
filename = 'onlypresence.tsv'
df = read.csv(filename, sep='\t', header=TRUE)
list_only_presence = df$entrezgene_id
length(list_only_presence)

filename = 'onlyhigh.tsv'
df = read.csv(filename, sep='\t', header=TRUE)
list_only_high = df$entrezgene_id
length(list_only_high)

geneLists = list('presence'=list_only_presence, 
                'high levels'=list_only_high)

#
# 3. run the analysis on different Ontologies
#
# this step takes surprisingly long time. It took xx in an M1 chip
tic()
#ck = compareCluster(geneLists, fun="enrichGO", pvalueCutoff=0.05, OrgDb='org.Dm.eg.db')
ck = compareCluster(geneLists, fun="enrichPathway", pvalueCutoff=0.05, organism='fly')
toc()

p1 = dotplot(ck, size='count', showCategory=5, font.size=8) 
print(p1)

# and I have a preference for cividis, but this is just personal preference
my_log_breaks = seq(from=round(log10(0.05)), to=round(log10(min(ck@compareClusterResult$p.adjust))), by=-3)
my_breaks = 10**my_log_breaks
p5 = p1 +  scale_fill_viridis(direction=-1, trans="log", breaks=my_breaks, option='cividis')
print(p5)

# importantly, store your fuctional enrichment in a form of table which will be a supplementary file of your paper
storage_file = 'clusterProfiler_enrichments.tsv'
write.table(ck@compareClusterResult, storage_file, quote=FALSE, sep='\t')
