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
setwd('/Users/adrian/research/016.saudarkrokur/results/deseq2/')

#
# 2. read files and generate lists of genes
#
filename = 'effect_ko_wt.tsv'
df = read.csv(filename, sep='\t', header=TRUE)
df_up = df[df$log2FoldChange > 0, ] 
df_down = df[df$log2FoldChange < 0, ] 
print(dim(df))

ensemblIDs = row.names(df_up)
convertedIDs = bitr(ensemblIDs, fromType = 'ENSEMBL', toType = 'ENTREZID', OrgDb='org.Dm.eg.db')
list_one_up = convertedIDs$ENTREZID
length(list_one_up)

ensemblIDs = row.names(df_down)
convertedIDs = bitr(ensemblIDs, fromType = 'ENSEMBL', toType = 'ENTREZID', OrgDb='org.Dm.eg.db')
list_one_down = convertedIDs$ENTREZID
length(list_one_down)

filename = 'effect_h1M8_ko.tsv'
df = read.csv(filename, sep='\t', header=TRUE)
df_up = df[df$log2FoldChange > 0, ] 
df_down = df[df$log2FoldChange < 0, ] 
print(dim(df))

ensemblIDs = row.names(df_up)
convertedIDs = bitr(ensemblIDs, fromType = 'ENSEMBL', toType = 'ENTREZID', OrgDb='org.Dm.eg.db')
list_two_up = convertedIDs$ENTREZID
length(list_two_up)

ensemblIDs = row.names(df_down)
convertedIDs = bitr(ensemblIDs, fromType = 'ENSEMBL', toType = 'ENTREZID', OrgDb='org.Dm.eg.db')
list_two_down = convertedIDs$ENTREZID
length(list_two_down)

filename = 'effect_h2F14_wt.tsv'
df = read.csv(filename, sep='\t', header=TRUE)
df_up = df[df$log2FoldChange > 0, ] 
df_down = df[df$log2FoldChange < 0, ] 
print(dim(df))

ensemblIDs = row.names(df_up)
convertedIDs = bitr(ensemblIDs, fromType = 'ENSEMBL', toType = 'ENTREZID', OrgDb='org.Dm.eg.db')
list_three_up = convertedIDs$ENTREZID
length(list_three_up)

ensemblIDs = row.names(df_down)
convertedIDs = bitr(ensemblIDs, fromType = 'ENSEMBL', toType = 'ENTREZID', OrgDb='org.Dm.eg.db')
list_three_down = convertedIDs$ENTREZID
length(list_three_down)

filename = 'effect_h2M8_wt.tsv'
df = read.csv(filename, sep='\t', header=TRUE)
df_up = df[df$log2FoldChange > 0, ] 
df_down = df[df$log2FoldChange < 0, ] 
print(dim(df))

ensemblIDs = row.names(df_up)
convertedIDs = bitr(ensemblIDs, fromType = 'ENSEMBL', toType = 'ENTREZID', OrgDb='org.Dm.eg.db')
list_four_up = convertedIDs$ENTREZID
length(list_four_up)

ensemblIDs = row.names(df_down)
convertedIDs = bitr(ensemblIDs, fromType = 'ENSEMBL', toType = 'ENTREZID', OrgDb='org.Dm.eg.db')
list_four_down = convertedIDs$ENTREZID
length(list_four_down)

geneLists = list('KO vs WT up'=list_one_up, 
                'KO vs WT down'=list_one_down,
                'h1M8 vs KO up'=list_two_up, 
                'h1M8 vs KO down'=list_two_down,
                'h2F14 vs WT up'=list_three_up, 
                'h2F14 vs WT down'=list_three_down,
                'h2M8 vs WT up'=list_four_up, 
                'h2M8 vs WT down'=list_four_down)


#
# 3. run the analysis on different Ontologies
#
# this step takes surprisingly long time. It took xx in an M1 chip
tic()
#ck = compareCluster(geneLists, fun="enrichGO", pvalueCutoff=0.05, OrgDb='org.Dm.eg.db')
ck = compareCluster(geneLists, fun="enrichPathway", pvalueCutoff=0.05, organism='fly')
#toc()

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
