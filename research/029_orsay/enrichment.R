rm(list = ls())
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
setwd('/Users/adrian/research/bmcbf/025_isafjordur/results/002_DEGs/')

#
# 2. read files and generate lists of genes
#
filename = 'effect_ko_vs_wt.skmel.for.tsv'
df = read.csv(filename, sep='\t', header=TRUE)
ensemblIDs = row.names(df)
df_up = df[df$log2FoldChange > 0, ] 
df_down = df[df$log2FoldChange < 0, ] 
ensemblIDs = row.names(df_up)
convertedIDs = bitr(ensemblIDs, fromType = 'ENSEMBL', toType = 'ENTREZID', OrgDb='org.Hs.eg.db')
list_one_up = convertedIDs$ENTREZID
ensemblIDs = row.names(df_down)
convertedIDs = bitr(ensemblIDs, fromType = 'ENSEMBL', toType = 'ENTREZID', OrgDb='org.Hs.eg.db')
list_one_down = convertedIDs$ENTREZID

filename = 'effect_ko_vs_wt.u2os.two.for.tsv'
df = read.csv(filename, sep='\t', header=TRUE)
ensemblIDs = row.names(df)
df_up = df[df$log2FoldChange > 0, ] 
df_down = df[df$log2FoldChange < 0, ] 
print(dim(df_up))
print(dim(df_down))
ensemblIDs = row.names(df_up)
convertedIDs = bitr(ensemblIDs, fromType = 'ENSEMBL', toType = 'ENTREZID', OrgDb='org.Hs.eg.db')
list_two_up = convertedIDs$ENTREZID
ensemblIDs = row.names(df_down)
convertedIDs = bitr(ensemblIDs, fromType = 'ENSEMBL', toType = 'ENTREZID', OrgDb='org.Hs.eg.db')
list_two_down = convertedIDs$ENTREZID

filename = 'effect_ko_vs_wt.u2os.group.for.tsv'
df = read.csv(filename, sep='\t', header=TRUE)
ensemblIDs = row.names(df)
df_up = df[df$log2FoldChange > 0, ] 
df_down = df[df$log2FoldChange < 0, ] 
ensemblIDs = row.names(df_up)
convertedIDs = bitr(ensemblIDs, fromType = 'ENSEMBL', toType = 'ENTREZID', OrgDb='org.Hs.eg.db')
list_three_up = convertedIDs$ENTREZID
ensemblIDs = row.names(df_down)
convertedIDs = bitr(ensemblIDs, fromType = 'ENSEMBL', toType = 'ENTREZID', OrgDb='org.Hs.eg.db')
list_three_down = convertedIDs$ENTREZID

geneLists = list('skmel up'=list_one_up, 
                 'skmel down'=list_one_down,
                 'u2os two up'=list_two_up, 
                 'u2os two down'=list_two_down,
                 'u2os group up'=list_three_up, 
                 'u2os group down'=list_three_down)

#
# 3. run the analysis on different Ontologies
#
# this step takes surprisingly long time. It took xx in an M1 chip
ck = compareCluster(geneLists, fun="enrichPathway", pvalueCutoff=0.05)

p1 = dotplot(ck, size='count', showCategory=10, font.size=6) + scale_size_area(max_size=9)
print(p1)
my_log_breaks = seq(from=round(log10(0.05)), to=round(log10(min(ck@compareClusterResult$p.adjust))), by=-4)
my_breaks = 10**my_log_breaks
p5 = p1 +  scale_fill_viridis(direction=-1, trans="log", breaks=my_breaks, option='cividis')
print(p5)

# importantly, store your fuctional enrichment in a form of table which will be a supplementary file of your paper
storage_file = 'clusterProfiler_enrichments.tsv'

all_names = c()
for (index in 1:dim(ck@compareClusterResult)[1]) {
  v = strsplit(ck@compareClusterResult[index, 12], '/')[[1]]
  results = bitr(v, fromType = 'ENTREZID', toType = 'SYMBOL', OrgDb='org.Hs.eg.db')
  names = results$SYMBOL
  sorted_names = sort(names)
  joined_names = paste(sorted_names, collapse=', ')
  print(joined_names)
  all_names = c(all_names, joined_names)
}

ck@compareClusterResult$transformed_ids = all_names
write.table(ck@compareClusterResult, storage_file, quote=FALSE, sep='\t')

#ggsave('/Users/adrian/scratch/h1.enrichment.svg')
#dev.off()
