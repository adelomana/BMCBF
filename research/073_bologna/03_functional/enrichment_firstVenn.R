#if (!require("BiocManager", quietly = TRUE))
# install.packages("BiocManager")
#BiocManager::install("org.Dm.eg.db")
#BiocManager::install("clusterProfiler")

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
setwd('/Users/adrian/research/bmcbf/073_bologna/results/degs/')

#
# 2. read files and generate lists of genes
#
df <- read.csv("first_sets.csv")          


# convert
l = bitr(df$left, fromType = 'ENSEMBL', toType = 'ENTREZID', OrgDb = 'org.Hs.eg.db')$ENTREZID
m = bitr(df$middle, fromType = 'ENSEMBL', toType = 'ENTREZID', OrgDb = 'org.Hs.eg.db')$ENTREZID
r = bitr(df$right, fromType = 'ENSEMBL', toType = 'ENTREZID', OrgDb = 'org.Hs.eg.db')$ENTREZID

length(l)
length(m)
length(r)

geneLists <- list(leftList = l, middleList = m, rightList = r)

ck = compareCluster(geneLists, fun="enrichPathway", pvalueCutoff=0.05, organism='human')

p1 = dotplot(ck, size='count', showCategory=10, font.size=6) + scale_size_area(max_size=9)
print(p1)

my_log_breaks = seq(from=round(log10(0.05)), to=round(log10(min(ck@compareClusterResult$p.adjust))), by=-4)
my_breaks = 10**my_log_breaks
p5 = p1 +  scale_fill_viridis(direction=-1, trans="log", breaks=my_breaks, option='cividis')
print(p5)

storage_file = '/Users/adrian/research/bmcbf/073_bologna/results/enrichment/clusterProfiler_enrichments.firstVenn.tsv'
write.table(ck@compareClusterResult, storage_file, quote=FALSE, sep='\t')

#ggsave('/Users/adrian/scratch/enrichment.svg')
#dev.off()

