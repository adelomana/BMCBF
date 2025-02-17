#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("clusterProfiler")
#BiocManager::install("org.Hs.eg.db")
#BiocManager::install("ReactomePA")
#BiocManager::install("tictoc")
#BiocManager::install("wesanderson")

#
# 0. load libraries
#
library(crayon)
library(clusterProfiler)
library(enrichplot)
library(tictoc)
library(viridis)
library(ggplot2)
library(wesanderson)

#
# 2. read files and generate lists of genes
#
filename = 'KO_VS_17-1.8.tsv'
df = read.table(filename, sep='\t', header=TRUE)
df_up = df[df$log2FoldChange > 0, ] 
df_down = df[df$log2FoldChange < 0, ] 

ensemblIDs = row.names(df_up)
convertedIDs = bitr(ensemblIDs, fromType = 'ENSEMBL', toType = 'ENTREZID', OrgDb='org.Hs.eg.db')
list_one_up = convertedIDs$ENTREZID
length(list_one_up)

ensemblIDs = row.names(df_down)
convertedIDs = bitr(ensemblIDs, fromType = 'ENSEMBL', toType = 'ENTREZID', OrgDb='org.Hs.eg.db')
list_one_down = convertedIDs$ENTREZID
length(list_one_down)

filename = 'KO_VS_17-4.2.tsv'
df = read.table(filename, sep='\t', header=TRUE)
df_up = df[df$log2FoldChange > 0, ] 
df_down = df[df$log2FoldChange < 0, ] 

ensemblIDs = row.names(df_up)
convertedIDs = bitr(ensemblIDs, fromType = 'ENSEMBL', toType = 'ENTREZID', OrgDb='org.Hs.eg.db')
list_two_up = convertedIDs$ENTREZID
length(list_two_up)

ensemblIDs = row.names(df_down)
convertedIDs = bitr(ensemblIDs, fromType = 'ENSEMBL', toType = 'ENTREZID', OrgDb='org.Hs.eg.db')
list_two_down = convertedIDs$ENTREZID
length(list_two_down)

geneLists = list('case A up'=list_one_up, 
                'case A down'=list_one_down,
                'case B up'=list_two_up, 
                'case B down'=list_two_down)

#
# 3. run the analysis on different Ontologies
#
# this step takes surprisingly long time
tic()
ck = compareCluster(geneLists, 
                    fun="enrichPathway", 
                    pvalueCutoff=0.05)
toc()
# this long step takes 42 seconds on an M1 chip
p1 = dotplot(ck, size='count', showCategory=5, font.size=8) 
print(p1)

# please don't use this
# https://github.com/karthik/wesanderson
pal <- wes_palette(name="Moonrise3", n=5, type="continuous")
p2 <- p1 + scale_fill_gradientn(colours = pal)
print(p2)

# please use viridis and a more intuitive direction for P values, where lower is better 
p3 = p1 +  scale_fill_viridis(direction=-1)
print(p3)

# but log scale is more appropriate
my_log_breaks = round(log10(0.05)):round(log10(min(ck@compareClusterResult$p.adjust)))
my_breaks = 10**my_log_breaks
p4 = p1 +  scale_fill_viridis(direction=-1, trans="log", breaks = my_breaks)
print(p4)

# and I have a preference for cividis, but this is just personal preference
p5 = p1 +  scale_fill_viridis(direction=-1, trans="log", breaks = my_breaks, option='cividis')
print(p5)
# arguably this plot communicates best data patterns, IMHO

# importantly, store your fuctional enrichment in a form of table which will be a supplementary file of your paper
storage_file = 'clusterProfiler_enrichments_RP.tsv'
write.table(ck@compareClusterResult, storage_file, quote=FALSE, sep='\t')
