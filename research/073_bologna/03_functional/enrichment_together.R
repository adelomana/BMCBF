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
library(ReactomePA)
#
setwd('/Users/adrian/research/bmcbf/073_bologna/results/degs/')
#
# 2. read files and generate lists of genes (no up/down split — full gene set per condition)
#
files <- c(
  'effect_timepoint_D7_vs_D0.responders.tsv',
  'effect_timepoint_D28_vs_D7.responders.tsv',
  'effect_timepoint_D21_vs_D7.responders.tsv',
  'effect_timepoint_D28_vs_D21.responders.tsv',
  'effect_timepoint_D35_vs_D28.responders.tsv'
)

get_gene_lists <- function(filename) {
  
  # extract the label between "effect_timepoint_" and ".responders.tsv"
  label <- sub("^effect_timepoint_(.*)\\.responders\\.tsv$", "\\1", basename(filename))
  
  df <- read.csv(filename, sep = '\t', header = TRUE, row.names = 1)
  print(dim(df))
  
  ids <- bitr(row.names(df), fromType = 'ENSEMBL', 
              toType = 'ENTREZID', OrgDb = 'org.Hs.eg.db')$ENTREZID
  
  message(sprintf("[%s] post-conversion: %d/%d genes mapped (%.1f%%)",
                  label, length(ids), nrow(df), 100 * length(ids) / nrow(df)))
  
  result <- list(ids)
  names(result) <- label
  
  return(result)
}

all_lists <- unlist(lapply(files, get_gene_lists), recursive = FALSE)
names(all_lists)
# D7_vs_D0
# D28_vs_D7
# D21_vs_D7
# D28_vs_D21
# D35_vs_D28

geneLists <- all_lists
names(geneLists) <- sub("_vs_", " vs ", names(all_lists))
#
# 3. run the analysis on different Ontologies
#
ck = compareCluster(geneLists, fun="enrichPathway", pvalueCutoff=0.05, organism='human')
p1 = dotplot(ck, size='count', showCategory=25, font.size=6) + scale_size_area(max_size=9)
print(p1)
# and I have a preference for cividis, but this is just personal preference
my_log_breaks = seq(from=round(log10(0.05)), to=round(log10(min(ck@compareClusterResult$p.adjust))), by=-4)
my_breaks = 10**my_log_breaks
p5 = p1 +  scale_fill_viridis(direction=-1, trans="log", breaks=my_breaks, option='cividis')
print(p5)
# importantly, store your fuctional enrichment in a form of table which will be a supplementary file of your paper
storage_file = 'clusterProfiler_enrichments.together.tsv'
write.table(ck@compareClusterResult, storage_file, quote=FALSE, sep='\t')
#ggsave('/Users/adrian/scratch/h2.enrichment.svg')
#dev.off()