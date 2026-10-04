rm(list = ls())

#
# -1. install libraries
# 
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("DESeq2")

#
# 0. load libraries s
#
library(DESeq2)
library(tximport)
library(biomaRt)
library(BiocParallel)
library(crayon) 
library(ggplot2)
library(stringr)
library(ramify) # this is for clip, for the volcano
library(dplyr)   # if not installed: install.packages("dplyr")


#
# 0. user-defined variables
#
setwd("~/scratch/")
kallisto_dir = "/Users/adrian/research/bmcbf/025_isafjordur/results/000_quantification"
results_dir = '/Users/adrian/research/bmcbf/025_isafjordur/results/002_DEGs/261002'
a_cases = c('SkMel28-MITFKO_ev_skmel28_rep1', 'SkMel28-MITFKO_ev_skmel28_rep2', 'SkMel28-MITFKO_ev_skmel28_rep3', 'SkMel28-MITFKO_ev_skmel28_rep4', 'SkMel28-MITFKO_mitf_x6_rep1', 'SkMel28-MITFKO_mitf_x6_rep3', 'SkMel28-MITFKO_mitf_x6_rep4')

#
# 1. generate gene to transcript mapping
#
df = read.csv('/Users/adrian/software/kallisto/human_index_standard/annotation.tsv', sep='\t')
t2g = df[, 2:3]
dim(t2g)

annotation_file = '/Users/adrian/software/kallisto/human_index_standard/annotation.tsv'
full_annotation = read.csv(annotation_file, sep='\t')
annotation <- full_annotation %>% distinct(ensembl_gene_id, .keep_all = TRUE)

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
threshold = 10
effect_size_threshold = log2(2)

#
# 3.1. contrast 
#
txi = tximport(metadata$path, type="kallisto", tx2gene=t2g, ignoreTxVersion = TRUE)
dds = DESeqDataSetFromTximport(txi, colData=metadata, design=~genotype) 
dds$genotype = relevel(dds$genotype, ref="wt")

# keep features with at least 10 counts median in one sample
cts <- counts(dds)
cat(blue(paste('size before counts filtering:', dim(dds)[1], sep=' ')), fill=TRUE)
med_g1 <- matrixStats::rowMedians(cts[, 1:4])
med_g2 <- matrixStats::rowMedians(cts[, 5:7])
keep <- pmax(med_g1, med_g2) >= threshold   # threshold = 10
dds  <- dds[keep, ]
cat(blue(paste('size after counts filtering:', dim(dds)[1], sep=' ')), fill=TRUE)

dds = DESeq(dds, test="LRT", reduced=~1)
resultsNames(dds)

res = results(dds, name = "genotype_ko_vs_wt", alpha = 0.05)
summary(res)

# 3.5 Shrink LFC for stable effect-size reporting
res_shr <- lfcShrink(dds, coef = "genotype_ko_vs_wt", type = "apeglm", res = res)
plotMA(res_shr, ylim = c(-5, 5), main=paste('shrunk egfl7'))

# 3.6. Build a final table with:
#    - padj from res0 (correct for the H0: LFC = 0 test)
#    - shrunken LFC from res0_shr (recommended for reporting/thresholding)
full_results <- data.frame(
  gene_id = rownames(res),
  baseMean = res$baseMean,
  log2FC_MLE = res$log2FoldChange,
  lfcSE = res$lfcSE,
  stat = res$stat,
  pvalue = res$pvalue,
  padj = res$padj,
  log2FC_shr = res_shr$log2FoldChange,
  stringsAsFactors = FALSE
)

# add counts differences 
design_name = 'genotype'
level_A = 'ko'         
level_B <- 'wt'
print('levels')
print(level_A)
print(level_B)

norm_counts <- counts(dds, normalized = TRUE)
idx_A <- which(colData(dds)[[design_name]] == level_A)
idx_B <- which(colData(dds)[[design_name]] == level_B)
print('indexes')
print(idx_A)
print(idx_B)

median_A <- rowMedians(norm_counts[, idx_A, drop = FALSE])
median_B <- rowMedians(norm_counts[, idx_B, drop = FALSE])
delta_counts <- median_B - median_A
names(delta_counts) <- rownames(norm_counts)
full_results$delta_counts <- delta_counts[match(full_results$gene_id, names(delta_counts))]

tpm_mat <- txi$abundance
median_TPM_A <- matrixStats::rowMedians(tpm_mat[, idx_A, drop = FALSE])
median_TPM_B <- matrixStats::rowMedians(tpm_mat[, idx_B, drop = FALSE])
names(median_TPM_A) <- rownames(tpm_mat)
names(median_TPM_B) <- rownames(tpm_mat)
full_results$median_TPM_A <- median_TPM_A[full_results$gene_id]
full_results$median_TPM_B <- median_TPM_B[full_results$gene_id]

# add gene names and descriptions
tempo = sapply(strsplit(rownames(full_results), split='.',fixed=TRUE), function(x) (x[1]))
rownames(full_results) = tempo
length(tempo)
sub = annotation[annotation$ensembl_gene_id %in% tempo, ]
dim(sub)
sub = sub[, c(3, 4, 5, 6)]
sub$description2 = sapply(strsplit(sub$description, split='[Source',fixed=TRUE), function(x) (x[1]))
rownames(sub) <- sub$ensembl_gene_id

full_results$ensembl_id = sub[rownames(full_results), 'ensembl_gene_id']
full_results$gene_name = sub[rownames(full_results), 'external_gene_name']
full_results$biotype = sub[rownames(full_results), 'gene_biotype']
full_results$description2 = sub[rownames(full_results), 'description2']


# 6.4 Storing full and final “biologically relevant” subset (you define thresholds)
responders <- subset(
  full_results,
  !is.na(padj) &
    padj < 0.05 &
    abs(log2FC_shr) >= 1 &
    abs(delta_counts) >= 50
)
responders <- responders[order(responders$padj), ]

# no response
no_responders <- full_results[!full_results$gene_id %in% responders$gene_id, ]

# report values
dim(full_results)
dim(no_responders)
sum(full_results$padj < 0.05, na.rm = TRUE)
dim(responders)

# store results
filename = paste(results_dir, '/effect_ko_vs_wt', '.full.tsv', sep='')
write.table(full_results, file=filename, quote=FALSE, sep='\t')

filename = paste(results_dir, '/effect_ko_vs_wt', '.responders.tsv', sep='')
write.table(responders, file=filename, quote=FALSE, sep='\t')

design_name = 'genotype'
sample_flag = 'ko'
control_flag = 'wt'

#
# 7. visualization
#

# 7.1. a simple PCA
plotPCA(rlog(dds), intgroup=c(design_name)) + ggtitle(paste(design_name, sample_flag, control_flag))

# 7.2. a simple volcano plot
plotting_x = responders$log2FC_shr
y = responders$padj
epsilon = min(y[y !=0])
plotting_y = -log10(y + epsilon) ## why???
df = data.frame(plotting_x=plotting_x, plotting_y=plotting_y)
reds = df[df$plotting_x > 0, ]
blues = df[df$plotting_x < 0, ]

plotting_x = no_responders$log2FC_shr
plotting_y = -log10(no_responders$padj)
blacks = data.frame(plotting_x=plotting_x, plotting_y=plotting_y)

p <- ggplot() + 
  geom_point(data=reds, aes(plotting_x, plotting_y), color = "red", size=1, shape=19, alpha=0.5, stroke=0) + 
  geom_point(data=blues, aes(plotting_x, plotting_y), color = "blue", size=1, shape=19, alpha=0.5, stroke=0) +
  geom_point(data=blacks, aes(plotting_x, plotting_y), color = "black", size=1, shape=19, alpha=0.1, stroke=0) +
  labs(x='log2FC', y='-log10 adj P') + 
  theme_linedraw() 
filename = paste(results_dir, '/', design_name, '_', sample_flag, '_vs_', control_flag, '.simple_volcano.pdf', sep='')
ggsave(filename, plot = p)

#
# 5.4. a rather elaborated volcano including TPM values
#
plotting_x = responders$log2FC_shr
y = responders$padj
epsilon = min(y[y !=0])
plotting_y = -log10(y + epsilon) 
plotting_z = log10(rowMeans(responders[, c('median_TPM_A', 'median_TPM_B')]))

df = data.frame(plotting_x=clip(plotting_x, .min=-3, .max=3), 
                plotting_y=clip(plotting_y, .min=0, .max=38),
                plotting_z=clip(plotting_z, .min=0, .max=3))

reds = df[df$plotting_x > 0, ]
blues = df[df$plotting_x < 0, ]

plotting_x = no_responders$log2FC_shr
plotting_y = -log10(no_responders$padj)
blacks = data.frame(plotting_x=clip(plotting_x, .min=-3, .max=3), plotting_y=clip(plotting_y, .min=0, .max=38))

p <- ggplot() +  
  geom_point(data=reds, aes(x=plotting_x, y=plotting_y, color=plotting_z), , size=3, shape=19, alpha=2/3, stroke=0) + 
  geom_point(data=blues, aes(plotting_x, plotting_y, color=plotting_z), size=3, shape=19, alpha=2/3, stroke=0) +
  geom_point(data=blacks, aes(plotting_x, plotting_y), color = "black", size=1, shape=19, alpha=1/3, stroke=0) +
  labs(x=expression('Expression difference [log'[2]~'FC]'), y=expression('Significance [log'[10]~'adjusted P]'), color=expression('Expression average [log'[10]~'TPM]'), title=paste(design_name, sample_flag, control_flag)) + 
  theme_linedraw() +
  geom_segment(aes(x=-1, xend=-1, y=-log10(0.05), yend=50), linetype=2) +
  geom_segment(aes(x=1, xend=1, y=-log10(0.05), yend=50), linetype=2) +
  geom_segment(aes(x=-8, xend=-1, y=-log10(0.05), yend=-log10(0.05)), linetype=2) +
  geom_segment(aes(x=1, xend=8, y=-log10(0.05), yend=-log10(0.05)), linetype=2) +
  xlim(-3, 3) + 
  ylim(-1, 38) +
  scale_color_viridis_c(option = "cividis") 
filename = paste(results_dir, '/', design_name, '_', sample_flag, '_vs_', control_flag, '.elaborated_volcano.pdf', sep='')
ggsave(filename, plot = p)



