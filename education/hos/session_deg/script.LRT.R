rm(list = ls())

#
# -1. packages installation
#

# use the following block of code if libraries are not installed in your computer

# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# 
# BiocManager::install("tximport")
# BiocManager::install("DESeq2")
# BiocManager::install('rhdf5')
# BiocManager::install('this.path')
# BiocManager::install('ramify')
# BiocManager::install('crayon')
# BiocManager::install('apeglm')

library(tximport)       # required to read input files
library(DESeq2)         # the library that will call DEGs
library(crayon)         # so the messages are blue
library(this.path)      # necessary to locate where this file is
library(ggplot2)        # useful for plotting
library(ramify)         # necessary for the clip function
library(rhdf5)          # necessary for reading the input files

#
# 0. user-defined variables
#

# set your working directory. This is an option, but you can change as you prefer.
# please be familiar with getwd() and setwd(), very useful commands to define your working directory, which is critical to know where the outputs will be
script_path = this.dir()
script_path
setwd(script_path) 

kallisto_dir = "kallisto_output"
results_dir = 'DEGs_DESeq2'

# 
# 1. get todays working data: kallisto output from two conditions
#
#system('wget https://ireigogn.hi.is/index.php/s/Lj9APBfmgXiTpWY/download/kallisto_output.tgz')
#untar('kallisto_output.tgz')

list.files('kallisto_output')
list.files('kallisto_output/WT_with_IFN_1/')
df = read.csv('kallisto_output/WT_with_IFN_1/abundance.tsv', sep='\t')
View(df)

#
# 2. get annotation mapping from transcript to genes
#
#system('wget https://github.com/pachterlab/kallisto-transcriptome-indices/releases/download/v1/human_index_standard.tar.xz')
#untar('human_index_standard.tar.xz')

t2g = read.csv('t2g.txt', sep='\t', header = FALSE)
View(t2g)

#
# 3. define metadata
#
dirnames = list.dirs(kallisto_dir, full.names=TRUE, recursive=FALSE)
paths = file.path(dirnames, 'abundance.h5')
labels = sapply(strsplit(paths, split='/',fixed=TRUE), function(x) (x[2]))
replicates = rep(c('A', 'B', 'C'), 2)
treatments = rep(c(rep('with', 3), rep('without', 3)))

metadata = data.frame(labels)
metadata$replicate = replicates
metadata$treatment = treatments
metadata$path = paths

View(metadata)

#
# 4. bring the expression profiles into DESeq2
#
txi = tximport(metadata$path, type="kallisto", tx2gene=t2g)

dds = DESeqDataSetFromTximport(txi, colData=metadata, design=~treatment) 
dds$treatment = relevel(dds$treatment, ref="without")
# minimal filter on undetected genes: at least 10 reads in three samples.
# additional filters should happen a posteriori, because removing a lot of genes would affect DESeq2 assumptions
cat(blue(paste('size before counts filtering:', dim(dds)[1], sep=' ')), fill=TRUE)
smallestGroupSize <- 3
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep,]
cat(blue(paste('size after counts filtering:', dim(dds)[1], sep=' ')), fill=TRUE)

#
# 6. Statistical test
#

# Fit model (Wald test by default)
dds <- DESeq(dds, test = "LRT", reduced = ~ 1)

# 6.1.get p-values / padj for H0: log2FC = 0
resultsNames(dds)
res <- results(dds, alpha = 0.05)

# 6.2 Shrink LFC for stable effect-size reporting
res_shr <- lfcShrink(dds, coef = "treatment_with_vs_without", type = "apeglm")
# visualize shrinkage
plotMA(res, ylim = c(-5, 5))
plotMA(res_shr, ylim = c(-5, 5))

# 6.3 Build a final table with:
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
condition_col <- "treatment"     # change if your metadata column differs
level_A <- "with"                # change to your level name
level_B <- "without"             # change to your level name
stopifnot(condition_col %in% colnames(colData(dds)))
stopifnot(all(c(level_A, level_B) %in% as.character(unique(colData(dds)[[condition_col]]))))
norm_counts <- counts(dds, normalized = TRUE)
idx_A <- which(colData(dds)[[condition_col]] == level_A)
idx_B <- which(colData(dds)[[condition_col]] == level_B)
median_A <- rowMedians(norm_counts[, idx_A, drop = FALSE])
median_B <- rowMedians(norm_counts[, idx_B, drop = FALSE])
delta_counts <- median_A - median_B
names(delta_counts) <- rownames(norm_counts)
full_results$delta_counts <- delta_counts[match(full_results$gene_id, names(delta_counts))]

tpm_mat <- txi$abundance
median_TPM_A <- matrixStats::rowMedians(tpm_mat[, idx_A, drop = FALSE])
median_TPM_B <- matrixStats::rowMedians(tpm_mat[, idx_B, drop = FALSE])
names(median_TPM_A) <- rownames(tpm_mat)
names(median_TPM_B) <- rownames(tpm_mat)
full_results$median_TPM_A <- median_TPM_A[full_results$gene_id]
full_results$median_TPM_B <- median_TPM_B[full_results$gene_id]

# 6.4 Example: final “biologically relevant” subset (you define thresholds)
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
dir.create(results_dir)
write.table(full_results, file=paste(results_dir, '/effect_IFN_vs_noIFN.tsv', sep=''), quote=FALSE, sep='\t')

#
# 5. visualization
#

# 5.1. a simple PCA
plotPCA(rlog(dds), intgroup=c('treatment')) + ggtitle('effect IFN vs control')

# 5.2. a simple volcano plot
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

ggplot() + 
  geom_point(data=reds, aes(plotting_x, plotting_y), color = "red", size=1, shape=19, alpha=0.5, stroke=0) + 
  geom_point(data=blues, aes(plotting_x, plotting_y), color = "blue", size=1, shape=19, alpha=0.5, stroke=0) +
  geom_point(data=blacks, aes(plotting_x, plotting_y), color = "black", size=1, shape=19, alpha=0.1, stroke=0) +
  labs(x='log2FC', y='-log10 adj P') + 
  theme_linedraw() 
       
#
# 5.4. a rather elaborated volcano including TPM values
#
plotting_x = responders$log2FC_shr
y = responders$padj
epsilon = min(y[y !=0])
plotting_y = -log10(y + epsilon) 
plotting_z = log10(rowMeans(responders[, c('median_TPM_A', 'median_TPM_B')]))

df = data.frame(plotting_x=clip(plotting_x, .min=-8, .max=8), 
                plotting_y=clip(plotting_y, .min=0, .max=50),
                plotting_z=clip(plotting_z, .min=0, .max=3))

reds = df[df$plotting_x > 0, ]
blues = df[df$plotting_x < 0, ]

plotting_x = no_responders$log2FC_shr
plotting_y = -log10(no_responders$padj)
blacks = data.frame(plotting_x=clip(plotting_x, .min=-8, .max=8), plotting_y=clip(plotting_y, .min=0, .max=50))

ggplot() + 
  geom_point(data=reds, aes(x=plotting_x, y=plotting_y, color=plotting_z), , size=3, shape=19, alpha=2/3, stroke=0) + 
  geom_point(data=blues, aes(plotting_x, plotting_y, color=plotting_z), size=3, shape=19, alpha=2/3, stroke=0) +
  geom_point(data=blacks, aes(plotting_x, plotting_y), color = "black", size=1, shape=19, alpha=1/3, stroke=0) +
  labs(x=expression('Expression difference [log'[2]~'FC]'), y=expression('Significance [log'[10]~'adjusted P]'), color=expression('Expression average [log'[10]~'TPM]')) + 
  theme_linedraw() +
  geom_segment(aes(x=-1, xend=-1, y=-log10(0.05), yend=50), linetype=2) +
  geom_segment(aes(x=1, xend=1, y=-log10(0.05), yend=50), linetype=2) +
  geom_segment(aes(x=-8, xend=-1, y=-log10(0.05), yend=-log10(0.05)), linetype=2) +
  geom_segment(aes(x=1, xend=8, y=-log10(0.05), yend=-log10(0.05)), linetype=2) +
  xlim(-8, 8) +
  scale_color_viridis_c(option = "cividis") 






