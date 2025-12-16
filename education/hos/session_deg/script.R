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
system('wget https://ireigogn.hi.is/index.php/s/Lj9APBfmgXiTpWY/download/kallisto_output.tgz')
untar('kallisto_output.tgz')

list.files('kallisto_output')
list.files('kallisto_output/WT_with_IFN_1/')
df = read.csv('kallisto_output/WT_with_IFN_1/abundance.tsv', sep='\t')
View(df)

#
# 2. get annotation mapping from transcript to genes
#
system('wget https://github.com/pachterlab/kallisto-transcriptome-indices/releases/download/v1/human_index_standard.tar.xz')
untar('human_index_standard.tar.xz')
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

#
# 5. minimal filter on undetected genes: at least 10 reads in three samples.
# additional filters should happen a posteriori, because removing a lot of genes would affect DESeq2 assumptions
#
cat(blue(paste('size before counts filtering:', dim(dds)[1], sep=' ')), fill=TRUE)
smallestGroupSize <- 3
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep,]
cat(blue(paste('size after counts filtering:', dim(dds)[1], sep=' ')), fill=TRUE)

#
# 6. Statisical test
#
# LRT and Wald tests are available, more info here: https://hbctraining.github.io/DGE_workshop_salmon/lessons/08_DGE_LRT.html
dds = DESeq(dds, test="LRT", reduced=~1)

# retrieve significant results
res = results(dds, parallel=TRUE, alpha=0.05) 
res1 = res[which(res$padj < 0.05), ]

#
# 7. filter DEGs that are very lowly expressed, below an average normalized counts of 50, baseMean >= 50
#
cat(blue(paste('size before filtering:', dim(res1)[1], sep=' ')), fill=TRUE)
norm_counts <- counts(dds, normalized = TRUE)[rownames(res1), ]
keep <- rowMeans(norm_counts) >= 50
res2 = res1[keep, ]
cat(blue(paste('size after filtering:', dim(res2)[1], sep=' ')), fill=TRUE)

#
# 8. filter DEGs that do not show great difference at the count level, even more strict than previous one
#
# filter results that do not have at least 50 counts difference
cat(blue(paste('size before filtering:', dim(res2)[1], sep=' ')), fill=TRUE)
a = counts(dds, normalize=TRUE)[rownames(res2), 1:3]
b = counts(dds, normalize=TRUE)[rownames(res2), 4:6]
c = rowMedians(a) - rowMedians(b)
keep = abs(c) >= 50
res3 = res2[keep, ]
cat(blue(paste('size after filtering:', dim(res3)[1], sep=' ')), fill=TRUE)

#
# 9. filter on log2FC of one for final results
#
cat(blue(paste('size before filtering:', dim(res3)[1], sep=' ')), fill=TRUE)
filtred_results = res3[abs(res3$log2FoldChange) >= 1, ]
sorted_filtred_results = filtred_results[order(filtred_results[["padj"]]),]
cat(blue(paste('size after filtering:', dim(sorted_filtred_results)[1], sep=' ')), fill=TRUE)

# get anti results for plotting
anti_results = res[which(res$padj > 0.05 | abs(res$log2FoldChange) < 1), ]

# store results
dir.create(results_dir)
write.table(sorted_filtred_results, file=paste(results_dir, '/effect_IFN_vs_noIFN.tsv', sep=''), quote=FALSE, sep='\t')
write.table(anti_results, file=paste(results_dir, '/effect_IFN_vs_noIFN.anti.tsv', sep=''), quote=FALSE, sep='\t')

#
# 5. visualization
#

# 5.1. a simple PCA
plotPCA(rlog(dds), intgroup=c('treatment')) + ggtitle('effect IFN vs control')

# 5.2.a volcano plot
plotting_x = sorted_filtred_results$log2FoldChange
y = sorted_filtred_results$padj
epsilon = min(y[y !=0])
plotting_y = -log10(y + epsilon) ## why???
df = data.frame(plotting_x=plotting_x, plotting_y=plotting_y)
reds = df[df$plotting_x > 0, ]
blues = df[df$plotting_x < 0, ]

plotting_x = anti_results$log2FoldChange
plotting_y = -log10(anti_results$padj)
blacks = data.frame(plotting_x=plotting_x, plotting_y=plotting_y)

ggplot() + 
  geom_point(data=reds, aes(plotting_x, plotting_y), color = "red", size=1, shape=19, alpha=0.5, stroke=0) + 
  geom_point(data=blues, aes(plotting_x, plotting_y), color = "blue", size=1, shape=19, alpha=0.5, stroke=0) +
  geom_point(data=blacks, aes(plotting_x, plotting_y), color = "black", size=1, shape=19, alpha=0.1, stroke=0) +
  labs(x='log2FC', y='-log10 adj P') + 
  theme_linedraw() 
       
#               
# 5.3. a rather centered volcano
#
plotting_x = sorted_filtred_results$log2FoldChange
y = sorted_filtred_results$padj
epsilon = min(y[y !=0])
plotting_y = -log10(y + epsilon) 
df = data.frame(plotting_x=clip(plotting_x, .min=-8, .max=8), plotting_y=clip(plotting_y, .min=0, .max=50))
reds = df[df$plotting_x > 0, ]
blues = df[df$plotting_x < 0, ]

plotting_x = anti_results$log2FoldChange
plotting_y = -log10(anti_results$padj)
blacks = data.frame(plotting_x=clip(plotting_x, .min=-8, .max=8), plotting_y=clip(plotting_y, .min=0, .max=50))

ggplot() + 
  geom_point(data=reds, aes(plotting_x, plotting_y), color = "red", size=3, shape=19, alpha=0.5, stroke=0) + 
  geom_point(data=blues, aes(plotting_x, plotting_y), color = "blue", size=3, shape=19, alpha=0.5, stroke=0) +
  geom_point(data=blacks, aes(plotting_x, plotting_y), color = "black", size=1, shape=19, alpha=0.2, stroke=0) +
  labs(x=expression('Expression [log'[2]~'FC]'), y=expression('Significance [log'[10]~'adjusted P]')) + 
  theme_linedraw() +
  geom_segment(aes(x=-1, xend=-1, y=-log10(0.05), yend=50), linetype=2) +
  geom_segment(aes(x=1, xend=1, y=-log10(0.05), yend=50), linetype=2) +
  geom_segment(aes(x=-8, xend=-1, y=-log10(0.05), yend=-log10(0.05)), linetype=2) +
  geom_segment(aes(x=1, xend=8, y=-log10(0.05), yend=-log10(0.05)), linetype=2) +
  xlim(-8.2, 8) 


#
# 5.4. a volcano including TPM values
#
plotting_x = sorted_filtred_results$log2FoldChange

y = sorted_filtred_results$padj
epsilon = min(y[y !=0])
plotting_y = -log10(y + epsilon) 

plotting_z = log10(rowMeans(txi$abundance[rownames(sorted_filtred_results), ]))

df = data.frame(plotting_x=clip(plotting_x, .min=-8, .max=8), 
                plotting_y=clip(plotting_y, .min=0, .max=50),
                plotting_z=clip(plotting_z, .min=0, .max=3))

reds = df[df$plotting_x > 0, ]
blues = df[df$plotting_x < 0, ]

plotting_x = anti_results$log2FoldChange
plotting_y = -log10(anti_results$padj)
blacks = data.frame(plotting_x=clip(plotting_x, .min=-8, .max=8), plotting_y=clip(plotting_y, .min=0, .max=50))

ggplot() + 
  geom_point(data=reds, aes(x=plotting_x, y=plotting_y, color=plotting_z), , size=3, shape=19, alpha=2/3, stroke=0) + 
  geom_point(data=blues, aes(plotting_x, plotting_y, color=plotting_z), size=3, shape=19, alpha=2/3, stroke=0) +
  geom_point(data=blacks, aes(plotting_x, plotting_y), color = "black", size=1, shape=19, alpha=0.2, stroke=0) +
  labs(x=expression('Expression difference [log'[2]~'FC]'), y=expression('Significance [log'[10]~'adjusted P]'), color=expression('Expression average [log'[10]~'TPM]')) + 
  theme_linedraw() +
  geom_segment(aes(x=-1, xend=-1, y=-log10(0.05), yend=50), linetype=2) +
  geom_segment(aes(x=1, xend=1, y=-log10(0.05), yend=50), linetype=2) +
  geom_segment(aes(x=-8, xend=-1, y=-log10(0.05), yend=-log10(0.05)), linetype=2) +
  geom_segment(aes(x=1, xend=8, y=-log10(0.05), yend=-log10(0.05)), linetype=2) +
  xlim(-8, 8) +
  scale_color_viridis_c(option = "cividis") 






