rm(list = ls())

#
# -1. install libraries
# 
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("DESeq2")

#
# 0. load libraries
#
library(DESeq2)
library(tximport)
library(crayon) 
library(ggplot2)
library(stringr)
library(ramify) # this is for clip, for the volcano
library(dplyr)

#
# 0. user-defined variables
#
setwd("/Users/adrian/scratch/")
kallisto_dir = "/Users/adrian/research/bmcbf/019_fossvogur/001_profiles"
results_dir = '/Users/adrian/research/bmcbf/019_fossvogur/03_degs'

#
# 1. generate gene to transcript mapping
#
t2g_file = '/Users/adrian/software/kallisto/human_index_standard/t2g.txt'
t2g = read.csv(t2g_file, sep='\t', header=FALSE)
dim(t2g)

annotation_file = '/Users/adrian/software/kallisto/human_index_standard/annotation.tsv'
full_annotation = read.csv(annotation_file, sep='\t')
annotation <- full_annotation %>% distinct(ensembl_gene_id, .keep_all = TRUE)

#
# 2. define full metadata
#
dirnames = list.dirs(kallisto_dir, full.names=TRUE, recursive=FALSE)
dirnames = dirnames[grep('_processed', dirnames)]
labels = sapply(strsplit(dirnames, split='/',fixed=TRUE), function(x) (x[8]))
labels = str_remove(labels, '_processed')
metadata = data.frame(labels)

paths = file.path(dirnames, paste0("kallisto_output_c"), 'abundance.h5')
metadata$path = paths

# drop the 501-background KO13 samples (duplicate KO13 label issue)
metadata <- metadata[!grepl("^501-", metadata$labels), ]
nrow(metadata)  # should be 28

# keep only siCTRL, KO24, KO13, C3_RD
keep <- c("siCTRL", "KO24", "KO13", "C3_RD")
pattern <- paste(keep, collapse = "|")
metadata <- metadata[grepl(pattern, metadata$labels), ]
nrow(metadata)  # should be 12

metadata$treatment = c(rep('C3', 3), rep('KO13', 3), rep('KO24', 3), rep('CTRL', 3))

View(metadata)

#
# 3. iterate contrasts
#
contrasts = list(
  c('KO13', 'C3'),
  c('KO24', 'CTRL')
)

for (contrast in contrasts) {
  sample_flag = contrast[1]
  control_flag = contrast[2]
  design_name = 'treatment'
  print(c('working with', sample_flag, control_flag))
  
  # 3.1. define working metadata
  rules = grepl(sample_flag, metadata$labels) | grepl(control_flag, metadata$labels)
  working_metadata <- metadata[rules, ]
  print('working metadata')
  print(working_metadata)
  
  # 3.2. read quantification files
  txi = tximport(working_metadata$path, type="kallisto", tx2gene=t2g)
  design_formula <- as.formula(paste("~", design_name))
  dds = DESeqDataSetFromTximport(txi, colData=working_metadata, design=design_formula) 
  
  print(c('and the reference is', control_flag))
  dds[[ design_name ]] <- relevel(dds[[ design_name ]], ref = control_flag)
  
  # 3.3. minimal filter on undetected genes: at least 10 reads in each biological group
  # additional filters should happen a posteriori, because removing a lot of genes would affect DESeq2 assumptions
  cat(blue(paste('size before counts filtering:', dim(dds)[1], sep=' ')), fill=TRUE)
  group <- droplevels(as.factor(colData(dds)$treatment))
  keep <- Reduce(`|`, lapply(levels(group), function(g)
    rowSums(counts(dds)[, group == g, drop = FALSE] >= 10) == sum(group == g)))
  dds <- dds[keep,]
  cat(blue(paste('size after counts filtering:', dim(dds)[1], sep=' ')), fill=TRUE)
  
  # 3.4. Fit model, either Wald or LRT test
  dds <- DESeq(dds, test = "LRT", reduced = ~ 1)
  #dds <- DESeq(dds, test = "Wald")
  resultsNames(dds)
  res = results(dds, alpha=0.05)
  
  # 3.5 Shrink LFC for stable effect-size reporting
  cat(blue('coefficient:', resultsNames(dds)[2]))
  res_shr <- lfcShrink(dds, coef = resultsNames(dds)[2], type = "apeglm")
  # visualize shrinkage
  plotMA(res, ylim = c(-5, 5), main=paste('no shrink', resultsNames(dds)[2]))
  plotMA(res_shr, ylim = c(-5, 5), main=paste('no shrink', resultsNames(dds)[2]))
  
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
  level_A <- unique(working_metadata[design_name])[[1]][1]   # typically sample         
  level_B <- unique(working_metadata[design_name])[[1]][2]   # typically control
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
  cat(blue(paste('full results dim:', dim(full_results)[1], sep=' ')), fill=TRUE)
  cat(blue(paste('no responders:', dim(no_responders)[1], sep=' ')), fill=TRUE)
  cat(blue(paste('DEGs:', sum(full_results$padj < 0.05, na.rm = TRUE),  sep=' ')), fill=TRUE)
  cat(blue(paste('responders:', dim(responders)[1],  sep=' ')), fill=TRUE)
  
  # store results
  responders_clean <- responders %>%
    select(
      `ENSEMBL ID`                    = gene_id,
      `Biotype`                       = biotype,
      `Gene symbol`                   = gene_name,
      `Description`                   = description2,
      `Log2 FC shr`                       = log2FC_shr,
      `Expression control [TPM]`      = median_TPM_A,
      `Expression treatment [TPM]`    = median_TPM_B,
      `Counts delta`                  = delta_counts,
      `Adjusted P`                    = padj
    )
  
  full_results_clean <- full_results %>%
    select(
      `ENSEMBL ID`                    = gene_id,
      `Biotype`                       = biotype,
      `Gene symbol`                   = gene_name,
      `Description`                   = description2,
      `Log2 FC`                       = log2FC_shr,
      `Expression control [TPM]`      = median_TPM_A,
      `Expression treatment [TPM]`    = median_TPM_B,
      `Counts delta`                  = delta_counts,
      `Adjusted P`                    = padj
    )
  
  dir.create(results_dir)
  filename = paste(results_dir, '/effect_', design_name, '_', sample_flag, '_vs_', control_flag, '.full.tsv', sep='')
  write.table(full_results_clean, file=filename, quote=FALSE, sep='\t', row.names = FALSE)
  
  filename = paste(results_dir, '/effect_', design_name, '_', sample_flag, '_vs_', control_flag, '.responders.tsv', sep='')
  write.table(responders_clean, file=filename, quote=FALSE, sep='\t', row.names = FALSE)
  
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
  
}









