# Differential H3K9me2 per gene, EpiLC vs ESC (2i/LIF), input-controlled.
#
# Model: ~ assay + cell + assay:cell
# The interaction coefficient (assayChIP.cellEpiLC) is the change in ChIP signal
# between EpiLC and ESC beyond any change in the input.
#
# Two normalizations, because H3K9me2 rises globally in EpiLC and there is no spike-in:
#   "median_ratio": DESeq2 default. Assumes most genes do not change, so it
#                   removes a global gain and can create artificial losses.
#   "library_size": size factors from total reads per sample. Keeps a global gain,
#                   assuming most reads are background.
# Genes significant under both are the robust set.
#
# Step 3.
# Input:  <workdir>/metadata/samples.tsv, <workdir>/results/counts_<regions>.tsv and counts_<regions>_totals.tsv
# Output: <workdir>/results/<regions>_H3K9me2_EpiLC_vs_ESC.tsv
#
# Usage: Rscript deseq2_h3k9me2.R gene_bodies
#        Rscript deseq2_h3k9me2.R promoters
#        Rscript deseq2_h3k9me2.R gene_bodies /other/path   (optional working directory)

suppressPackageStartupMessages({
  library(DESeq2)
  library(apeglm)
})

args <- commandArgs(trailingOnly = TRUE)
regions <- args[1]
workdir <- if (length(args) >= 2) args[2] else "/Users/adrian/research/bmcbf/102"
results_dir <- file.path(workdir, "results")
counts_file <- file.path(results_dir, paste0("counts_", regions, ".tsv"))
totals_file <- file.path(results_dir, paste0("counts_", regions, "_totals.tsv"))
samples_file <- file.path(workdir, "metadata", "samples.tsv")
prefix <- file.path(results_dir, regions)

counts <- as.matrix(read.delim(counts_file, row.names = 1, check.names = FALSE))
coldata <- read.delim(samples_file, row.names = 1)
coldata <- coldata[colnames(counts), ]
coldata$cell <- factor(coldata$cell, levels = c("ESC", "EpiLC"))
coldata$assay <- factor(coldata$assay, levels = c("Input", "ChIP"))
totals <- read.delim(totals_file, row.names = 1)[colnames(counts), "total_reads"]

# Keep genes with at least 10 reads in at least 2 ChIP samples
chip <- coldata$assay == "ChIP"
keep <- rowSums(counts[, chip] >= 10) >= 2
counts <- counts[keep, ]
message(sum(keep), " genes kept after filtering")

coef_name <- "assayChIP.cellEpiLC"

run <- function(norm) {
  dds <- DESeqDataSetFromMatrix(counts, coldata, design = ~ assay + cell + assay:cell)
  if (norm == "library_size") {
    sizeFactors(dds) <- totals / exp(mean(log(totals)))
  } else {
    dds <- estimateSizeFactors(dds)
  }
  dds <- DESeq(dds, quiet = TRUE)
  res <- results(dds, name = coef_name, alpha = 0.05)
  shr <- lfcShrink(dds, coef = coef_name, type = "apeglm", quiet = TRUE)
  out <- data.frame(
    gene = rownames(res),
    baseMean = res$baseMean,
    log2FC = res$log2FoldChange,
    log2FC_shrunk = shr$log2FoldChange,
    pvalue = res$pvalue,
    padj = res$padj
  )
  colnames(out)[-1] <- paste0(colnames(out)[-1], "_", norm)
  list(table = out, dds = dds)
}

a <- run("median_ratio")
b <- run("library_size")
tab <- merge(a$table, b$table, by = "gene")

# Normalized ChIP / input per cell type (library-size normalization) for reading the table
nc <- counts(b$dds, normalized = TRUE)
for (cl in c("ESC", "EpiLC")) {
  ch <- rowMeans(nc[, coldata$cell == cl & coldata$assay == "ChIP", drop = FALSE])
  inp <- rowMeans(nc[, coldata$cell == cl & coldata$assay == "Input", drop = FALSE])
  tab[[paste0("log2_ChIP_over_input_", cl)]] <- log2((ch + 1) / (inp + 1))[tab$gene]
}

sig <- function(x, n) !is.na(x[[paste0("padj_", n)]]) & x[[paste0("padj_", n)]] < 0.05 &
  abs(x[[paste0("log2FC_shrunk_", n)]]) > 1
tab$robust <- sig(tab, "median_ratio") & sig(tab, "library_size") &
  sign(tab$log2FC_shrunk_median_ratio) == sign(tab$log2FC_shrunk_library_size)
tab$direction <- ifelse(tab$log2FC_shrunk_library_size > 0, "gain_in_EpiLC", "loss_in_EpiLC")
tab <- tab[order(tab$padj_library_size), ]

write.table(tab, paste0(prefix, "_H3K9me2_EpiLC_vs_ESC.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

cat("\nSignificant genes (padj < 0.05, |shrunk log2FC| > 1):\n")
for (n in c("median_ratio", "library_size")) {
  s <- sig(tab, n)
  up <- sum(s & tab[[paste0("log2FC_shrunk_", n)]] > 0)
  cat(sprintf("  %-13s gain: %5d  loss: %5d\n", n, up, sum(s) - up))
}
cat(sprintf("  robust (both, same direction): %d\n", sum(tab$robust)))

# Sanity checks from the paper (Fig 5C): expect gains in EpiLC
cat("\nPositive controls:\n")
print(tab[grepl("\\|(Asz1|Cdkn1a)$", tab$gene),
          c("gene", "log2FC_shrunk_median_ratio", "padj_median_ratio",
            "log2FC_shrunk_library_size", "padj_library_size")])
