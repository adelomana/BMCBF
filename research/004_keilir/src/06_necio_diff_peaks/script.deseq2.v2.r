## ------------------------------------------------------------------
## Session hygiene: start from a clean environment
## ------------------------------------------------------------------

rm(list = ls(all.names = TRUE))  # remove all objects, including hidden ones
gc()                             # trigger garbage collection

## Optional but recommended for reproducibility
options(stringsAsFactors = FALSE)


## Differentially bound peaks (DBPs) with DESeq2
## Input: featureCounts output: counts_featureCounts.black.txt
## Peaks: black.AM.union.saf / <peak_label>.AM.union.bed

library(DESeq2)
library(ggplot2)
library(crayon) 
library(ramify) # this is for clip, for the volcano

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(GenomicFeatures)
  library(GenomeInfoDb)
  library(AnnotationDbi)
  library(ChIPseeker)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  library(org.Hs.eg.db)
})


annotate_peaks_with_promoter_genes <- function(res_annot_df,
                                               txdb,
                                               orgdb,
                                               upstream = 1000,
                                               downstream = 500,
                                               genes_colname = "promoter_genes_s1kb") {
  
  ## Full peaks GRanges (do NOT drop peaks here)
  peaks_gr_all <- GenomicRanges::GRanges(
    seqnames = res_annot_df$Chr,
    ranges   = IRanges::IRanges(start = res_annot_df$Start,
                                end   = res_annot_df$End),
    strand   = "*"
  )
  names(peaks_gr_all) <- res_annot_df$Geneid
  
  ## Output vector aligned to res_annot_df (default NA)
  out <- rep(NA_character_, nrow(res_annot_df))
  names(out) <- res_annot_df$Geneid
  
  ## Transcript promoters
  tx <- GenomicFeatures::transcripts(txdb)
  tx_names <- as.character(mcols(tx)$tx_name)
  
  tss <- GenomicFeatures::promoters(tx, upstream = upstream, downstream = downstream)
  names(tss) <- tx_names
  tss <- GenomicRanges::trim(tss)
  
  ## Transcript -> gene (Entrez)
  tx2gene <- AnnotationDbi::select(
    txdb,
    keys    = tx_names,
    keytype = "TXNAME",
    columns = "GENEID"
  )
  
  ## Harmonize seqlevels style (peaks -> tss style)
  if (length(GenomeInfoDb::seqlevelsStyle(tss)) > 0) {
    suppressWarnings(
      GenomeInfoDb::seqlevelsStyle(peaks_gr_all) <- GenomeInfoDb::seqlevelsStyle(tss)[1]
    )
  }
  
  ## Filter to standard chromosomes for overlap computation ONLY
  standard <- GenomeInfoDb::standardChromosomes(txdb)
  
  standard_peaks <- intersect(GenomeInfoDb::seqlevels(peaks_gr_all), standard)
  standard_tss   <- intersect(GenomeInfoDb::seqlevels(tss),         standard)
  
  peaks_gr <- GenomeInfoDb::keepSeqlevels(peaks_gr_all, standard_peaks, pruning.mode = "coarse")
  tss_use  <- GenomeInfoDb::keepSeqlevels(tss,          standard_tss,   pruning.mode = "coarse")
  
  ## Keep only common seqlevels
  common <- intersect(GenomeInfoDb::seqlevels(peaks_gr), GenomeInfoDb::seqlevels(tss_use))
  if (length(common) == 0) {
    res_annot_df[[genes_colname]] <- out
    return(res_annot_df)
  }
  
  peaks2 <- GenomeInfoDb::keepSeqlevels(peaks_gr, common, pruning.mode = "coarse")
  tss2   <- GenomeInfoDb::keepSeqlevels(tss_use,  common, pruning.mode = "coarse")
  
  ## Overlaps
  hits <- GenomicRanges::findOverlaps(peaks2, tss2, ignore.strand = TRUE)
  if (length(hits) == 0) {
    res_annot_df[[genes_colname]] <- out
    return(res_annot_df)
  }
  
  tx_hits <- names(tss2)[subjectHits(hits)]
  pk_hits <- names(peaks2)[queryHits(hits)]  # these are Geneid strings
  
  ## Transcript -> gene
  m <- match(tx_hits, tx2gene$TXNAME)
  gene_hits <- tx2gene$GENEID[m]
  
  ## Filter gene_hits and pk_hits together (critical)
  keep <- !is.na(gene_hits) & gene_hits != "" & !is.na(pk_hits) & pk_hits != ""
  gene_hits <- gene_hits[keep]
  pk_hits   <- pk_hits[keep]
  
  if (length(gene_hits) == 0) {
    res_annot_df[[genes_colname]] <- out
    return(res_annot_df)
  }
  
  ## Group Entrez gene IDs by peak (multi-hit)
  entrez_by_peak <- split(gene_hits, pk_hits)
  
  ## Entrez -> SYMBOL
  all_entrez <- unique(unname(unlist(entrez_by_peak)))
  
  sym_map <- AnnotationDbi::mapIds(
    orgdb,
    keys      = all_entrez,
    column    = "SYMBOL",
    keytype   = "ENTREZID",
    multiVals = "first"
  )
  
  ## Collapse per peak and write into 'out' by peak Geneid
  for (pk in names(entrez_by_peak)) {
    ent  <- unique(entrez_by_peak[[pk]])
    syms <- unique(unname(sym_map[ent]))
    syms <- syms[!is.na(syms) & syms != ""]
    if (length(syms) > 0) {
      out[pk] <- paste(sort(unique(syms)), collapse = ";")
    }
  }
  
  ## Assign aligned vector back to data frame
  res_annot_df[[genes_colname]] <- out
  return(res_annot_df)
}

### end annotate peaks function

## ------------------------------------------------------------------
## CONFIGURATION: change only here for purple / orange and directories
## ------------------------------------------------------------------

peak_label <- "black" 

base_dir   <- "/Users/adrian/research/bmcbf/004_keilir/results/06_necio_diff_peaks"
setwd(base_dir)

## all outputs will go here, e.g. results_purple, results_orange
output_dir <- file.path(base_dir, paste0("results_", peak_label))
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

cat(blue(paste0("Running analysis for label: ", peak_label)), fill = TRUE)
cat(blue(paste0("Base directory: ", base_dir)), fill = TRUE)
cat(blue(paste0("Output directory: ", output_dir)), fill = TRUE)

seta_indexes    = 1:3
setb_indexes    = 4:6

## 1. Read featureCounts table
fc_file <- "/Users/adrian/research/bmcbf/004_keilir/results/05_elja_count_retriever/counts_featureCounts.txt"
cat(blue(paste0("Reading featureCounts file: ", fc_file)), fill = TRUE)

fc <- read.table(
  fc_file,
  header           = TRUE,
  sep              = "\t",
  comment.char     = "#",
  stringsAsFactors = FALSE,
  check.names      = FALSE  # keep full paths as column names initially
)

cat("Columns in raw featureCounts table:\n")
print(colnames(fc))

## 2. Clean sample names: keep the directory name just before the BAM filename
## Example:
## /hpcdata/.../MITF_A_Untreated_FLAG_1/human.30_120.bam
## becomes: MITF_A_Untreated_FLAG_1

count_col_idx    <- 7:ncol(fc)
raw_sample_cols  <- colnames(fc)[count_col_idx]

extract_sample_name <- function(path) {
  parts <- strsplit(path, "/", fixed = TRUE)[[1]]
  if (length(parts) >= 2) {
    return(parts[length(parts) - 1])  # directory just before BAM file
  } else {
    return(path)
  }
}

sample_names <- vapply(raw_sample_cols, extract_sample_name, character(1))

cat("\nRenaming count columns to:\n")
print(sample_names)

colnames(fc)[count_col_idx] <- sample_names

## 3. Separate annotation and count matrix
peak_annot <- fc[, c("Geneid","Chr","Start","End","Strand")]

count_mat <- as.matrix(fc[, count_col_idx])
storage.mode(count_mat) <- "integer"
rownames(count_mat) <- fc$Geneid

cat("\nDimension of count matrix:\n")
print(dim(count_mat))
#View(count_mat)

## 4. Build sample metadata (A vs M)
## Here: MITF_A_Untreated_FLAG_* = condition A
##       MITF_M_Untreated_FLAG_* = condition M

conditions <- ifelse(grepl("MITF_A", sample_names), "A", "M")

## Make M the reference level
condition_factor <- factor(conditions, levels = c("M","A"))

coldata <- data.frame(
  row.names = sample_names,
  condition = condition_factor
)

cat("\nSample metadata:\n")
print(coldata)
#View(coldata)

## 5. Construct DESeqDataSet
dds <- DESeqDataSetFromMatrix(
  countData = count_mat,
  colData   = coldata,
  design    = ~ condition
)

## 7. Run DESeq2
dds <- DESeq(dds)

## 8. Extract results: condition A vs M (M is reference)
res <- results(dds, contrast = c("condition","A","M"))

# manipulate results for later
length(rowMedians(counts(dds)[ , 1:3]))

res_df           <- as.data.frame(res)
res_df$Geneid    <- rownames(res_df)
res_df$countsa   <- rowMedians(counts(dds)[ , 1:3])
res_df$countsb   <- rowMedians(counts(dds)[ , 4:6])
res_df$counts_diff <- res_df$countsa - res_df$countsb

## 9. Attach genomic coordinates
res_annot <- merge(
  res_df,
  peak_annot,
  by   = "Geneid",
  sort = FALSE
)

## Reorder columns for readability and sort on adj P
res_annot <- res_annot[, c(
  "Geneid","Chr","Start","End","Strand","log2FoldChange","padj", 'countsa', 'countsb', 'counts_diff'
)]

res_annot <- res_annot[order(res_annot$padj), ]

cat("\nSummary of DESeq2 results:\n")
print(summary(res))


# --- Add promoter gene annotation (s1kb: -1000/+500 around TSS) ---
txdb  <- TxDb.Hsapiens.UCSC.hg38.knownGene
orgdb <- org.Hs.eg.db

res_annot <- annotate_peaks_with_promoter_genes(
  res_annot_df = res_annot,
  txdb         = txdb,
  orgdb        = orgdb,
  upstream     = 1000,
  downstream   = 500,
  genes_colname = "promoter_genes_s1kb"
)


## 10. Define strict DBPs
## Here: padj < 0.01 and |log2FC| > 1
dbp_strict <- subset(
  res_annot,
  !is.na(padj) &
    padj < 0.01 &
    abs(log2FoldChange) > 1 & counts_diff > 20
)

cat("\nNumber of strict DBPs (padj < 0.01 & |log2FC| > 1 & counts_diff > 20 :\n")
print(nrow(dbp_strict))

## Split by direction
dbp_A_up <- subset(dbp_strict, log2FoldChange > 0)  # gained in A
dbp_M_up <- subset(dbp_strict, log2FoldChange < 0)  # gained in M

cat("\nGained in A (log2FC > 0):\n")
print(nrow(dbp_A_up))

cat("Gained in M (log2FC < 0):\n")
print(nrow(dbp_M_up))

## 11. Write outputs (all labeled and inside output_dir)

## Full annotated table
write.table(
  res_annot,
  file      = file.path(output_dir, paste0("DESeq2_results_all_peaks.", peak_label, ".tsv")),
  sep       = "\t",
  quote     = FALSE,
  row.names = FALSE
)


## Strict DBPs: gained in A (log2FC > 0)
write.table(
  dbp_A_up[, c("Chr","Start","End","Geneid","log2FoldChange","padj","promoter_genes_s1kb")],
  file      = file.path(output_dir, paste0("DBPs_strict_A_gained.", peak_label, ".bed")),
  sep       = "\t",
  quote     = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

## Strict DBPs: gained in M (log2FC < 0)
write.table(
  dbp_M_up[, c("Chr","Start","End","Geneid","log2FoldChange","padj","promoter_genes_s1kb")],
  file      = file.path(output_dir, paste0("DBPs_strict_M_gained.", peak_label, ".bed")),
  sep       = "\t",
  quote     = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

##
## Volcano plot (A vs M, M reference)
##
plotting_x <- res_annot$log2FoldChange
y          <- res_annot$padj
plotting_y <- -log10(y) 
z          <- log10(rowMedians(as.matrix(res_annot[, c("countsa", "countsb")])) + 1)
print(c(min(z), max(z)))

df    <- data.frame(
  plotting_x = clip(plotting_x, .min=-2.5, .max=2.5),
  plotting_y = clip(plotting_y, .min=0, .max=60),
  plotting_z = clip(z, .min=0, .max=4.5)
)
reds  <- df[(df$plotting_x > 1)  & (df$plotting_y > -log10(0.01)), ]
blues <- df[(df$plotting_x < -1) & (df$plotting_y > -log10(0.01)), ]
blacks <- df[((df$plotting_x > -1) & (df$plotting_x < 1)) | (df$plotting_y < -log10(0.01)), ]

print(c(dim(reds)[1], dim(blues)[1], dim(dbp_strict)[1]))



ggplot() + 
  geom_point(data=reds,  aes(x=plotting_x, y=plotting_y, color=plotting_z),
             size=3.5, shape=19, alpha=1/3, stroke=0) + 
  geom_point(data=blues, aes(plotting_x, plotting_y, color=plotting_z),
             size=3.5, shape=19, alpha=1/3, stroke=0) +
  geom_point(data=blacks, aes(plotting_x, plotting_y),
             size=1.2, shape=19, alpha=1/6, stroke=0, color="black") + 
  labs(
    x = expression('Peak difference [log'[2]~'FC]'),
    y = expression('Significance [log'[10]~'adjusted P]')
  ) +
  theme_linedraw(base_size = 22) +
  theme(
    axis.title.x = element_text(size = 34, face = "bold"),
    axis.title.y = element_text(size = 34, face = "bold"),
    axis.text.x  = element_text(size = 28),
    axis.text.y  = element_text(size = 28),
    axis.ticks   = element_line(size = 1.2),
    axis.ticks.length = unit(0.35, "cm"),
    legend.title = element_text(size = 28),
    legend.text  = element_text(size = 24)
  ) +
  geom_segment(aes(x=-1,  xend=-1,  y=-log10(0.01), yend=60), linetype=2) +
  geom_segment(aes(x=1,   xend=1,   y=-log10(0.01), yend=60), linetype=2) +
  geom_segment(aes(x=-2.5,xend=-1,  y=-log10(0.01), yend=-log10(0.01)), linetype=2) +
  geom_segment(aes(x=1,   xend=2.5, y=-log10(0.01), yend=-log10(0.01)), linetype=2) +
  scale_x_continuous(breaks = c(-2, -1, 0, 1, 2), limits = c(-2.5, 2.5)) +
  scale_color_viridis_c(option = "viridis", name = "log10 Counts", limits = c(1.8, 4.1))

