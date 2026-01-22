#!/usr/bin/env Rscript

rm(list = ls())

suppressPackageStartupMessages({
  library(ChIPseeker)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  library(org.Hs.eg.db)
  library(GenomicRanges)
  library(AnnotationDbi)
})

txdb  <- TxDb.Hsapiens.UCSC.hg38.knownGene
orgdb <- org.Hs.eg.db

input_dir  <- "/Users/adrian/research/bmcbf/004_keilir/results/03_necio_viz"
output_dir <- "/Users/adrian/research/bmcbf/004_keilir/results/04_necio_f2_score"

tss_windows_to_try <- list(
  "c150b" = c(-150, 25),
  "c300b" = c(-300, 50),
  "c500b" = c(-500, 0),
  "s500b" = c(-500, 100),
  "c1kb" = c(-1000, 100),
  "s1kb" = c(-1000, 500),
  "2kb"  = c(-2000, 2000),
  "5kb" = c(-5000, 5000),
  "10kb" = c(-10000, 10000),
  "25kb" = c(-25000, 25000),
  "50kb" = c(-50000, 50000)
)

bed_files <- list.files(path = input_dir, pattern = "\\.bed$", full.names = TRUE)
if (length(bed_files) == 0) stop("No .bed files found in input_dir: ", input_dir)

## Precompute gene TSS windows per window size (faster, and easier to debug)
g <- genes(txdb)  # names(g) are Entrez IDs in knownGene TxDb

for (bed in bed_files) {
  base <- sub("\\.bed$", "", basename(bed))
  
  peaks <- readPeakFile(bed)
  cat("\n=== File:", basename(bed), "peaks:", length(peaks), "===\n")
  cat("Peak seqlevels (first 20):", paste(head(seqlevels(peaks), 20), collapse = ","), "\n")
  
  for (wname in names(tss_windows_to_try)) {
    tss_window <- tss_windows_to_try[[wname]]
    up <- abs(tss_window[1])
    dn <- abs(tss_window[2])
    
    ## Gene-centered TSS windows (gene IDs are names(tss))
    tss <- promoters(g, upstream = up, downstream = dn)
    
    ## Ensure common seqlevels (avoids silent 0 overlaps if seqlevels differ)
    common <- intersect(seqlevels(peaks), seqlevels(tss))
    if (length(common) == 0) {
      cat("Window", wname, "common seqlevels: 0 -> genes: 0\n")
      out <- file.path(output_dir, paste0(base, ".genes.multi.", wname, ".txt"))
      writeLines(character(0), con = out)
      next
    }
    peaks2 <- keepSeqlevels(peaks, common, pruning.mode = "coarse")
    tss2   <- keepSeqlevels(tss,   common, pruning.mode = "coarse")
    
    hits <- findOverlaps(peaks2, tss2, ignore.strand = TRUE)
    cat("Window", wname, "hits:", length(hits), "\n")
    
    if (length(hits) == 0) {
      genes <- character(0)
    } else {
      entrez <- names(tss2)[subjectHits(hits)]
      entrez <- unique(entrez)
      entrez <- entrez[!is.na(entrez) & entrez != ""]
      cat("Window", wname, "unique Entrez:", length(entrez), "\n")
      
      if (length(entrez) == 0) {
        genes <- character(0)
      } else {
        symbol <- AnnotationDbi::mapIds(
          orgdb,
          keys      = entrez,
          column    = "SYMBOL",
          keytype   = "ENTREZID",
          multiVals = "first"
        )
        genes <- unique(unname(symbol[!is.na(symbol) & symbol != ""]))
      }
    }
    
    genes <- sort(unique(genes))
    out <- file.path(output_dir, paste0(base, ".genes.multi.", wname, ".txt"))
    writeLines(genes, con = out)
    
    cat("Wrote", out, "genes:", length(genes), "\n")
  }
}
