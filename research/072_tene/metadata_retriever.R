if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("TCGAbiolinks")

library(TCGAbiolinks)
query <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)
GDCdownload(query)
data <- GDCprepare(query)

library(SummarizedExperiment)
counts <- assay(data, "unstranded")

# protein coding
BiocManager::install("AnnotationHub")
BiocManager::install("ensembldb")
library(AnnotationHub)
ah <- AnnotationHub()
query(ah, c("EnsDb", "Homo sapiens", "102"))

edb <- ah[["AH89180"]]  # use whatever ID is returned by query above

rownames(counts) <- gsub("\\..*", "", rownames(counts))
gene_biotype <- mapIds(
  edb,
  keys = rownames(counts),
  keytype = "GENEID",
  column = "GENEBIOTYPE",
  multiVals = "first"
)
counts <- counts[which(gene_biotype == "protein_coding"), ]
dim(counts)
View(counts)

write.table(counts, file = "counts.tsv", sep = "\t", row.names = TRUE, quote = FALSE)
