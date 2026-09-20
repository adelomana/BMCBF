if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("TCGAbiolinks")

library(TCGAbiolinks)
library(SummarizedExperiment)
library(AnnotationHub)
library(dplyr)
library(cBioPortalData)


query <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)
#GDCdownload(query)
data <- GDCprepare(query)
counts <- assay(data, "unstranded")

# protein coding
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


#### metadata
clin <- as.data.frame(colData(data))
df_clinical <- data.frame(
  sample      = rownames(clin),
  patient     = clin$patient,
  age         = clin$age_at_diagnosis,
  subtype     = clin$paper_BRCA_Subtype_PAM50,
  sample_type = clin$shortLetterCode,
  row.names   = rownames(clin)
)
df_clinical <- df_clinical[!is.na(df_clinical$age) & df_clinical$sample_type == "TP", ]

cbio <- cBioPortal()
clinical_cbio <- clinicalData(cbio, studyId = "brca_tcga")

df_receptor <- clinical_cbio[, c("patientId", "ER_STATUS_BY_IHC")]
df_receptor <- df_receptor[df_receptor$ER_STATUS_BY_IHC %in% c("Positive", "Negative"), ]

df_receptor %>%
  group_by(patientId) %>%
  summarise(n_unique_ER = n_distinct(ER_STATUS_BY_IHC)) %>%
  filter(n_unique_ER > 1)

df_receptor <- df_receptor %>%
  distinct(patientId, .keep_all = TRUE)

table(df_receptor$ER_STATUS_BY_IHC)

df_clinical <- left_join(df_clinical, df_receptor, by = c("patient" = "patientId"))
table(df_clinical$subtype, df_clinical$ER_STATUS_BY_IHC, useNA = "always")

# Check discordance
df_clinical$discordant <- ifelse(
  (df_clinical$ER_STATUS_BY_IHC == "Positive" & df_clinical$subtype %in% c("Basal")) |
    (df_clinical$ER_STATUS_BY_IHC == "Negative" & df_clinical$subtype %in% c("LumA", "LumB")),
  TRUE, FALSE
)
df_clinical_clean <- df_clinical[
  df_clinical$discordant == FALSE &
    !is.na(df_clinical$discordant) &
    df_clinical$subtype %in% c("LumA", "LumB", "Basal"), ]

table(df_clinical_clean$subtype, df_clinical_clean$ER_STATUS_BY_IHC)
