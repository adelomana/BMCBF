# Or more modernly, use the cBioPortalData package
BiocManager::install("cBioPortalData")
library(cBioPortalData)

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
