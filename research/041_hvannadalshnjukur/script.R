#
# installation
#
# if (!requireNamespace("BiocManager", quietly = TRUE)){
#   install.packages("BiocManager")
#   BiocManager::install()
# }
# BiocManager::install("IsoformSwitchAnalyzeR")

#
# loading
#
library(IsoformSwitchAnalyzeR)
packageVersion('IsoformSwitchAnalyzeR')

#
# loading expression quantification
#
quantification <- importIsoformExpression(parentDir = '/Users/adrian/research/bmcbf/041_clem/kallisto_output')
head(quantification$abundance, 2)
head(quantification$counts, 2)

#
# define experimental design
#
myDesign <- data.frame(
  sampleID = colnames(quantification$abundance)[-1],
  condition = rep(c('wt', 'ko'), 3)
)
myDesign

#
# load annotation
#
# wget https://ftp.ensembl.org/pub/release-108/gtf/homo_sapiens/Homo_sapiens.GRCh38.108.gtf.gz
# wget https://ftp.ensembl.org/pub/release-108/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz
#
aSwitchList <- importRdata(
  isoformCountMatrix   = quantification$counts,
  isoformRepExpression = quantification$abundance,
  designMatrix         = myDesign,
  isoformExonAnnoation = '/Users/adrian/research/bmcbf/041_clem/Homo_sapiens.GRCh38.108.gtf.gz',
  isoformNtFasta       = '/Users/adrian/research/bmcbf/041_clem/Homo_sapiens.GRCh38.cdna.all.fa.gz',
  fixStringTieAnnotationProblem = TRUE,
  showProgress = TRUE,
  ignoreAfterPeriod = TRUE
)
summary(aSwitchList)
head(aSwitchList$isoformFeatures,2)

#
# filtering
#
filteredlist <- preFilter(
  switchAnalyzeRlist = aSwitchList,
  geneExpressionCutoff = 10,
  isoformExpressionCutoff = 3,
  removeSingleIsoformGenes = TRUE
)

#
# identify isoform changes
#
switchListAnalyzed <- isoformSwitchTestDEXSeq(
  switchAnalyzeRlist = filteredlist,
  reduceToSwitchingGenes=TRUE
)
extractSwitchSummary(switchListAnalyzed)

