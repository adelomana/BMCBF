#
# 1. installation
#
#install.packages("devtools")
#library(devtools)
#devtools::install_github("cit-bioinfo/mMCP-counter")

#
# 2. load libraries
#
library("mMCPcounter")
library(pheatmap)
library(readxl)

#
# 3. read data
#
expressionDataFile = '/Users/adrian/research/bmcbf/029_orsay/results/000.quantification/results/DESeq2_TPM_values.tsv'
expressionData = read.table(expressionDataFile, header = TRUE, sep = "\t", row.names = 1)
rownames(expressionData) <- sub("\\..*", "", rownames(expressionData))
expressionData <- as.matrix(expressionData)
View(expressionData)

# 
# 4. read metadata
#
filename = '/Users/adrian/research/bmcbf/029_orsay/metadata/SampleDescription-annotated.xlsx'
metadata <- read_excel(filename)
View(metadata)

# 5. change column names
name_map <- setNames(as.character(metadata$"mouse Nb"), as.character(metadata$"NAME"))
colnames(expressionData) <- name_map[colnames(expressionData)]
View(expressionData)
dim(expressionData)

#
# 6. remove samples as they are outliers or second extractions
#

# remove blue bc it looks away in TILs
dim(expressionData)
drop_these <- c("390B")
expressionData <- expressionData[, !(colnames(expressionData) %in% drop_these)]
dim(expressionData)

# remove two reds, 2105B looks very weird in TILs, 2239 was already flagged
dim(expressionData)
drop_these <- c("2239", '2105B')
expressionData <- expressionData[, !(colnames(expressionData) %in% drop_these)]
dim(expressionData)

# remove three yellows, 2409 seems close to red. Also 2657B without clear reasons
dim(expressionData)
drop_these <- c("2409", '2728', '2657B')
expressionData <- expressionData[, !(colnames(expressionData) %in% drop_these)]
dim(expressionData)

# remove five blacks. 376, 376B and 2727 are simply far way in PCA
# 2773B is close to red, as 2423
dim(expressionData)
drop_these <- c("376", '376B', '2727', '2773B', '2423')
expressionData <- expressionData[, !(colnames(expressionData) %in% drop_these)]
dim(expressionData)


#
# 6. estimate
#
# index in kallisto indexes github is 108 which is GRCm39. Seems that the flag name has a typo. gCr? It should be Genome Reference Consortium
immunoProfiles = mMCPcounter.estimate(expressionData, features = "ENSEMBL.ID", genomeVersion = "GCRm39")
View(immunoProfiles)

# 5. visualize
pheatmap(
  immunoProfiles,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "complete",
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize_row = 9,
  fontsize_col = 9,
  main = "mMCP-counter (row Z-scores)"
)

# lets drop some non immune cells to avoid clustering biases
drop_rows <- c("Vessels", "Endothelial cells", "Fibroblasts", 'Lymphatics')
immunoProfiles <- immunoProfiles[!(rownames(immunoProfiles) %in% drop_rows), ]

# create z score
mat_z <- t(scale(t(immunoProfiles), center = TRUE, scale = TRUE))

# define the colorbar
limit <- 4
breaks <- seq(-limit, limit, length.out = 100)
cols <- colorRampPalette(c("blue", "white", "red"))(length(breaks)-1)

# Define distance measures and clustering methods


# define group colors
metadata$"mouse Nb" <- as.character(metadata$"mouse Nb")
colnames(mat_z) <- as.character(colnames(mat_z))
idx <- match(colnames(mat_z), metadata$"mouse Nb")
annot <- data.frame(
  genotype = metadata$genotype[idx],
  row.names = colnames(mat_z)
)
ann_colors <- list(
  genotype = c(
    wt  = "grey70",
    hom = "red",
    het = "gold",     # if you have these
    "Vga/F" = "skyblue"   # if you have these
  )
)

# Loop through all combinations

distances <- c("euclidean", "maximum", "manhattan", "canberra", "minkowski",
               "correlation", "binary")
methods   <- c("single", "complete", "average", "mcquitty",
               "median", "centroid", "ward.D", "ward.D2")



d = 'euclidean'
m = 'ward.D2'
pheatmap(
  mat_z,
  color = cols,
  breaks = breaks,
  clustering_distance_rows = d,
  clustering_distance_cols = d,
  clustering_method = m,
  annotation_col = annot,
  annotation_colors = ann_colors,
  angle_col = 90,
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize_row = 12,
  fontsize_col = 12,
  main = paste("mMCP-counter z-scores [", d, " + ", m, "]")
)

#################### LOW GRANULAROTU

## 1. Define your groups of cell types (by row name in mat_z)
groups <- list(
  T_lineage = c("T cells", "CD8 T cells"),
  B_lineage = c("B derived", "Memory B cells"),
  NK        = c("NK cells"),
  Myeloid_Gran = c(
    "Monocytes / macrophages",
    "Monocytes",
    "Granulocytes",
    "Mast cells",
    "Eosinophils",
    "Neutrophils",
    "Basophils"
  )
)

## 2. For each group, compute the median across rows for each column
group_medians <- lapply(groups, function(rows) {
  # be robust in case some rows were dropped earlier
  keep_rows <- intersect(rows, rownames(mat_z))
  # if only one row, apply() still works
  apply(mat_z[keep_rows, , drop = FALSE], 2, mean, na.rm = TRUE)
})

## 3. Bind into a new matrix with 4 rows
mat_z_grouped <- do.call(rbind, group_medians)

## 4. Check result
View(mat_z_grouped)


pheatmap(
  mat_z_grouped,
  color = cols,
  breaks = breaks,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "ward.D2",
  annotation_col = annot,
  annotation_colors = ann_colors,
  angle_col = 90,
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize_row = 12,
  fontsize_col = 12,
  main = "Grouped mMCP-counter z-scores (median per lineage)"
)







