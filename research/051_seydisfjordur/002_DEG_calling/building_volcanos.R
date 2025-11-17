# Load library
library(ggplot2)


# Read data
file = '/Users/adrian/research/bmcbf/051_seydisfjordur/results/degs_deseq/effect_RES_vs_SEN.for.tsv'
response_df <- read.table(file, header = TRUE, sep = "\t")

# Basic volcano plot
ggplot(response_df, aes(x = log2FoldChange, y = clip(-log10(padj), .min=-6, .max=6))) +
  geom_point(alpha = 0.6) +
  theme_minimal() +
  labs(x = "log2(Fold Change)", y = "-log10(p-value)")

