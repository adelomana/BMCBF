# Load library
library(ggplot2)

# Read data
# getwd() and setwd() are good functions to locate the path to the data in your computer
response_df <- read.table("effect_RES_vs_SEN.raw.tsv", header = TRUE, sep = "\t")
View(response_df)

# Basic volcano plot
ggplot(response_df, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(alpha = 0.6) +
  theme_minimal() +
  labs(x = "log2(Fold Change)", y = "-log10(p-value)")

# think. How would you improve this visualization to convey a clear scientific message?

# this is a deeper modification approach to represent the same effect
library(ramify) # this is for clip, for the volcano
library(matrixStats) # this is required to call rowMedians()

noresponse_df <- read.table("effect_RES_vs_SEN.anti.tsv", header = TRUE, sep = "\t")

plotting_x = response_df$log2FoldChange
y = response_df$padj
epsilon = min(y[y !=0])
plotting_y = -log10(y + epsilon) 
df = data.frame(plotting_x=clip(plotting_x, .min=-6, .max=6), plotting_y=clip(plotting_y, .min=0, .max=20))
reds = df[df$plotting_x > 0, ]
blues = df[df$plotting_x < 0, ]

plotting_x = noresponse_df$log2FoldChange
plotting_y = -log10(noresponse_df$padj)
blacks = data.frame(plotting_x=clip(plotting_x, .min=-6, .max=6), plotting_y=clip(plotting_y, .min=0, .max=20))

ggplot() + 
  geom_point(data=reds, aes(x=plotting_x, y=plotting_y), , size=3, shape=19, alpha=2/3, stroke=0, color='red') + 
  geom_point(data=blues, aes(plotting_x, plotting_y), size=3, shape=19, alpha=2/3, stroke=0, color='blue') +
  geom_point(data=blacks, aes(plotting_x, plotting_y), color = "black", size=1, shape=19, alpha=0.2, stroke=0) +
  labs(x=expression('Expression [log'[2]~'FC]'), y=expression('Significance [log'[10]~'adjusted P]')) + 
  theme_linedraw() +
  geom_segment(aes(x=-1, xend=-1, y=-log10(0.05), yend=20), linetype=2) +
  geom_segment(aes(x=1, xend=1, y=-log10(0.05), yend=20), linetype=2) +
  geom_segment(aes(x=-6, xend=-1, y=-log10(0.05), yend=-log10(0.05)), linetype=2) +
  geom_segment(aes(x=1, xend=6, y=-log10(0.05), yend=-log10(0.05)), linetype=2) +
  xlim(-6, 6)
  
