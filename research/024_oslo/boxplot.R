# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("ggsignif")

library(ggplot2)
library(ggsignif)

control = c(5.02180304970648, 5.00880907480489, 5.01228320293699, 4.98936393050915, 4.97936008652491, 4.97121824693022,  4.97317314061993, 4.95493241732533,4.9609540735458,4.9600110087741,4.95072258764067)
treated = c(5.31457127822945, 5.30612736990571, 5.3041823790328,  5.29129560895088, 5.27806744847277, 5.28014710489925,  5.26175248601452, 5.24739641819883,5.22987147454147,5.25117588091659,5.23276471363254) 

wilcox.test(control, treated)
t.test(control, treated)

df = data.frame(group=rep(c("control", "treated"), each = length(control)), values = c(control, treated))
df

p1 = ggplot(df, aes(group, values, fill=group)) + 
  geom_boxplot() +
  geom_signif(comparisons=list(c('control', 'treated')), map_signif_level = TRUE) +
  theme_minimal() +
  scale_fill_manual(values = c("skyblue", "gold"))
p1
