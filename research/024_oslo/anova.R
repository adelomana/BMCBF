v1 = c(4.8, 4.9, 5, 5.1, 5.2)
v2 = c(9.8, 9.9, 10, 10.1, 10.2)
v3 = c(3.8, 3.9, 4, 4.1, 4.2)
v4 = c(14.8, 14.9, 15, 15.1, 15.2)

df = data.frame(ocr=c(v1, v2, v3, v4), 
                 genotype=c(rep('wt', 10), rep('mut', 10)), 
                 treatment=c(rep('no', 5), rep('yes', 5), rep('no', 5), rep('yes', 5)))
df

result = aov(ocr ~ genotype + treatment + genotype:treatment, data=df)
summary(result)