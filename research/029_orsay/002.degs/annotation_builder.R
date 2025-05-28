#
# 1. generate gene to transcript mapping
#
mart = biomaRt::useMart(biomart="ENSEMBL_MART_ENSEMBL", 
                        dataset="mmusculus_gene_ensembl",
                        #host = 'https://oct2022.archive.ensembl.org', 
                        # Ensembl 108 Oct 2022 https://oct2022.archive.ensembl.org     108
                        # last time I ran this I had to call fist host = 'https://www.ensembl.org', then it worked
                        verbose = TRUE)
# attributes = listAttributes(mart)
# hgnc_symbol gives less than external_gene_name
working_attributes = c('ensembl_gene_id', 
                       'external_gene_name',
                       'gene_biotype',
                       'description')
annotation = biomaRt::getBM(attributes=working_attributes, 
                     mart=mart,
                     verbose=TRUE)
dim(annotation) # Version 108 gives 274081 entries. Version 113 gives 412034
View(annotation)

store = paste('/Users/adrian/software/kallisto/mouse_index_standard', '/annotation.tsv', sep='')
write.table(annotation, file=store, quote=FALSE, sep='\t', col.names=NA)
