# get the list of mouse innate immune genes

library(readxl)

mouse_innate_genes <- read_excel('innatedb_curated_genes.xls')

#species 10090 is mouse #species 9606 is human
mouse_innate_genes <- mouse_innate_genes %>%
    filter(Species == '10090')

write.csv(mouse_innate_genes, 'mouse_innate_genes.csv') 
