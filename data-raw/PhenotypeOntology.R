setwd('data-raw/')
#04032019
getwd()
pheno2gene<- read.delim("ALL_SOURCES_ALL_FREQUENCIES_phenotype_to_genes_28022019_12022019.txt", header = FALSE, sep = "\t")
#HPO-ID to HPO_ID, HPO-Name to HPO_Name, Gene-ID to GeneID - 13022019
colnames(pheno2gene) <- c("HPO_ID","HPO_Name","GeneID","GeneSymbol")
#remove first row of #Format: HPO-ID<tab>HPO-Name<tab>Gene-ID<tab>Gene-Name
head(pheno2gene)
pheno2gene <- pheno2gene[-1,]
head(pheno2gene)
#removal of HP terms Clinical modifier, Clinical course, Mode of inheritance, Frequency,Blood group
HPO_removal_list <- read.csv("HPO_removal_list_05042019.csv")
head(HPO_removal_list,5)
summary(HPO_removal_list) #Summary check is okay - 170
library(tibble)
library(dplyr)
pheno2gene_tbl <- as_tibble(pheno2gene)
pheno2gene_tbl
#483,449 rows - only phenotypic abnormality terms kept
pheno2gene_tbl <- pheno2gene_tbl %>%
  filter(!HPO_ID %in% HPO_removal_list$value)
#472,515 rows
#return to dataframe
pheno2gene <- as.data.frame(pheno2gene_tbl)
summary(pheno2gene) #Summary - no NA's
nrow(pheno2gene) #472515
head(pheno2gene,5)

usethis::use_data(pheno2gene, overwrite = TRUE)

summary(pheno2gene)

#Session_Info
sessionInfo()
