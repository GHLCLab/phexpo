library(ontologyIndex)
library(tibble)
library(magrittr)
library(dplyr)
setwd('data-raw/')
#data(hpo)
#####################################################################################
#HPO - OBO 28/02/2019
ontology <- get_ontology("hp28022019_12022019.obo")
ontology
check(ontology, stop_if_invalid = FALSE) #consistency checks - no error

#####################################################################################
#Clinical modifier
#Should be 82
Clinical_modifier <- get_descendants(ontology, roots="HP:0012823")
length(Clinical_modifier) #83 including HP:0012823
#Clinical course
Clinical_course <- get_descendants(ontology, roots="HP:0031797")
length(Clinical_course) #49 including HP:0031797
#Mode of inheritance
M_o_inheritance <- get_descendants(ontology, roots="HP:0000005")
length(M_o_inheritance) #29 including HP:0000005
#Frequency
Freq <- get_descendants(ontology, roots="HP:0040279")
length(Freq) #7 including HP:0040279
#Blood group
Bld_grp <- get_descendants(ontology, roots="HP:0032223")
length(Bld_grp) #2 including HP:0032223

#####Total

clinical_modifier_tbl <- enframe(Clinical_modifier)
clinical_modifier_tbl <- clinical_modifier_tbl %>%
  mutate(parent = "Clinical modifier")

Clinical_course_tbl <- enframe(Clinical_course)
Clinical_course_tbl <- Clinical_course_tbl %>%
  mutate(parent="Clinical course")

M_o_inheritance_tbl <- enframe(M_o_inheritance)
M_o_inheritance_tbl <- M_o_inheritance_tbl %>%
  mutate(parent="Mode of inheritance")

Freq_tbl <- enframe(Freq)
Freq_tbl <- Freq_tbl %>%
  mutate(parent="Frequency")

Bld_grp_tbl <- enframe(Bld_grp)
Bld_grp_tbl <- Bld_grp_tbl %>%
  mutate(parent ="Blood group")

HPO_removal_list <-bind_rows(clinical_modifier_tbl,Clinical_course_tbl,M_o_inheritance_tbl,Freq_tbl,Bld_grp_tbl)
length(HPO_removal_list$parent) #Should be 168 +2 for blood group - 170
HPO_removal_list <- HPO_removal_list %>%
  select(value,parent)

write.csv(HPO_removal_list, file="HPO_removal_list_05042019.csv")

#Session_Info
sessionInfo()


#########################################################
#Test out

pheno_gene <- read.delim("ALL_SOURCES_ALL_FREQUENCIES_phenotype_to_genes_28022019_12022019.txt", header = FALSE, sep = "\t")
pheno_gene <- as_tibble(pheno_gene)
colnames(pheno_gene) <- c("HPO_ID","HPO_Name","GeneID","GeneSymbol")
pheno_gene <- pheno_gene[-1,]
#pheno_gene - 483449 to 472515
pheno_gene %>%
  filter(!HPO_ID %in% HPO_removal_list$value) #HPO_ID that is not in the removal_list keep
#
