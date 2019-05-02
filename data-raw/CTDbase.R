setwd('data-raw/')
library(tibble)
library(dplyr)
library(readr)
library(tidyr)
library(magrittr)

#
chem2gene <- read.csv("CTD_chem_gene_ixns_11022019_05022019.csv")
head(chem2gene,5)
#remove first empty row
chem2gene <- chem2gene[-1,]
#chem2gene
head(chem2gene,5)
#summary
summary(chem2gene) #There is a NA within GeneSymbol and OrganismID. Empty fields in Organism, GeneForms

#tibble
chem2gene_tb <- as_tibble(chem2gene)
chem2gene_tb
#####Look up the NA
chem2gene_tb %>%
  filter(is.na(GeneSymbol))
#NA Genesymbol with GeneID 45338 for chemical entinostat
chem2gene_tb %>%
  filter(is.na(OrganismID)) #No organism - 23,811

#unique Gene Symbols and GeneID
length(unique(chem2gene$GeneSymbol)) #47588
length(unique(chem2gene$GeneID)) #47588
nrow(chem2gene) #1749548



#Dataset cleaning #####################################################
chem2gene_tb$GeneID

#Homo sapien gene info - NCBI (ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia)
human_geneinfo <- read_tsv(file = "Homo_sapiens_19032019_19032019.gene_info")
human_geneinfo
summary(human_geneinfo)
human_geneinfo  <- human_geneinfo %>%
  dplyr::select(GeneID,Symbol)

human_geneinfo #61,588
# distinct(human_geneinfo) - no duplicate rows


#########################################################################################
#CTD Dataset
summary(chem2gene_tb$Organism)
# chem2gene_tb %>%
#   filter(Organism == "Rattus rattus")

#CTD - filtered by selected organisms
chem2gene_tb_org  <- chem2gene_tb %>%
  filter(Organism %in% c("Homo sapiens",
                         "Mus musculus",
                         "Rattus norvegicus",
                         "Rattus rattus",
                         "Danio rerio",
                         "Macaca fascicularis",
                         "Macaca mulatta"))
#1,624,760
summary(chem2gene_tb_org$Organism)
#CTD- filtered by geneIDs belonging to humans - 1,539,606
Human_IDs <- chem2gene_tb_org %>%
  filter(GeneID %in% human_geneinfo$GeneID)

#Orthologs - filter
#Filter out humans and geneIDs not belonging to humans - 81,699
ORTHOLOGS <- chem2gene_tb_org %>%
  filter(!GeneID %in% human_geneinfo$GeneID) %>%
  filter(!Organism == "Homo sapiens")

##Does the Organism of the Gene ID match the OrganismID? Remove if not.

#Mammalia gene info from NCBI /gene/DATA/GENE_INFO/Mammalia/
mam_genes <-read.delim("All_Mammalia_15032019_15032019.gene_info",header=TRUE)
mam_genes <- as_tibble(mam_genes)
mam_genes
summary(mam_genes) #Summary - okay
mam_genes <- mam_genes %>%
  dplyr::select(X.tax_id,GeneID,Symbol,Synonyms,dbXrefs,description,type_of_gene,Symbol_from_nomenclature_authority,Modification_date)
mam_genes
# mam_genes %>%
#   distinct() - no change - 3,843,545

#Animal ids are present
#https://www.ncbi.nlm.nih.gov/Taxonomy
#"Mus musculus" - 10090
#"Rattus norvegicus" - 10116
#"Rattus rattus" -  10117
#"Macaca fascicularis" - 9541
#"Macaca mulatta" - 9544
# 3,843,545 -> 184,157
mam_genes %>%
  filter(X.tax_id %in% c("10090","10116","10117","9541","9544")) %>%
  select(X.tax_id) %>%
  distinct() #All present

# Org_Gene_test <- filter(mam_genes, GeneID == 11947)
# tax_test<- Org_Gene_test$X.tax_id

#Lookup Function - to find the organism ID of the gene ID
LookupFunc <- function(a){
  Org_Gene <- filter(mam_genes, GeneID == a)
  tax <- Org_Gene$X.tax_id
  if (length(tax) == 0) {
    tax <- NA
  }
  tax
}

#Filter the orthologs by mammalia - no Danio rerio
#81,699 to 58,328
ORTHOLOGS_MAM <- ORTHOLOGS %>%
  filter(Organism %in% c("Mus musculus",
                         "Rattus norvegicus",
                         "Rattus rattus",
                         "Macaca fascicularis",
                         "Macaca mulatta"))
#58,328
# ORTHOLOGS %>%
#   filter(OrganismID %in% c("10090","10116","10117","9541","9544"))

#Reduce size by looking at GeneSymbol, GeneID, Organism and OrganismID - 8,666
ORTHOLOGS_MAM_DISTINCT <- ORTHOLOGS_MAM %>%
  dplyr::select(GeneSymbol,GeneID,Organism,OrganismID) %>%
  distinct()

#Lookup the OrganismID of the GeneID - new column GeneID_OrganismID
Orth_Mam_Dis_ID <- ORTHOLOGS_MAM_DISTINCT %>%
  mutate(RowNo = 1:n()) %>%
  group_by(RowNo) %>%
  mutate(GeneID_OrganismID = LookupFunc(GeneID)) %>%
  ungroup()

# Orth_Mam_Dis_ID_row <- ORTHOLOGS_MAM_DISTINCT %>%
#   mutate(RowNo = 1:n())
# tail(Orth_Mam_Dis_ID_row)

#GeneID and GeneID_OrganismID distinct dataframe
Orth_Mam_Dis_ID_DISTINCT <- Orth_Mam_Dis_ID %>%
  dplyr::select(GeneID,GeneID_OrganismID) %>%
  distinct()

#NA's in GeneID_OrganismID - as GeneID belongs to an Organism that is not Mammalia
 Orth_Mam_Dis_ID_DISTINCT %>%
  filter(is.na(GeneID_OrganismID))

#Join Mammalia Ortholog dataframe with our distinct ID dataframe by GeneID
#Left join - "return all rows from x, and all columns from x and y.
#Rows in x with no match in y will have NA values in the new columns.
#If there are multiple matches between x and y, all combinations of the matches are returned."
#https://dplyr.tidyverse.org/reference/join.html
Orthologs_MAM_2 <- ORTHOLOGS_MAM %>%
  left_join(Orth_Mam_Dis_ID_DISTINCT, by = "GeneID")

#Compare OrganismID with GeneID OrganismID
Orthologs_MAM_Compare <- Orthologs_MAM_2 %>%
  mutate(Compare = OrganismID == GeneID_OrganismID)

summary(Orthologs_MAM_Compare$Compare)
#
Orthologs_MAM_Compare %>%
  filter(Compare == FALSE) #10,831

Orthologs_MAM_Compare %>%
  filter(is.na(Compare)) #2,409

#Filter by Compare is TRUE to find those who match both OrganismID and GeneID_OrganismID
#45,088
Orthologs_MAM_Compare_TRUE <- Orthologs_MAM_Compare %>%
  filter(Compare == TRUE)

#Quick summary check
summary(Orthologs_MAM_Compare_TRUE)


#Danio_rerio gene info - Danio rerio - 7955
ZEBFISH_genes <-read.delim("Danio_rerio_15032019_15032019.gene_info",header=TRUE)
ZEBFISH_genes <- as_tibble(ZEBFISH_genes)
ZEBFISH_genes #43,988
summary(ZEBFISH_genes)
#Select
ZEBFISH_genes <- ZEBFISH_genes %>%
  dplyr::select(X.tax_id,GeneID,Symbol,Synonyms,dbXrefs,description,type_of_gene,Symbol_from_nomenclature_authority,Modification_date)

# ZEBFISH_genes %>%
#  distinct() - 43,988

#Filter Orthologs by Danio rerio orthologs
ORTHOLOGS_ZEBFISH <- ORTHOLOGS %>%
  filter(Organism == "Danio rerio")

summary(ORTHOLOGS_ZEBFISH)

#Reduce and distinct - GeneSymbol, GeneID, Organism and OrganismID
ORTHOLOGS_ZEBFISH_DISTINCT <- ORTHOLOGS_ZEBFISH %>%
  dplyr::select(GeneSymbol,GeneID,Organism,OrganismID) %>%
  distinct()
#5,776

#ZebraFish Look-up
LookupFunc_FISH <- function(a){
  Org_Gene <- filter(ZEBFISH_genes, GeneID == a)
  tax <- Org_Gene$X.tax_id
  if (length(tax) == 0) {
    tax <- NA
  }
  tax
}

#Lookup the Organism ID that the Gene ID belongs to
Orth_ZF_Dis_ID <- ORTHOLOGS_ZEBFISH_DISTINCT %>%
  mutate(RowNo = 1:n()) %>%
  group_by(RowNo) %>%
  mutate(GeneID_OrganismID = LookupFunc_FISH(GeneID))%>%
  ungroup()

tail(Orth_ZF_Dis_ID)
#GeneID and GeneID_OrganismID dataframe distinct
Orth_ZF_Dis_ID_Dis <- Orth_ZF_Dis_ID %>%
  dplyr::select(GeneID,GeneID_OrganismID) %>%
  distinct()
#Summary
summary(Orth_ZF_Dis_ID_Dis)

#GeneID_OrganismID is not Zebra Fish
Orth_ZF_Dis_ID_Dis %>%
  filter(is.na(GeneID_OrganismID))

#Left join by GeneID to match to OrganismID_GeneID
Orthologs_ZEBFISH_2 <- ORTHOLOGS_ZEBFISH  %>%
  left_join(Orth_ZF_Dis_ID_Dis, by = "GeneID")

#NA - Because the organism is not Danio rerio
Orthologs_ZEBFISH_2 %>%
  filter(is.na(GeneID_OrganismID)) %>%
  dplyr::select(GeneID,GeneID_OrganismID)

#Comparison of OrganismID and GeneID_OrganismID
Orthologs_ZEBFISH_COMPARE <- Orthologs_ZEBFISH_2 %>%
  mutate(Compare = OrganismID == GeneID_OrganismID)

summary(Orthologs_ZEBFISH_COMPARE$Compare)
#No GeneID_OrganismID
Orthologs_ZEBFISH_COMPARE %>%
  filter(is.na(Compare))

#Select for Comparison is TRUE - 22,842
Orthologs_ZEBFISH_COMPARE_TRUE <- Orthologs_ZEBFISH_COMPARE %>%
  filter(Compare == TRUE)

#Join the Mammalian and Fish Orthologs
Orthologs_MAMFISH <- dplyr::bind_rows(Orthologs_MAM_Compare_TRUE,Orthologs_ZEBFISH_COMPARE_TRUE)
Orthologs_MAMFISH #67,930
summary(Orthologs_MAMFISH)

#Load annotationTools and homologene file
library(annotationTools)
homologene<-read.delim("homologene_11032019_06052014.data",header=FALSE)
summary(homologene)
head(homologene,3)
#getHOMOLOG(geneid, targetspecies, homol)

#Greatest number of orthologs - 5
unique(Orthologs_MAMFISH$GeneID)[1]
#give list
entire_ortho <- getHOMOLOG(unique(Orthologs_MAMFISH$GeneID),9606,homologene)
#list
entire_ortho
#give list names of the gene IDs
names(entire_ortho) <- unique(Orthologs_MAMFISH$GeneID)
entire_ortho
#make table
entire_ortho_tb <- enframe(entire_ortho)
#Look at greatest number of orthologs - 5
entire_ortho_tb %>%
  unnest() %>%
  group_by(name) %>%
  summarise(ort = n()) %>%
  arrange(desc(ort))

#OrthFunc - get homolog, if length is 2 as in there is more than one ortholog assign 0.2, if 3 assign 0.3 etc.
OrthFunc_2 <- function(a,b){
  Or <- getHOMOLOG(a,b,homologene)
  Or <- unlist(Or)
  if (length(Or) == 2) {
    Or <- 0.2
  } else if (length(Or) == 3) {
    Or <- 0.3
  } else if (length(Or) == 4) {
    Or <- 0.4
  } else if (length(Or) == 5) {
    Or <- 0.5
  }
  Or
}

#Lookup the Orthologs of the GeneID in humans 9606
#warnings refer to no homolog - getHOMOLOG(226143,9606,homologene)
new_ID_annotationtools <- Orthologs_MAMFISH %>%
  mutate(RowNo = 1:n()) %>%
  group_by(RowNo) %>%
  mutate(new_ID = OrthFunc_2(GeneID,9606)) %>%
  ungroup()

tail(new_ID_annotationtools)
#new_ID_annotationtools - 67,930

#no homolog - 33,444
new_ID_annotationtools %>%
  filter(is.na(new_ID))

summary(new_ID_annotationtools) #min 0 from the 0.2 - 0.5



#Has one ortholog only -
new_ID_one_ortholog <- new_ID_annotationtools %>%
  filter(!is.na(new_ID)) %>% #remove NA - 34,486
  filter(new_ID > 0.6) #remove more than one ortholog - 33,885
#34486 - 33885 = 601

summary(new_ID_one_ortholog)
#More than one Ortholog - 601
More1_Ortholog <- new_ID_annotationtools %>%
  filter(!is.na(new_ID)) %>%
  filter(new_ID < 1) %>%
  filter(new_ID > 0.05)
More1_Ortholog
summary(More1_Ortholog)

#OrthFunc_3, select ortholog i
OrthFunc_3 <- function(a,b,i){
  Or <- getHOMOLOG(a,b,homologene)
  Or <- unlist(Or)
  Or[i]
}

#5 orthologs is the greatest number
More1_Ortholog %>%
  arrange(desc(new_ID))

#New columns with ids
More1_Ortholog_2 <- More1_Ortholog  %>%
  mutate(RowNo = 1:n()) %>%
  group_by(RowNo) %>%
  mutate(new_ID_1 = OrthFunc_3(GeneID,9606,1))%>%
  mutate(new_ID_2 = OrthFunc_3(GeneID,9606,2))%>%
  mutate(new_ID_3 = OrthFunc_3(GeneID,9606,3))%>%
  mutate(new_ID_4 = OrthFunc_3(GeneID,9606,4))%>%
  mutate(new_ID_5 = OrthFunc_3(GeneID,9606,5))%>%
  ungroup()


#New_ID_1:5
More1_Ortholog_2 %>%
  dplyr::select(15:20)%>%
  arrange(desc(new_ID))

#Gather the new_ID's - remove NAs
More1_Ortholog_COMPLETE <- More1_Ortholog_2 %>%
  gather(key = Category, value = ID, 16:20, na.rm = TRUE)

##Check - new IDs
More1_Ortholog_COMPLETE %>%
  dplyr::group_by(RowNo) %>%
  dplyr::filter(ChemicalID == "C028474")

#Remove compare
new_ID_one_ortholog_2 <- new_ID_one_ortholog %>%
  dplyr::select(-Compare)
#Rename new_ID to New_GeneID
names(new_ID_one_ortholog_2)[names(new_ID_one_ortholog_2) == "new_ID"] <- "New_GeneID"

#Remove Category, Compare, new_ID
More1_Ortholog_COMPLETE_2 <- More1_Ortholog_COMPLETE %>%
  dplyr::select(-Compare,-Category,-new_ID)

#Rename ID to New_GeneID
names(More1_Ortholog_COMPLETE_2 )[names(More1_Ortholog_COMPLETE_2 ) == "ID"] <- "New_GeneID"
More1_Ortholog_COMPLETE_2


########Bind rows for one ortholog and more than one ortholog
#33,885 + 1,544 = 35429
Orthologs_COMPLETE <- bind_rows(new_ID_one_ortholog_2,More1_Ortholog_COMPLETE_2)
Orthologs_COMPLETE

# write.csv(Orthologs_COMPLETE, "orthologs_complete.csv")

######################### UPDATE GENE IDs OF ORTHOLOGS #################
#gene history from ftp.ncbi.nih.gov/gene/DATA
gene_history <-read.delim("gene_history_15032019_15032019",header=TRUE)
head(gene_history)
summary(gene_history)

#human gene history
gene_history_human <- gene_history %>%
  filter(X.tax_id == 9606)

summary(gene_history_human)
length(unique(gene_history_human$Discontinued_GeneID))

#GeneIDs that need updating
Orthologs_COMPLETE_GID_UPDATE <- Orthologs_COMPLETE %>%
  filter(New_GeneID %in% gene_history_human$Discontinued_GeneID)
#32 rows

#Rename New_GeneID to Discontinued_GeneID
names(Orthologs_COMPLETE_GID_UPDATE)[names(Orthologs_COMPLETE_GID_UPDATE) == "New_GeneID"] <- "Discontinued_GeneID"


gene_history_human <- as_tibble(gene_history_human)
gene_history_human <- gene_history_human %>%
  select(GeneID,Discontinued_GeneID)
distinct(gene_history_human) #163,013 - none
#Rename_column name - GeneID here to New_GeneID
gene_history_human
names(gene_history_human)[names(gene_history_human) == "GeneID"] <- "New_GeneID"

#Check
gene_history_human %>%
  filter(Discontinued_GeneID == 100506627)

#Left_join by Discontinued_Gene to include updated New_GeneID to dataframe
Orthologs_COMPLETE_GID_UPDATE_2 <- Orthologs_COMPLETE_GID_UPDATE %>%
  left_join(gene_history_human, by = "Discontinued_GeneID")

summary(Orthologs_COMPLETE_GID_UPDATE_2) #Every GeneID has a new GeneID

#Remove Discontinued_GeneID
Orthologs_COMPLETE_GID_UPDATE_2 <- Orthologs_COMPLETE_GID_UPDATE_2 %>%
  select(-Discontinued_GeneID)

####Recombine new updated GeneIDs with dataframe where GeneID's did not need updating
#Orthologs_COMPLETE - 35,429
##GeneIDs that don't need updating
Orthologs_COMPLETE_GID_OK <- Orthologs_COMPLETE %>%
  filter(!New_GeneID %in% gene_history_human$Discontinued_GeneID)

#35,397 + 32 = 35,429
Orthologs_COMPLETE_GID_OK #35,397
Orthologs_COMPLETE_GID_UPDATE_2 #32

#as.numeric for row binding
#When changing from factor, needs to be as.numeric(as.character())
#Orthologs_COMPLETE_GID_UPDATE_2$PubMedIDs <- as.integer(Orthologs_COMPLETE_GID_UPDATE_2$PubMedIDs)
Orthologs_COMPLETE_GID_UPDATE_2$New_GeneID <- as.numeric(as.character(Orthologs_COMPLETE_GID_UPDATE_2$New_GeneID))
Orthologs_COMPLETE_GID_UPDATE_2
#Bind rows
Orthologs_COMPLETE_GID_UPDATED <- dplyr::bind_rows(Orthologs_COMPLETE_GID_OK,Orthologs_COMPLETE_GID_UPDATE_2) #35,429

######################### NEW GENE SYMBOL ###########################################################################

Orthologs_COMPLETE_GID_UPDATED

#Assign to symbol_tb
symbol_tb <- human_geneinfo
#Rename
colnames(symbol_tb) <- c("New_GeneID","New_GeneSymbol")

#As.integer - not needed
# symbol_tb$New_GeneID <- as.integer(symbol_tb$New_GeneID)

#Join by new GeneID to get new GeneSymbol
Orthologs_COMPLETE_NG <- Orthologs_COMPLETE_GID_UPDATED %>%
  left_join(symbol_tb, by = "New_GeneID")

#Check
summary(Orthologs_COMPLETE_NG)
Orthologs_COMPLETE_NG  %>%
  filter(is.na(New_GeneSymbol))


#Remove old GeneID, GeneSymbol, GeneID_OrganismID and RowNo
Orthologs_COMPLETE_NG <- Orthologs_COMPLETE_NG %>%
  select(-GeneID,-GeneSymbol,-GeneID_OrganismID,-RowNo)

#Rename New_GeneID to GeneID and New_GeneSymbol to GeneSymbol
names(Orthologs_COMPLETE_NG)[names(Orthologs_COMPLETE_NG) == "New_GeneID"] <- "GeneID"
names(Orthologs_COMPLETE_NG)[names(Orthologs_COMPLETE_NG) == "New_GeneSymbol"] <- "GeneSymbol"

#Bind orthologs to chem2gene dataset
Human_IDs #1,539,606
Human_IDs$GeneSymbol <- as.character(Human_IDs$GeneSymbol) #To character
Orthologs_COMPLETE_NG #35,429
chem2gene_fin <- dplyr::bind_rows(Orthologs_COMPLETE_NG,Human_IDs)
chem2gene_fin #1,539,606 + 35,429 = 1575032 #1,575,035

#summary
summary(chem2gene_fin)

#All ID's are human ID's
chem2gene_fin %>%
  filter(!GeneID %in% human_geneinfo$GeneID)

#Previous dataframe
nrow(chem2gene) #1749548
chem2gene <- chem2gene_fin
nrow(chem2gene) #1575035
#Unique GeneID's - 23754
length(unique(chem2gene$GeneID))
#write to a csv to load for use in SYM.R
write.csv(chem2gene,"chem2gene_usedata_05042019.csv")
#Use this data
usethis::use_data(chem2gene,overwrite=TRUE)

#Session_Info
sessionInfo()
