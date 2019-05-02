setwd('data-raw/')
library(dplyr)
library(tidyr)
library(stringr)
library(magrittr)
library(tibble)
#read in file and replace "" with NA
sym <- read.csv('CTD_chemicals_11022019_05022019.csv',na.strings = "")
#look at head and summary
head(sym,5)
summary(sym)
#need to text edit the file
sym_tbl <- as_tibble(sym)
#remove first row as NA
sym_tbl  <- sym_tbl[-1,]
#
summary(sym_tbl) #Summary check is fine

#Dextrin has two synonyms caloreen and dextrin
sym_tbl %>%
  filter(Synonyms == "Dextrin")

sym_tbl %>%
  filter(Synonyms == "Dextrin") %>%
  select(ChemicalID, ParentIDs)

##### Load cleaned chem2gene
 chem2gene <- read.csv("chem2gene_usedata_05042019.csv")
 head(chem2gene,5)
 summary(chem2gene) #Summary check is fine


###################################

sym_tbl
###REMOVE THE MESH: FROM CHEMICALID
sym_tbl_1 <- sym_tbl %>%
  mutate(ChemID = str_remove(ChemicalID,"MESH:")) %>%
  dplyr::select(X..ChemicalName,Synonyms,ChemID) %>%
  rename(ChemicalID = ChemID)

#sym_tbl_1$Synonyms <- as.character(sym_tbl_1$Synonyms)
sym_tbl_1
#find max number of synonyms that a chemical has, to properly separate the dataset
sym_tbl_1 %>%
  mutate(count = str_count(Synonyms, "\\|")) %>%
  group_by(ChemicalID) %>%
  arrange(desc(count))
#78 is our biggest number of synonyms for Acetylcysteine - 11022019
sym_tbl_2 <- sym_tbl_1 %>%
  mutate(count = str_count(Synonyms, "\\|")) %>%
   separate(Synonyms,paste0("Synonym",1:78) , sep = "\\|", extra="merge", fill="right") %>%
   rename(Synonym0 = X..ChemicalName)

#Synonym0 as a character - to remove warning of "Warning message:
#attributes are not identical across measure variables;
#they will be dropped
sym_tbl_2$Synonym0 <- as.character(sym_tbl_2$Synonym0)

sym_tbl_3 <- sym_tbl_2 %>%
  gather(key = synonym_no ,value = synonym,Synonym0:Synonym78 ,na.rm = TRUE) %>%
  group_by(ChemicalID) %>%
  #mutate(number = str_extract(synonym_no,"\\d"))%>%
  arrange(desc(ChemicalID))

###FILTER CHEMICAL VOCABULARY BY CHEMICALS WITH GENE INTERACTIONS - because some chemicals do not have gene interactions
sym_tbl_4 <- sym_tbl_3 %>%
  ungroup()

#just needed to ungroup ! - 430512 rows to 63,868 rows
#unselect count
sym_tbl_5 <- sym_tbl_4 %>%
  filter(ChemicalID %in% unique(chem2gene$ChemicalID)) %>%
  dplyr::select(-count)


############Check for duplicate synonyms
sym_tbl_5 %>%
  group_by(synonym) %>%
  summarise(Number = n()) %>%
  arrange(desc(Number))

# Tetrachlorodibenzodioxin + Polychlorinated Dibenzodioxins
# Duplicates of
# 2,3,7,8-Tetrachlorodibenzo-p-dioxin, Dibenzo(b,e)(1,4)dioxin, 2,3,7,8-tetrachloro-, TCDD, Tetrachlorodibenzodioxin

sym_tbl_5 %>%
  filter(synonym == "2,3,7,8-Tetrachlorodibenzo-p-dioxin")

sym_tbl_5 %>%
  filter(synonym == "Dibenzo(b,e)(1,4)dioxin, 2,3,7,8-tetrachloro-")

sym_tbl_5 %>%
  filter(synonym == "TCDD")

sym_tbl_5 %>%
  filter(synonym == "Tetrachlorodibenzodioxin")

#
PCDD <- sym_tbl_5 %>%
  filter(ChemicalID == "D000072317")

View(PCDD)

chem2gene %>%
  filter(ChemicalID == "D013749")

chem2gene %>%
  filter(ChemicalID == "D000072317")


#################################################################################
#Filter out duplicate synonyms in D000072317 - Polychlorinated Dibenzodioxins
# sym_tbl_5 - 63,868
# sym_tbl_6 - should be 63,864
sym_removal <- sym_tbl_5 %>%
  filter(ChemicalID == "D000072317" & c(synonym == "2,3,7,8-Tetrachlorodibenzo-p-dioxin"  |
                                          synonym == "Dibenzo(b,e)(1,4)dioxin, 2,3,7,8-tetrachloro-"  |
                                          synonym == "TCDD"  |
                                          synonym == "Tetrachlorodibenzodioxin"))

sym_removal

#anti_join - to remove the rows
sym_tbl_6 <- anti_join(sym_tbl_5, sym_removal, by = c("ChemicalID", "synonym_no", "synonym"))
#63,864
############Re-check for duplicate synonyms
sym_tbl_6 %>%
  group_by(synonym) %>%
  summarise(Number = n()) %>%
  arrange(desc(Number))
#No duplicate synonyms

sym_tbl_6 %>%
  filter(ChemicalID == "D013749")

sym_tbl_6 %>%
  filter(ChemicalID == "D000072317")

# Change synonyms to all lower case
sym_tbl_6 <- sym_tbl_6 %>%
  mutate(synonym = tolower(synonym))

tail(sym_tbl_6)
summary(sym_tbl_6)

#Re-check for duplicate synonyms
sym_tbl_6 %>%
  group_by(synonym) %>%
  summarise(Number = n()) %>%
  arrange(desc(Number))

#Remove duplicates - 63864 to 63818
sym_tbl_7 <- sym_tbl_6 %>%
  select(ChemicalID, synonym) %>%
  distinct()

#No duplicate synonyms
sym_tbl_7 %>%
  group_by(synonym) %>%
  summarise(Number = n()) %>%
  arrange(desc(Number))

#
sym_tbl <- as.data.frame(sym_tbl_7)
summary(sym_tbl) #Summary - Looks fine
glimpse(sym_tbl)
nrow(sym_tbl) #63818
head(sym_tbl)

#use-data
usethis::use_data(sym_tbl, overwrite = TRUE)

#Session_Info
sessionInfo()
