#SINGLE HPO TERMS

#' Enriched chemicals for a single phenotype term's associated gene list
#'
#' @param HPO Human Phenotype Ontology term character.
#' @param enrich_1S Enrichment one-sided. Default set to true for enrichment, where the one-sided Fisher's exact is set to alternative = "greater."
#' If enrich_1S = FALSE, Fisher's exact is set to alternative = "two.sided" for two-sided.
#' @return Tibble data frame with enrichment.
#' @note HPO terms used should be phenotypic abnormality terms within the HPO release 2019-02-12.
#' @examples
#' perfFishTestHPOSingle("Rickets")
#' @export

perfFishTestHPOSingle <- function(HPO, enrich_1S = TRUE){

#Load as tibble
chem2gene_tbl <- tibble::as_tibble(chem2gene)
pheno2gene_tbl <- tibble::as_tibble(pheno2gene)

#Genes in both datasets
genesannotate <- Reduce(base::intersect, list(pheno2gene$GeneID,chem2gene$GeneID))

#Filter pheno2gene with genes in both datasets
pheno <- pheno2gene_tbl %>%
  dplyr::filter(.data$GeneID %in% genesannotate)

#Change human phenotype term to lower case
HPO_low <- tolower(HPO)
#Name_Given
HPO_Name_Given <- HPO

#Create dataframe of lowercase and HPO terms
pheno_HPO_Cap <- pheno %>%
  dplyr::select(.data$HPO_Name) %>%
  dplyr::distinct() %>%
  dplyr::mutate(lowercase = tolower(.data$HPO_Name))

#Look up proper HPO term
pheno_HPO_Cap <- pheno_HPO_Cap %>%
  dplyr::filter(.data$lowercase == HPO_low)

#If HPO term is not in dataset return that it is not supported
if ( length(pheno_HPO_Cap$HPO_Name) == 0) {
  stop(paste0(HPO_Name_Given," is not supported"))
} else {
  #Assign HPO Name
  HPO <- as.character(pheno_HPO_Cap$HPO_Name)
  print("HPO term accepted")
}


#Enrichment - True or False
if (enrich_1S == TRUE) {
  side <- "greater"
} else if (enrich_1S == FALSE) {
  side <- "two.sided"
} else
  stop("enrich_1S should be TRUE or FALSE")

#Filter by selected HPO term
selectedpheno <-pheno %>%
  dplyr::filter(.data$HPO_Name == HPO)

#Create HPO term genelist
genelist  <-  unique(selectedpheno$GeneID)

#Genelist
gl <- chem2gene_tbl %>%
  dplyr::filter(.data$GeneID %in% genelist)%>%
  dplyr::group_by(.data$X..ChemicalName) %>%
  dplyr::select(.data$X..ChemicalName,.data$GeneID) %>%
  dplyr::distinct() %>%  #remove duplicate genes
  dplyr::summarise(GL_YES = dplyr::n() ) %>%
  dplyr::mutate(TOTAL_GL = length(genelist)) %>%
  dplyr::arrange(dplyr::desc(.data$GL_YES)) %>%
  dplyr::mutate(GL_NO = .data$TOTAL_GL - .data$GL_YES)


#Not Genelist
Not_Genelist <- chem2gene_tbl %>%
  dplyr::filter(!.data$GeneID %in% genelist)

#Not Genelist
nl <- Not_Genelist %>%
  dplyr::group_by(.data$X..ChemicalName) %>%
  dplyr::select(.data$X..ChemicalName,.data$GeneID) %>%
  dplyr::distinct() %>%
  dplyr::summarise(NOT_GL_YES = dplyr::n() ) %>%
  dplyr::mutate(TOTAL_NOT_GL = length(unique(Not_Genelist$GeneID)))%>%
  dplyr::arrange(dplyr::desc(.data$NOT_GL_YES)) %>%
  dplyr::mutate(NOT_GL_NO = .data$TOTAL_NOT_GL - .data$NOT_GL_YES)

#Chemicals that have been filtered out because they do not have any genes in the not gene list
#need to be re-added with a NOT_GL of 0 in nl
chemical_addition <- chem2gene_tbl$X..ChemicalName[!chem2gene_tbl$X..ChemicalName %in% nl$X..ChemicalName]
chemical_addition <- unique(chemical_addition) #just need the term once
chemical_addition_2 <- tibble::enframe(chemical_addition)
chemical_addition_2 <-   chemical_addition_2[,-1]
colnames(chemical_addition_2) <- "X..ChemicalName"

chemical_addition_3 <- chemical_addition_2 %>%
  dplyr::mutate(NOT_GL_YES = 0,
                TOTAL_NOT_GL = length(unique(Not_Genelist$GeneID)),
                NOT_GL_NO = .data$TOTAL_NOT_GL - .data$NOT_GL_YES)

#Bind the nl(Not_Genelist) and chemical_addition_3
nl_2 <- dplyr::bind_rows(nl,chemical_addition_3 )

#Join dataframes
fl <- gl %>%
  dplyr::left_join(nl_2, by = "X..ChemicalName")


#Fisher's exact test function
#Attribution: https://stackoverflow.com/questions/28966840/fishers-test-stat-on-data-frame-using-dplyrmutate
fishertest <- function(a,b,c,d){
  m <- matrix(c(a,b,c,d),nrow=2)
  p <- stats::fisher.test(m,alternative = side)$p.value
}

#Fisher's exact test p_value
fl2 <- fl %>%
  dplyr::group_by(.data$X..ChemicalName)%>%
  dplyr::mutate(p_value = fishertest(.data$GL_YES,.data$NOT_GL_YES,.data$GL_NO,.data$NOT_GL_NO))%>%
  dplyr::arrange(dplyr::desc(.data$p_value))%>%
  dplyr::ungroup()

#p value adjustment
fl2 <- fl2 %>%
  dplyr::mutate(bonf = stats::p.adjust(.data$p_value,"bonferroni")) %>%
  dplyr::mutate(FDR = stats::p.adjust(.data$p_value,"fdr")) %>%
  dplyr::arrange(.data$p_value)
fl3 <- fl2
#Change first column name
colnames(fl3)[1] <- "Chemical_Name"
fl3
}




###############

HPOlistfunction <- function(i){
  chem2gene_tbl <- tibble::as_tibble(chem2gene)
  pheno2gene_tbl <- tibble::as_tibble(pheno2gene)
  #Genes in both datasets
  genesannotate <- Reduce(base::intersect, list(pheno2gene$GeneID,chem2gene$GeneID))

  pheno <- pheno2gene_tbl %>%
    dplyr::filter(.data$GeneID %in% genesannotate)

  #Change human phenotype term to lower case
  HPO_low <- tolower(i)
  #Name_Given
  HPO_Name_Given <- i

  #Create dataframe of lowercase and HPO terms
  pheno_HPO_Cap <- pheno %>%
    dplyr::select(.data$HPO_Name) %>%
    dplyr::distinct() %>%
    dplyr::mutate(lowercase = tolower(.data$HPO_Name))

  #Look up proper HPO term
  pheno_HPO_Cap <- pheno_HPO_Cap %>%
    dplyr::filter(.data$lowercase == HPO_low)

   #If HPO term is not in dataset return that it is not supported
  if ( length(pheno_HPO_Cap$HPO_Name) == 0) {
    stop(paste0(HPO_Name_Given," is not supported"))
  } else {
    #Assign HPO Name
    HPO <- as.character(pheno_HPO_Cap$HPO_Name)
    print("HPO term accepted")
  }

  #Filter by selected HPO term
  selectedpheno <- pheno %>%
    dplyr::filter(.data$HPO_Name == HPO)
  #Create HPO term genelist
  genelist <- unique(selectedpheno$GeneID)

}


#' Chemical enrichment for multiple phenotype terms
#'
#' @param HPO List of Human Phenotype Ontology terms characters
#' @param enrich_1S Enrichment one-sided. Default set to true for enrichment, where the one-sided Fisher's exact is set to alternative = "greater."
#' If enrich_1S = FALSE, Fisher's exact is set to alternative = "two.sided" for two.sided.
#' @return Tibble data frame with enrichment.
#' @examples
#' perfFishTestHPOMultiple(list("Increased circulating ACTH level",
#' "Androgen insufficiency",
#' "Decreased circulating aldosterone level"))
#' @note Chemical enrichment for multiple HPO terms is a simplistic exploratory function, where we take an additive approach.
#' HPO terms used should be phenotypic abnormality terms within the HPO release 2019-02-12.
#' @export

perfFishTestHPOMultiple <- function(HPO, enrich_1S = TRUE){

  #Enrichment - True or False
  if (enrich_1S == TRUE) {
    side <- "greater"
  } else if (enrich_1S == FALSE) {
    side <- "two.sided"
  } else
    stop("enrich_1S should be TRUE or FALSE")

  #Create a genelist
  agenelist <- lapply(HPO, HPOlistfunction)
  unlist_agenelist <- unlist(agenelist) #deduplicate
  genelist <-  unique(unlist_agenelist)

  #Genes contained in both datasets
  chem2gene_tbl <- tibble::as_tibble(chem2gene)
  genesannotate <- Reduce(base::intersect, list(pheno2gene$GeneID,chem2gene$GeneID))

   #Genelist
  gl <- chem2gene_tbl %>%
    dplyr::filter(.data$GeneID %in% genelist)%>%
    dplyr::group_by(.data$X..ChemicalName) %>%
    dplyr::select(.data$X..ChemicalName,.data$GeneID) %>%
    dplyr::distinct() %>%  #remove duplicate genes
    dplyr::summarise(GL_YES = dplyr::n() ) %>%
    dplyr::mutate(TOTAL_GL = length(genelist)) %>%
    dplyr::arrange(dplyr::desc(.data$GL_YES)) %>%
    dplyr::mutate(GL_NO = .data$TOTAL_GL - .data$GL_YES)

  #Not Genelist
  Not_Genelist <- chem2gene_tbl %>%
    dplyr::filter(!.data$GeneID %in% genelist)

  #Not Genelist
  nl <- Not_Genelist %>%
    dplyr::group_by(.data$X..ChemicalName) %>%
    dplyr::select(.data$X..ChemicalName,.data$GeneID) %>%
    dplyr::distinct() %>%
    dplyr::summarise(NOT_GL_YES = dplyr::n() ) %>%
    dplyr::mutate(TOTAL_NOT_GL = length(unique(Not_Genelist$GeneID)))%>%
    dplyr::arrange(dplyr::desc(.data$NOT_GL_YES)) %>%
    dplyr::mutate(NOT_GL_NO = .data$TOTAL_NOT_GL - .data$NOT_GL_YES)

  #Chemicals that have been filtered out because they do not have any genes in the not gene list
  #need to be re-added with a NOT_GL of 0 in nl
  chemical_addition <- chem2gene_tbl$X..ChemicalName[!chem2gene_tbl$X..ChemicalName %in% nl$X..ChemicalName]
  chemical_addition <- unique(chemical_addition) #just need the term once
  chemical_addition_2 <- tibble::enframe(chemical_addition)
  chemical_addition_2 <-   chemical_addition_2[,-1]
  colnames(chemical_addition_2) <- "X..ChemicalName"

  chemical_addition_3 <- chemical_addition_2 %>%
    dplyr::mutate(NOT_GL_YES = 0,
                  TOTAL_NOT_GL = length(unique(Not_Genelist$GeneID)),
                  NOT_GL_NO = .data$TOTAL_NOT_GL - .data$NOT_GL_YES)

  #Bind the nl(Not_Genelist) and chemical_addition_3
  nl_2 <- dplyr::bind_rows(nl,chemical_addition_3 )

  #Join dataframes
  fl <- gl %>%
    dplyr::left_join(nl_2, by = "X..ChemicalName")


  #Fishers exact test
  #Attribution: https://stackoverflow.com/questions/28966840/fishers-test-stat-on-data-frame-using-dplyrmutate
  fishertest <- function(a,b,c,d){
    m <- matrix(c(a,b,c,d),nrow=2)
    p <- stats::fisher.test(m,alternative = side)$p.value
  }

  #Fisher's exact test p_value
  fl2 <- fl %>%
    dplyr::group_by(.data$X..ChemicalName)%>%
    dplyr::mutate(p_value = fishertest(.data$GL_YES,.data$NOT_GL_YES,.data$GL_NO,.data$NOT_GL_NO))%>%
    dplyr::arrange(dplyr::desc(.data$p_value))%>%
    dplyr::ungroup()

  #P-value adjustment
  fl2 <- fl2 %>%
    dplyr::mutate(bonf = stats::p.adjust(.data$p_value,"bonferroni")) %>%
    dplyr::mutate(FDR = stats::p.adjust(.data$p_value,"fdr")) %>%
    dplyr::arrange(.data$p_value)
  fl3 <- fl2
  #Change first column name
  colnames(fl3)[1] <- "Chemical_Name"
  fl3
}


