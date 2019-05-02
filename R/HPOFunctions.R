#SINGLE HPO TERMS

#' Enriched Chemicals for a single phenotype term's associated gene list
#'
#' @param HPO Human Phenotype Ontology term
#' @param enrich_1S Enrichment one-sided. Automatically set to true for enrichment, where the one-sided Fisher's exact is set to alternative = 'greater.'
#' If enrich_1S = FALSE, Fisher's exact is set to alternative = 'two.sided' for two-sided.
#' @return Dataframe with enrichment.
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

#tail(pheno_HPO_Cap)

#Look up proper HPO term
pheno_HPO_Cap <- pheno_HPO_Cap %>%
  dplyr::filter(.data$lowercase == HPO_low)

#If HPO term is not in dataset return that it is not supported
if ( length(pheno_HPO_Cap$HPO_Name) == 0) {
  stop(paste0(HPO_Name_Given," is not supported"))
} else {
  #Assign HPO Name
  HPO <- as.character(pheno_HPO_Cap$HPO_Name)
  print(" ")
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
#Filter by homosapiens removed (  #filter(Organism == "Homo sapiens") %>%) and genes contained in both datasets
# chems <- chem2gene_tbl %>%
#   dplyr::filter(.data$GeneID %in% genesannotate)

#Create HPO term genelist
genelist  <-  unique(selectedpheno$GeneID)

#Genelist
gl <- chem2gene_tbl %>%
  dplyr::filter(.data$GeneID %in% genelist)%>%
  dplyr::group_by(.data$X..ChemicalName) %>%
  dplyr::select(.data$X..ChemicalName,.data$ChemicalID,.data$GeneID) %>%
  dplyr::distinct() %>%  #remove duplicate genes
  dplyr::summarise(GL_YES = dplyr::n() ) %>%
  dplyr::mutate(TOTAL_GL = length(genelist)) %>%
  dplyr::arrange(dplyr::desc(.data$GL_YES)) %>%
  dplyr::mutate(GL_NO = .data$TOTAL_GL - .data$GL_YES)

#Universe list
# filter(Organism == "Homo sapiens") %>% removed and Organism no longer selected
#filter(GeneID %in% genesannotate ) %>% removed

ul <- chem2gene_tbl %>%
  dplyr::group_by(.data$X..ChemicalName) %>%
  dplyr::select(.data$X..ChemicalName,.data$ChemicalID,.data$GeneID) %>%
  dplyr::distinct() %>%
  dplyr::summarise(UL_YES = dplyr::n() ) %>%
  dplyr::mutate(TOTAL_UL = length(unique(chem2gene_tbl$GeneID)))%>%
  dplyr::arrange(dplyr::desc(.data$UL_YES)) %>%
  dplyr::mutate(UL_NO = .data$TOTAL_UL - .data$UL_YES)

#Join dataframes
fl <- gl %>%
  dplyr::left_join(ul, by = "X..ChemicalName")


#Fisher's exact test function
fishertest <- function(a,b,c,d){
  m <- matrix(c(a,b,c,d),nrow=2)
  p <- stats::fisher.test(m,alternative = side)$p.value
}

#Fisher's exact test p_value
fl2 <- fl %>%
  dplyr::group_by(.data$X..ChemicalName)%>%
  dplyr::mutate(p_value = fishertest(.data$GL_YES,.data$GL_NO,.data$UL_YES,.data$UL_NO))%>%
  dplyr::arrange(dplyr::desc(.data$p_value))%>%
  dplyr::ungroup()

#p value adjustment
fl2 <- fl2 %>%
  dplyr::mutate(bonf = stats::p.adjust(.data$p_value,"bonferroni")) %>%
  dplyr::mutate(FDR = stats::p.adjust(.data$p_value,"fdr")) %>%
  dplyr::arrange(.data$p_value)
fl3 <- fl2
fl3
}




###############

HPOlistfunction <- function(i){
  chem2gene_tbl <- tibble::as_tibble(chem2gene)
  pheno2gene_tbl <- tibble::as_tibble(pheno2gene)
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

  #tail(pheno_HPO_Cap)

  #Look up proper HPO term
  pheno_HPO_Cap <- pheno_HPO_Cap %>%
    dplyr::filter(.data$lowercase == HPO_low)

   #If HPO term is not in dataset return that it is not supported
  if ( length(pheno_HPO_Cap$HPO_Name) == 0) {
    stop(paste0(HPO_Name_Given," is not supported"))
  } else {
    #Assign HPO Name
    HPO <- as.character(pheno_HPO_Cap$HPO_Name)
    print(" ")
  }


  selectedpheno <-pheno %>%
    dplyr::filter(.data$HPO_Name == HPO)

  genelist <- unique(selectedpheno$GeneID)

}


#' Chemical enrichment for multiple phenotype terms
#'
#' @param HPO List of Human Phenotype Ontology terms
#' @param enrich_1S Enrichment one-sided. Automatically set to true for enrichment, where the one-sided Fisher's exact is set to alternative = 'greater.'
#' If enrich_1S = FALSE, Fisher's exact is set to alternative = 'two.sided' for two.sided.
#' @return Dataframe with enrichment.
#' @examples
#' \dontrun{
#' x <- list("HPO term1","HPO term2","HPO term3")
#' perfFishTestHPOMultiple(x)
#' }
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

  unlist_agenelist <- unlist(agenelist)

  #deduplicate
  genelist <-  unique(unlist_agenelist)
  #show(genelist)

  #Genes contained in both datasets
  chem2gene_tbl <- tibble::as_tibble(chem2gene)
  genesannotate <- Reduce(base::intersect, list(pheno2gene$GeneID,chem2gene$GeneID))

  #remove "Homo sapiens" -  #filter(Organism == "Homo sapiens") %>%

  #Filter
  # chems <- chem2gene_tbl %>%
  #   dplyr::filter(.data$GeneID %in% genesannotate )

   #Genelist
  gl <- chem2gene_tbl %>%
    dplyr::filter(.data$GeneID %in% genelist)%>%
    dplyr::group_by(.data$X..ChemicalName) %>%
    dplyr::select(.data$X..ChemicalName,.data$ChemicalID,.data$GeneID) %>%
    dplyr::distinct() %>%  #remove duplicate genes
    dplyr::summarise(GL_YES = dplyr::n() ) %>%
    dplyr::mutate(TOTAL_GL = length(genelist)) %>%
    dplyr::arrange(dplyr::desc(.data$GL_YES)) %>%
    dplyr::mutate(GL_NO = .data$TOTAL_GL - .data$GL_YES)

  #Universelist
  #filter(Organism == "Homo sapiens") %>% and all organisms
  #filter(GeneID %in% genesannotate ) %>%

  ul <- chem2gene_tbl %>%
    dplyr::group_by(.data$X..ChemicalName) %>%
    dplyr::select(.data$X..ChemicalName,.data$ChemicalID,.data$GeneID) %>%
    dplyr::distinct() %>%
    dplyr::summarise(UL_YES = dplyr::n() ) %>%
    dplyr::mutate(TOTAL_UL = length(unique(chem2gene_tbl$GeneID)))%>%
    dplyr::arrange(dplyr::desc(.data$UL_YES)) %>%
    dplyr::mutate(UL_NO = .data$TOTAL_UL - .data$UL_YES)

  #Join dataframes
  fl <- gl %>%
    dplyr::left_join(ul, by = "X..ChemicalName")


  #Fishers exact test
  fishertest <- function(a,b,c,d){
    m <- matrix(c(a,b,c,d),nrow=2)
    p <- stats::fisher.test(m,alternative = side)$p.value
  }

  #P-value
  fl2 <- fl %>%
    dplyr::group_by(.data$X..ChemicalName)%>%
    dplyr::mutate(p_value = fishertest(.data$GL_YES,.data$GL_NO,.data$UL_YES,.data$UL_NO))%>%
    dplyr::arrange(dplyr::desc(.data$p_value))%>%
    dplyr::ungroup()

  #P-value adjustment
  fl2 <- fl2 %>%
    dplyr::mutate(bonf = stats::p.adjust(.data$p_value,"bonferroni")) %>%
    dplyr::mutate(FDR = stats::p.adjust(.data$p_value,"fdr")) %>%
    dplyr::arrange(.data$p_value)
  fl3 <- fl2
  fl3
}


