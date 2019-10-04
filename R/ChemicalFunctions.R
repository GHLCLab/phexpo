#Global variables - chem2gene, pheno2gene and sym_tbl - lazy data
#https://github.com/tidyverse/magrittr/issues/29
if(getRversion() >= "2.15.1")  utils::globalVariables(c("chem2gene", "pheno2gene","sym_tbl","."))




#' Enriched HPO terms for a single chemical's associated gene list
#'
#' @param cname Chemical character.
#' @param enrich_1S Enrichment one-sided. Default set to true for enrichment, where the one-sided Fisher's exact is set to alternative = "greater."
#' If enrich_1S = FALSE, Fisher's exact is set to alternative = "two.sided" for two-sided.
#' @return Tibble data frame with enrichment.
#' @note Chemical names and synonyms used are dependent on the chemical names within CTD.
#' @examples
#' perfFishTestChemSingle("Iodine")
#' @import magrittr shiny
#' @importFrom rlang .data
#' @export

perfFishTestChemSingle <- function(cname, enrich_1S = TRUE){

  #Change chemical name to lower case
  cname <- tolower(cname)

  #Enrichment - True or False
  if (enrich_1S == TRUE) {
    side <- "greater"
  } else if (enrich_1S == FALSE) {
    side <- "two.sided"
  } else
    stop("enrich_1S should be TRUE or FALSE")

  #as tibble
  sym_tbl <- tibble::as_tibble(sym_tbl)
  #Look up synonyms to chemical ID
  sel_chem <- sym_tbl %>%
    dplyr::filter(.data$synonym == cname)
  sel_chem2 <- sel_chem$ChemicalID
  chem2gene_tbl <- tibble::as_tibble(chem2gene)
  pheno2gene_tbl <- tibble::as_tibble(pheno2gene)
  #Genes in both datasets
  genesannotate <- Reduce(base::intersect, list(pheno2gene$GeneID,chem2gene$GeneID))

  #Filter
  chems <- chem2gene_tbl %>%
      dplyr::filter(.data$GeneID %in% genesannotate )


  #Check chemical is supported
  if (cname %in% unique(sym_tbl$synonym) == FALSE) {
    stop(paste0(cname," is not supported"))
  } else if (sel_chem2 %in% chems$ChemicalID == FALSE){
    stop(paste0(cname," is not supported"))
  } else {
    print("Chemical accepted")
  }

  #Filter
  selectedchem <- chems %>%
       dplyr::filter(.data$ChemicalID == sel_chem2)

 #Unique genelist
  genelist <-  unique(selectedchem$GeneID)

  #Genelist
  gl <- pheno2gene_tbl %>%
    dplyr::filter(.data$GeneID %in% genelist)%>%
    dplyr::group_by(.data$HPO_Name) %>%
    dplyr::summarise(GL_YES = dplyr::n() ) %>%
    dplyr::mutate(TOTAL_GL = length(genelist)) %>%
    dplyr::arrange(dplyr::desc(.data$GL_YES)) %>%
    dplyr::mutate(GL_NO = .data$TOTAL_GL - .data$GL_YES)

  #NOT Genelist
  Not_Genelist <- pheno2gene_tbl %>%
    dplyr::filter(!.data$GeneID %in% genelist)

  #Not genelist
  nl <- Not_Genelist %>%
    dplyr::group_by(.data$HPO_Name) %>%
    dplyr::summarise(NOT_GL_YES = dplyr::n() ) %>%
    dplyr::mutate(TOTAL_NOT_GL = length(unique(Not_Genelist$GeneID)))%>%
    dplyr::arrange(dplyr::desc(.data$NOT_GL_YES)) %>%
    dplyr::mutate(NOT_GL_NO = .data$TOTAL_NOT_GL - .data$NOT_GL_YES)

  #The HPO terms that have been filtered out because they do not have any genes in the not gene list
  #need to be re-added with a NOT_GL of 0 in nl
  HPOterms_addition <- pheno2gene_tbl$HPO_Name[!pheno2gene_tbl$HPO_Name %in% nl$HPO_Name]
  HPOterms_addition <- unique(HPOterms_addition) #just need the term once
  HPOterms_addition_2 <- tibble::enframe(HPOterms_addition)
  HPOterms_addition_2 <- HPOterms_addition_2[,-1]
  colnames(HPOterms_addition_2) <- "HPO_Name"


  HPOterms_addition_3 <- HPOterms_addition_2 %>%
    dplyr::mutate(NOT_GL_YES = 0,
                  TOTAL_NOT_GL = length(unique(Not_Genelist$GeneID)),
                  NOT_GL_NO = .data$TOTAL_NOT_GL - .data$NOT_GL_YES)

  #Bind the nl(Not_Genelist) and the HPOterms_addition_3
  nl_2 <- dplyr::bind_rows(nl,HPOterms_addition_3)

  #Join the dataframes - left_join
  fl <- gl %>%
    dplyr::left_join(nl_2, by = "HPO_Name")


  #Fisher's exact test
  #Attribution: https://stackoverflow.com/questions/28966840/fishers-test-stat-on-data-frame-using-dplyrmutate
  fishertest <- function(a,b,c,d){
    m <- matrix(c(a,b,c,d),nrow=2)
    p <- stats::fisher.test(m,alternative = side)$p.value
  }

  #P value
  fl2 <- fl %>%
    dplyr::group_by(.data$HPO_Name)%>%
    dplyr::mutate(p_value = fishertest(.data$GL_YES,.data$NOT_GL_YES,.data$GL_NO,.data$NOT_GL_NO))%>%
    dplyr::arrange(dplyr::desc(.data$p_value)) %>%
    dplyr::ungroup()

  #P value adjustment
  fl2 <- fl2 %>%
    dplyr::mutate(bonf = stats::p.adjust(.data$p_value,"bonferroni")) %>%
    dplyr::mutate(FDR = stats::p.adjust(.data$p_value,"fdr")) %>%
    dplyr::arrange(.data$p_value)
  fl3 <- fl2

  #The HPO terms that have been filtered out, need to be re-added with a GL of 0
  #terms that are only in pheno2gene and not in fl3, hence the terms that are filtered out
  terms_addition <- pheno2gene_tbl$HPO_Name[!pheno2gene_tbl$HPO_Name %in% fl3$HPO_Name]
  terms_addition <- unique(terms_addition) #just need the term once
  terms_addition_2 <- tibble::enframe(terms_addition)
  terms_addition_2 <- terms_addition_2[,-1]
  colnames(terms_addition_2) <- "HPO_Name"

  terms_addition_3 <- terms_addition_2 %>%
    dplyr::mutate(GL_YES = 0,
                  TOTAL_GL = length(genelist),
                  GL_NO = length(genelist),
                  p_value = NA,
                  bonf = NA,
                  FDR = NA)

  terms_addition_nl <- terms_addition_3 %>%
  dplyr::left_join(nl_2, by = "HPO_Name")

  #Bind rows
  fl4 <- dplyr::bind_rows(fl3,terms_addition_nl)

  fl4

}


#### Multiple chemicals

genelistfunction <- function(i){
  chem2gene_tbl <- tibble::as_tibble(chem2gene)
  pheno2gene_tbl <- tibble::as_tibble(pheno2gene)
  #Genes in both datasets
  genesannotate <- Reduce(base::intersect, list(pheno2gene$GeneID,chem2gene$GeneID))

  #Filter
  chems <- chem2gene_tbl %>%
    dplyr::filter(.data$GeneID %in% genesannotate )

  #Change chemical name to lower case
  i <- tolower(i)

  #tibble
  sym_tbl <- tibble::as_tibble(sym_tbl)
  #Synonym
  sel_chem <- sym_tbl %>%
    dplyr::filter(.data$synonym == i)

  sel_chem2 <- sel_chem$ChemicalID


  #Check chemical is supported
  if (i %in% unique(sym_tbl$synonym) == FALSE) {
    stop(paste0(i," is not supported"))
  } else if (sel_chem2 %in% chems$ChemicalID == FALSE){
    stop(paste0(i," is not supported"))
  } else {
    print("Chemical accepted")
  }

  #Filter
  selectedchem  <- chems %>%
    dplyr::filter(.data$ChemicalID == sel_chem2)

  #Genelist
  genelist <-  unique(selectedchem$GeneID)


}


#' Enriched HPO terms for multiple chemicals
#'
#' @param cname List of chemical characters.
#' @param enrich_1S Enrichment one-sided. Default set to true for enrichment, where the one-sided Fisher's exact is set to alternative = "greater."
#' If enrich_1S = FALSE, Fisher's exact is set to alternative = "two.sided" for two-sided.
#' @return Tibble data frame with enrichment.
#' @examples
#' perfFishTestChemMultiple(list("bisphenol A","bisphenol F","bisphenol B"))
#' @note Enrichment for multiple chemicals is a simplistic exploratory function, where we take an additive approach.
#' Chemical names and synonyms used are dependent on the chemical names within CTD.
#' @export
perfFishTestChemMultiple<- function(cname, enrich_1S = TRUE){

  #Enrichment - True or False
  if (enrich_1S == TRUE) {
    side <- "greater"
  } else if (enrich_1S == FALSE) {
    side <- "two.sided"
  } else
    stop("enrich_1S should be TRUE or FALSE")

  #Create genelist
  agenelist <- lapply(cname, genelistfunction)
  unlist_agenelist <- unlist(agenelist)
  genelist <-  unique(unlist_agenelist)


  pheno2gene_tbl <- tibble::as_tibble(pheno2gene)

  #Genelist
  gl <- pheno2gene_tbl %>%
    dplyr::filter(.data$GeneID %in% genelist)%>%
    dplyr::group_by(.data$HPO_Name) %>%
    dplyr::summarise(GL_YES = dplyr::n() ) %>%
    dplyr::mutate(TOTAL_GL = length(genelist)) %>%
    dplyr::arrange(dplyr::desc(.data$GL_YES)) %>%
    dplyr::mutate(GL_NO = .data$TOTAL_GL - .data$GL_YES)

  #NOT Genelist
  Not_Genelist <- pheno2gene_tbl %>%
    dplyr::filter(!.data$GeneID %in% genelist)

  #NOT Genelist
  nl <- Not_Genelist %>%
    dplyr::group_by(.data$HPO_Name) %>%
    dplyr::summarise(NOT_GL_YES = dplyr::n() ) %>%
    dplyr::mutate(TOTAL_NOT_GL = length(unique(Not_Genelist$GeneID)))%>%
    dplyr::arrange(dplyr::desc(.data$NOT_GL_YES)) %>%
    dplyr::mutate(NOT_GL_NO = .data$TOTAL_NOT_GL - .data$NOT_GL_YES)

  #The HPO terms that have been filtered out because they do not have any genes in the not gene list
  #need to be re-added with a NOT_GL of 0 in nl
  HPOterms_addition <- pheno2gene_tbl$HPO_Name[!pheno2gene_tbl$HPO_Name %in% nl$HPO_Name]
  HPOterms_addition <- unique(HPOterms_addition) #just need the term once
  HPOterms_addition_2 <- tibble::enframe(HPOterms_addition)
  HPOterms_addition_2 <- HPOterms_addition_2[,-1]
  colnames(HPOterms_addition_2) <- "HPO_Name"

  HPOterms_addition_3 <- HPOterms_addition_2 %>%
    dplyr::mutate(NOT_GL_YES = 0,
                  TOTAL_NOT_GL = length(unique(Not_Genelist$GeneID)),
                  NOT_GL_NO = .data$TOTAL_NOT_GL - .data$NOT_GL_YES)

  #Bind the nl(Not_Genelist) and the HPOterms_addition_3
  nl_2 <- dplyr::bind_rows(nl,HPOterms_addition_3)

  #Join the dataframes - left_join
  fl <- gl %>%
    dplyr::left_join(nl_2,by = "HPO_Name")


  #Fisher's exact test
  #Attribution: https://stackoverflow.com/questions/28966840/fishers-test-stat-on-data-frame-using-dplyrmutate
  fishertest <- function(a,b,c,d){
    m <- matrix(c(a,b,c,d),nrow=2)
    p <- stats::fisher.test(m,alternative = side)$p.value
  }

  #P value
  fl2 <- fl %>%
    dplyr::group_by(.data$HPO_Name)%>%
    dplyr::mutate(p_value = fishertest(.data$GL_YES,.data$NOT_GL_YES,.data$GL_NO,.data$NOT_GL_NO))%>%
    dplyr::arrange(dplyr::desc(.data$p_value)) %>%
    dplyr::ungroup()

  #P value adjustment
  fl2 <- fl2 %>%
    dplyr::mutate(bonf = stats::p.adjust(.data$p_value,"bonferroni")) %>%
    dplyr::mutate(FDR = stats::p.adjust(.data$p_value,"fdr")) %>%
    dplyr::arrange(.data$p_value)
  fl3 <- fl2

  #The HPO terms that have been filtered out, need to be re-added with a GL of 0
  #terms that are only in pheno2gene and not in fl3, hence the terms that are filtered out
  terms_addition <- pheno2gene_tbl$HPO_Name[!pheno2gene_tbl$HPO_Name %in% fl3$HPO_Name]
  terms_addition <- unique(terms_addition) #just need the term once
  terms_addition_2 <- tibble::enframe(terms_addition)
  terms_addition_2 <- terms_addition_2[,-1]
  colnames(terms_addition_2) <- "HPO_Name"

  terms_addition_3 <- terms_addition_2 %>%
    dplyr::mutate(GL_YES = 0,
                  TOTAL_GL = length(genelist),
                  GL_NO = length(genelist),
                  p_value = NA,
                  bonf = NA,
                  FDR = NA)

  terms_addition_nl <- terms_addition_3 %>%
    dplyr::left_join(nl_2, by = "HPO_Name")

  #Bind rows
  fl4 <- dplyr::bind_rows(fl3,terms_addition_nl)

  fl4


}


