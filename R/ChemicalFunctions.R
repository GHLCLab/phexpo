#Global variables - chem2gene, pheno2gene and sym_tbl - lazy data
#https://github.com/tidyverse/magrittr/issues/29
if(getRversion() >= "2.15.1")  utils::globalVariables(c("chem2gene", "pheno2gene","sym_tbl","."))




#' Enriched HPO terms for a single chemical's associated gene list
#'
#' @param cname Chemical character
#' @param enrich_1S Enrichment one-sided. Automatically set to true for enrichment, where the one-sided Fisher's exact is set to alternative = 'greater.'
#' If enrich_1S = FALSE, Fisher's exact is set to alternative = "two.sided" for two-sided.
#' @return Dataframe with enrichment.
#' @note Chemical names and synonyms used are dependent on the chemical names within CTD.
#' @examples
#' perfFishTestChemSingle("Iodine")
#' @import magrittr shiny
#' @importFrom rlang .data
#' @export

perfFishTestChemSingle <- function(cname, enrich_1S = TRUE){

  #Change chemical name to lower case
  cname <- tolower(cname)

  #Check chemical is supported
  if (cname %in% unique(sym_tbl$synonym) == FALSE) {
    stop(paste0(cname," is not supported"))
  } else {
    print(" ")
  }

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
  ##filter(Organism == "Homo sapiens") %>% in chems removed
  #Filter
  chems <- chem2gene_tbl %>%
      dplyr::filter(.data$GeneID %in% genesannotate )

  selectedchem <-chems %>%
       dplyr::filter(.data$ChemicalID == sel_chem2)
  #head(selectedchem$GeneID)
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


  #Universe list
  ul <- pheno2gene_tbl %>%
    dplyr::group_by(.data$HPO_Name) %>%
    dplyr::summarise(UL_YES = dplyr::n() ) %>%
    dplyr::mutate(TOTAL_UL = length(unique(pheno2gene$GeneID)))%>%
    dplyr::arrange(dplyr::desc(.data$UL_YES)) %>%
    dplyr::mutate(UL_NO = .data$TOTAL_UL - .data$UL_YES)

  #Join the dataframes - left_join
  fl <- gl %>%
    dplyr::left_join(ul, by = "HPO_Name")


  #Fisher's exact test
  fishertest <- function(a,b,c,d){
    m <- matrix(c(a,b,c,d),nrow=2)
    p <- stats::fisher.test(m,alternative = side)$p.value
  }

  #P value
    fl2 <- fl %>%
    dplyr::group_by(.data$HPO_Name)%>%
    dplyr::mutate(p_value = fishertest(.data$GL_YES,.data$GL_NO,.data$UL_YES,.data$UL_NO))%>%
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

  # print(terms_addition)
  terms_addition <- unique(terms_addition) #just need the term once
  #terms_addition_2 <- as_tibble(terms_addition)
  terms_addition_2 <- tibble::enframe(terms_addition)
  terms_addition_2 <- terms_addition_2[,-1]
  colnames(terms_addition_2) <- "HPO_Name"
  #print(terms_addition2)

  terms_addition_3 <- terms_addition_2 %>%
       dplyr::mutate(GL_YES = 0,
                  TOTAL_GL = length(genelist),
                  GL_NO = length(genelist),
                  p_value = NA,
                  bonf = NA,
                  FDR = NA)

  terms_addition_ul <- terms_addition_3 %>%
  dplyr::left_join(ul, by = "HPO_Name")

  #Bind rows
  fl4 <- dplyr::bind_rows(fl3,terms_addition_ul)

fl4

}


#### Multiple chemicals

genelistfunction <- function(i){
  chem2gene_tbl <- tibble::as_tibble(chem2gene)
  pheno2gene_tbl <- tibble::as_tibble(pheno2gene)
  genesannotate <- Reduce(base::intersect, list(pheno2gene$GeneID,chem2gene$GeneID))
  #filter(Organism == "Homo sapiens") %>% removed in chems
  chems <- chem2gene_tbl %>%
    dplyr::filter(.data$GeneID %in% genesannotate )

  #Change chemical name to lower case
  i <- tolower(i)

  #Check chemical is supported
  if (i %in% unique(sym_tbl$synonym) == FALSE) {
    stop(paste0(i," is not supported"))
  } else {
    print(" ")
  }
  #
  sym_tbl <- tibble::as_tibble(sym_tbl)
  #Synonym
  sel_chem <- sym_tbl %>%
    dplyr::filter(.data$synonym == i)

  sel_chem2 <- sel_chem$ChemicalID
  #show(sel_chem2)

  selectedchem  <- chems %>%
    dplyr::filter(.data$ChemicalID == sel_chem2)

  #show(selectedchem)
  #head(selectedchem$GeneID)

  #Genelist
  genelist <-  unique(selectedchem$GeneID)
  #print(genelist)

}


#' Enriched HPO terms for multiple chemicals
#'
#' @param cname List of chemical characters
#' @param enrich_1S Enrichment one-sided. Automatically set to true for enrichment, where the one-sided Fisher's exact is set to alternative = 'greater.'
#' If enrich = FALSE, Fisher's exact is set to alternative = 'two.sided' for two-sided.
#' @return Dataframe with enrichment
#' @examples
#' \dontrun{
#' x <- list("Chemical1","Chemical2","Chemical3")
#' perfFishTestChemMultiple(x)
#' }
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

  agenelist <- lapply(cname, genelistfunction)
  #print("genelists")
  #show(agenelist)
  #print("combinedlists")
  unlist_agenelist <- unlist(agenelist)
  #show(unlist_agenelist)
  #print("removed duplicates")
  genelist <-  unique(unlist_agenelist)
  #show(genelist)

  pheno2gene_tbl <- tibble::as_tibble(pheno2gene)

  #Genelist
  gl <- pheno2gene_tbl %>%
    dplyr::filter(.data$GeneID %in% genelist)%>%
    dplyr::group_by(.data$HPO_Name) %>%
    dplyr::summarise(GL_YES = dplyr::n() ) %>%
    dplyr::mutate(TOTAL_GL = length(genelist)) %>%
    dplyr::arrange(dplyr::desc(.data$GL_YES)) %>%
    dplyr::mutate(GL_NO = .data$TOTAL_GL - .data$GL_YES)

#print(gl)

  ul <- pheno2gene_tbl %>%
    dplyr::group_by(.data$HPO_Name) %>%
    dplyr::summarise(UL_YES = dplyr::n() ) %>%
    dplyr::mutate(TOTAL_UL = length(unique(pheno2gene$GeneID)))%>%
    dplyr::arrange(dplyr::desc(.data$UL_YES)) %>%
    dplyr::mutate(UL_NO = .data$TOTAL_UL - .data$UL_YES)

  # #join them
  fl <- gl %>%
    dplyr::left_join(ul,by = "HPO_Name")


  #Fishertest function
  fishertest <- function(a,b,c,d){
    m <- matrix(c(a,b,c,d),nrow=2)
    p <- stats::fisher.test(m,alternative = side)$p.value
  }

  fl2 <- fl %>%
    dplyr::group_by(.data$HPO_Name)%>%
    dplyr::mutate(p_value = fishertest(.data$GL_YES,.data$GL_NO,.data$UL_YES,.data$UL_NO))%>%
    dplyr::arrange(dplyr::desc(.data$p_value)) %>%
    dplyr::ungroup()

  fl2 <- fl2 %>%
    dplyr::mutate(bonf = stats::p.adjust(.data$p_value,"bonferroni")) %>%
    dplyr::mutate(FDR = stats::p.adjust(.data$p_value,"fdr")) %>%
    dplyr::arrange(.data$p_value)
  fl3 <- fl2

  #the HPO terms that have been filtered out, need to be re-added with a GL of 0
  #terms that are only in pheno2gene and not in fl3, hence the terms that are filtered out
  terms_addition <- pheno2gene_tbl$HPO_Name[!pheno2gene_tbl$HPO_Name %in% fl3$HPO_Name]

  # print(terms_addition)
  terms_addition <- unique(terms_addition) #just need the term once
  #terms_addition_2 <- as_tibble(terms_addition)
  terms_addition_2 <- tibble::enframe(terms_addition)
  terms_addition_2 <- terms_addition_2[,-1]
  colnames(terms_addition_2) <- "HPO_Name"
  #print(terms_addition2)

  terms_addition_3 <- terms_addition_2 %>%
    dplyr::mutate(GL_YES = 0,
                  TOTAL_GL = length(genelist),
                  GL_NO = length(genelist),
                  p_value = NA,
                  bonf = NA,
                  FDR = NA)

  terms_addition_ul <- terms_addition_3 %>%
    dplyr::left_join(ul, by = "HPO_Name")

  #Bind rows
  fl4 <- dplyr::bind_rows(fl3,terms_addition_ul)

  fl4


}


