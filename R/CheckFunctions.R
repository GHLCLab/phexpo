#supportChemical function
supportChemical <- function(i){
  chem2gene_tbl <- tibble::as_tibble(chem2gene)
  #Genes in both datasets
  genesannotate <- Reduce(base::intersect, list(pheno2gene$GeneID,chem2gene$GeneID))

  #Filter
  chems <- chem2gene_tbl %>%
    dplyr::filter(.data$GeneID %in% genesannotate)

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
    support <- "No"
  } else if (sel_chem2 %in% chems$ChemicalID == FALSE){
    support <- "No"
  } else {
    support <- "Yes"
  }
  support
}

#' Check chemicals are supported for analysis
#'
#' @param cname List of chemical characters.
#'
#' @return Tibble data frame with chemicals and if they are supported or not.
#' @examples checkChemical(list("Ethanol","Iodine","Zinc"))
#' @export
checkChemical <- function(cname){
  #Create list of chemicals
  chemicallist <- lapply(cname, supportChemical)
  #Add chemical names
  names(chemicallist) <- cname
  #Make tibble and unnest the list
  chemicallist <- tibble::enframe(chemicallist) %>% tidyr::unnest(.data$value)
  colnames(chemicallist) <- c("Chemical_Name", "Supported")
  chemicallist
}


#supportHPO function
supportHPO <- function(i){
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
    support <- "No"
  } else {
    support <- "Yes"
  }
  support
}

#' Check HPO terms are supported for analysis
#'
#' @param HPO List of Human Phenotype Ontology terms characters.
#'
#' @return Tibble data frame with HPO terms and if they are supported or not.
#' @examples checkHPO(list("Rickets","Breast carcinoma","Preeclampsia"))
#' @export
#'
checkHPO <- function(HPO){
  #Create list of HPO terms
  HPOlist <- lapply(HPO, supportHPO)
  #Add HPO names
  names(HPOlist) <- HPO
  #Make tibble and unnest the list
  HPOlist <- tibble::enframe(HPOlist) %>% tidyr::unnest(.data$value)
  colnames(HPOlist) <- c("HPO_Name", "Supported")
  HPOlist
}


