#' Chemical and gene interactions
#'
#' Dataset of chemicals and their gene interactions. Original dataset downloaded on 11-Feb-2019. Version Feb-05-2019.
#' Dataset has been filtered to contain only the organisms Homo sapiens, Mus musculus, Rattus norvegicus,
#' Rattus rattus, Danio rerio, Macaca fascicularis, Macaca mulatta and human gene IDs. Genes were checked that
#' they belonged to humans using the Homo sapien gene info from National Center for Biotechnology Information (NCBI) (downloaded 19-March-2019).
#' For gene IDs not belonging to Homo sapiens, gene IDs were checked that they belonged to their respective organism
#' using All Mammalia gene info from NCBI (15-March-2019) and Danio rerio gene info from NCBI (15-March-2019).
#' All orthologues for these gene IDs were found using NCBI Homologene (version 06-May-2014) and outdated orthologue gene IDs
#' were updated to the latest gene ID using NCBI gene_history (15-March-2019) and updated gene symbols were retrieved from
#' Homo sapien gene info. Tibble, dplyr, readr, tidyr, magrittr and annotationTools were used in the generation of this dataset.
#'
#'  @format A tibble data frame of 1,575,035 rows and 11 columns.
#'
#'
#'
#'  @source \url{http://ctdbase.org/downloads/#cg}
"chem2gene"


#' Phenotype and gene association
#'
#' Dataset of phenotype and their gene associations. Original dataset downloaded on 28-Feb-2019. Version 2019-02-12.
#' Dataset has been filtered to only contain phenotypic abnormality terms. Generation of this dataset used the Human Phenotype Ontology OBO file (2019-02-12 version)
#' and the R packages ontologyIndex, tibble, magrittr and dplyr.
#'
#'  @format A data frame of 472,515 rows and 4 columns.
#'
#'
#'
#'  @source \url{https://hpo.jax.org/app/download/annotation}
"pheno2gene"

#' Chemical vocabulary
#'
#' Dataset of chemicals and their synonyms. Data has been modified from original data to contain synonyms only for chemicals with gene interactions.
#' Duplicated synonyms for tetrachlorodibenzodioxin and polychlorinated dibenzodioxins were removed from polychlorinated dibenzodioxins,
#' synonyms have been changed to lowercase and duplicate synonyms have been removed. Downloaded on 11-Feb-2019. Version Feb-05-2019. The R packages dplyr, tidyr,
#' stringr, magrittr and tibble were used in the generation of this dataset.
#'
#'  @format A data frame of 63,819 rows and 2 columns.
#'
#'
#'
#'  @source \url{http://ctdbase.org/downloads/#allchems}
"sym_tbl"
