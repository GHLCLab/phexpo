#' Chemical and gene interactions
#'
#' Dataset of chemicals and their gene interactions. Original data set downloaded on 11-Feb-2019. Version Feb-05-2019.
#' Dataset has been filtered to contain only the organisms Homo sapiens, Mus musculus, Rattus norvegicus,
#' Rattus rattus, Danio rerio, Macaca fascicularis, Macaca mulatta and human gene IDs. Genes were checked that
#' they belonged to humans using the Homo sapien gene info from NCBI (downloaded 19-March-2019).
#' For gene IDs not belonging to Homo sapiens, gene IDs were checked that they belonged to their respective organism
#' using All Mammalia gene info from NCBI (15-March-2019) and Danio rerio gene info from NCBI (15-March-2019).
#' All orthologs for these gene IDs were found using NCBI Homologene (version 06-May-2014) and outdated ortholog gene IDs
#' were updated to the latest gene ID using NCBI gene_history (15-March-2019) and updated gene symbols were taken from
#' Homo sapien gene info.
#'
#'  @format A data frame of 1,575,035 rows and 11 columns.
#'
#'
#'
#'  @source \url{http://ctdbase.org/downloads/#cg}
"chem2gene"


#' Phenotype and gene association
#'
#' Dataset of phenotype and their gene associations. Original data set downloaded on 28-Feb-2019. Version 2019-02-12.
#' Dataset has been filtered to only contain Phenotypic abnormality terms.
#'
#'  @format A data frame of 472,515 rows and 4 columns.
#'
#'
#'
#'  @source \url{https://hpo.jax.org/app/download/annotation}
"pheno2gene"

#' Chemical vocabularly
#'
#' Dataset of chemicals and their synonyms. Data has been modified from original data, to contain synonyms only for chemicals with gene interactions,
#' synonyms have been changed to lower-case and duplicate synonyms have been removed. Downloaded on 11-Feb-2019. Version Feb-05-2019.
#'
#'  @format A data frame of 63,818 rows and 2 columns.
#'
#'
#'
#'  @source \url{http://ctdbase.org/downloads/#allchems}
"sym_tbl"
