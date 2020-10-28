#' Get gene annotation from biomaRt Ensembl database
#'
#' Get gene annotation based on identifier (e.g. Ensembl gene id) using biomaRt Ensembl database
#'
#' @param identifierList values of filters that are applied to biomaRt query
#' @param identifierType types of filters that are applied to biomaRt query
#' @param species should be specified as "hsapiens" or "mmusculus"
#'
#' @return Annotation retrieved from biomaRt database
#'
#' @export

getBiomaRtAnnotation <- function(identifierList, identifierType, species){

  ensembl <- biomaRt::useMart(biomart="ensembl",dataset=paste0(species,"_gene_ensembl"))

  geneSymbol <- if (species == "hsapiens") "hgnc_symbol" else if (species == "mmusculus") "mgi_symbol"

  att <- c("ensembl_gene_id", geneSymbol, "chromosome_name",
           "start_position", "end_position", "strand", "description")

  ann <- biomaRt::getBM(attributes=att, filters=identifierType, values=identifierList, mart=ensembl)

  return(ann)
}
