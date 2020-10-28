#' Wrapper function for UCSC's liftOver tool
#'
#' Converts coordinates in genomic range object from the given assembly to the target assembly
#'
#' @param grObject GRanges object with genomic coordinates (see GenomicRanges documentation for more info)
#' @param givenAssembly genome assembly of the supplied coordinates (e.g. "hg19")
#' @param targetAssembly genome assembly to which the coordinates will be converted
#' @param libPath directory where the liftOver chain is stored (the chain will be downloaded if not already present in this directory)
#'
#' @return GRanges object with genomic coordinates after liftOver conversion
#'
#' @export

liftOverGRange <- function(grObject, givenAssembly, targetAssembly, libPath){

  targetAssembly <- (paste0(toupper(substring(targetAssembly,1,1)), substring(targetAssembly,2)))

  chain <- paste0(givenAssembly, "To", targetAssembly, ".over.chain")
  chain_path <- paste0(libPath, "/", chain)

  chain_url <- paste0("https://hgdownload.soe.ucsc.edu/goldenPath/",
                      givenAssembly,"/liftOver/", chain, ".gz")

  if (!chain %in% list.files(libPath)){

    tryCatch({
      download.file(chain_url, destfile = paste0(chain_path, ".gz"), method = "wget", quiet = TRUE)
      R.utils::gunzip(paste0(chain_path, ".gz"), remove = TRUE)
    },
    error = function(cond){
      message("Error in downloading liftOver chain:")
      message(paste0("  given genome assembly (= \"", givenAssembly, "\") cannot be found"))
    }
    )
  }

  #seqlevelsStyle(grObject) = "UCSC"  # necessary to set "chr" format

  try({
    chain <- rtracklayer::import.chain(chain_path)
    grObject <- unlist(rtracklayer::liftOver(grObject, chain))
    #genome(grObject) <- tolower(targetAssembly)
    grObject$Assembly <- tolower(targetAssembly)
    return(grObject)
  }, silent = TRUE )

  return(NULL)
}
