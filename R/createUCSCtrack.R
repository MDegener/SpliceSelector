#' UCSC custom track generator
#'
#' Creates a custom UCSC track from a genomic range object to view supplied coordinates in the UCSC genome browser
#'
#' @param grObject GRanges object with genomic coordinates (see GenomicRanges documentation for more info)
#' @param outFile filepath and name of output UCSC track
#' @param trackName name of UCSC track
#' @param trackDescription description of UCSC track
#'
#' @return None
#'
#' @export

createUCSCtrack <- function(grObject, outFile, trackName, trackDescription){

  # create output file
  file.create(outFile)

  # open connection to file
  fileConn<- file(outFile)

  # export UCSC track as .bed file and add trackLine
  rtracklayer::export.ucsc(grObject, con = fileConn, subformat = "BED",
                           name = trackName,
                           description = trackDescription,
                           visibility = "2",
                           priority = 1)
}

