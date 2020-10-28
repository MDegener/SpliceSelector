#' Overlap between genomic and exon coordinates
#'
#' Compares supplied coordinates with coordinates from an exon annotation file (e.g. GENCODE .gtf)
#'
#' @param grObject GRanges object with genomic coordinates (see GenomicRanges documentation for more info)
#' @param exonAnnotation annotation file filtered for exon coordinates (e.g. from GENCODE)
#'
#' @return Table listing all overlapping coordinates with their corresponding overlap type
#'
#' @export

getOverlap <- function(grObject, exonAnnotation){

  # find coordinates that overlap with an annotated exon in any way
  anyOverlap <- GenomicRanges::findOverlaps(grObject, exonAnnotation, type = "any")

  # find coordinates that are equal to start/end position of an annotated exon (= exact matches)
  exactOverlap <- GenomicRanges::findOverlaps(grObject, exonAnnotation, type = "equal")

  # find coordinates that fall within an annotated exon
  withinOverlap <- GenomicRanges::findOverlaps(grObject, exonAnnotation, type = "within")

  # find coordinates that span an annotated exon (equals to exons that fall within the given coordinates)
  # note: result is transposed to have same structure as the other results
  spanOverlap <- t(GenomicRanges::findOverlaps(exonAnnotation, grObject, type = "within"))

  # find coordinates that are not overlapping with the exon annotation
  noOverlap <- grObject[ grObject %outside% exonAnnotation ]

  # create empty dataframe
  overlap <- data.frame()

  # check if there is any overlap
  if( length(anyOverlap) != 0) {

    if( length(exactOverlap) != 0 ){

      # remove redundant exact matches from withinOverlap, spanOverlap and anyOverlap
      withinOverlap <- withinOverlap[ -which(withinOverlap %in% exactOverlap) ]
      spanOverlap <- spanOverlap[ -which(spanOverlap %in% exactOverlap) ]
      anyOverlap <- anyOverlap[ -which(anyOverlap %in% exactOverlap) ]

      # save matches in dataframe and annotate with its specific type
      overlap <- dplyr::bind_rows(overlap, data.frame(overlap = "equal",
                                                      query.index = from(exactOverlap),
                                                      query = grObject[ from(exactOverlap) ],
                                                      annotation.index = to(exactOverlap),
                                                      annotation = exonAnnotation[ to(exactOverlap) ],
                                                      stringsAsFactors = FALSE))
    }

    if( length(withinOverlap) != 0 ){

      # remove redundant within matches from anyOverlap
      anyOverlap <- anyOverlap[ -which(anyOverlap %in% withinOverlap) ]

      overlap <- dplyr::bind_rows(overlap, data.frame(overlap = "within",
                                                      query.index = from(withinOverlap),
                                                      query = grObject[ from(withinOverlap) ],
                                                      annotation.index = to(withinOverlap),
                                                      annotation = exonAnnotation[ to(withinOverlap) ],
                                                      stringsAsFactors = FALSE))
    }

    if( length(spanOverlap) != 0 ){

      # remove redundant spanning matches from anyOverlap
      anyOverlap <- anyOverlap[ -which(anyOverlap %in% spanOverlap) ]

      # first store spanOverlap in dataframe to allow easy access to annotation
      spanOverlap <- data.table::data.table(overlap = "",
                                            query.index = from(spanOverlap),
                                            query = grObject[ from(spanOverlap) ],
                                            annotation.index = to(spanOverlap),
                                            annotation = exonAnnotation[ to(spanOverlap) ],
                                            stringsAsFactors = FALSE)

      # distinguish coordinates that span one or multiple exons
      spanMatches <- spanOverlap %>%
        dplyr::group_by(query.index) %>%
        dplyr::count(annotation.exon_id) %>%
        dplyr::summarise(sum(n))
      spanOverlap$overlap <- ifelse(spanMatches[ match(spanOverlap$query.index, spanMatches$query.index), 2] == 1,
                                    "spans one", "spans multiple")

      overlap <- dplyr::bind_rows(overlap, spanOverlap)
    }

    # save matches that are neither exact, within or spanning
    if( length(anyOverlap) != 0 ){
      overlap <- dplyr::bind_rows(overlap, data.frame(overlap = "other",
                                                      query.index = from(anyOverlap),
                                                      query = grObject[ from(anyOverlap) ],
                                                      annotation.index = to(anyOverlap),
                                                      annotation = exonAnnotation[ to(anyOverlap) ],
                                                      stringsAsFactors = FALSE))
    }
  }

  if( length(noOverlap) != 0 ){
    overlap <- dplyr::bind_rows(overlap, data.table::data.table(overlap = "none",
                                                                query.index = which(grObject %in% noOverlap),
                                                                query = noOverlap,
                                                                stringsAsFactors = FALSE))
  }

  # remove duplicate rows
  overlap <- unique(overlap)

  return(overlap)
}
