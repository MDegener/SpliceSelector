# TODO: split barplot (first two bars / last four bars)

getOverlap <- function(grObject, exonAnnotation){
  
  library(rtracklayer)
  library(dplyr)
  
  # find coordinates that overlap with an annotated exon in any way
  anyOverlap <- findOverlaps(grObject, exonAnnotation, type = "any")
  
  # find coordinates that are equal to start/end position of an annotated exon (= exact matches)
  exactOverlap <- findOverlaps(grObject, exonAnnotation, type = "equal")
  
  # find coordinates that fall within an annotated exon
  withinOverlap <- findOverlaps(grObject, exonAnnotation, type = "within")
  
  # find coordinates that span an annotated exon (equals to exons that fall within the given coordinates)
  # note: result is transposed to have same structure as the other results
  spanOverlap <- t(findOverlaps(exonAnnotation, grObject, type = "within"))

  # find coordinates that are not overlapping with the exon annotation
  noOverlap <- grObject[ grObject %outside% exonAnnotation ]
  
  # create empty dataframe
  overlap <- data.frame()
  
  # check if there is any overlap
  if( length(anyOverlap) != 0) {
    
    if( length(exactOverlap) != 0 ){
      
      # save matches in dataframe and annotate with its specific type
      overlap <- bind_rows(overlap, data.frame(overlap = "equal",
                                               query.index = from(exactOverlap),
                                               query = grObject[ from(exactOverlap) ],
                                               annotation.index = to(exactOverlap),
                                               annotation = exonAnnotation[ to(exactOverlap) ],
                                               stringsAsFactors = FALSE))
      
      # remove redundant exact matches from withinOverlap, spanOverlap and anyOverlap
      withinOverlap <- withinOverlap[ -which(withinOverlap %in% exactOverlap) ]
      spanOverlap <- spanOverlap[ -which(spanOverlap %in% exactOverlap) ] 
      anyOverlap <- anyOverlap[ -which(anyOverlap %in% exactOverlap) ]
    }
    
    if( length(withinOverlap) != 0 ){
      overlap <- bind_rows(overlap, data.frame(overlap = "within",
                                               query.index = from(withinOverlap),
                                               query = grObject[ from(withinOverlap) ], 
                                               annotation.index = to(withinOverlap),
                                               annotation = exonAnnotation[ to(withinOverlap) ],
                                               stringsAsFactors = FALSE))
      
      # remove redundant within matches from anyOverlap
      anyOverlap <- anyOverlap[ -which(anyOverlap %in% withinOverlap) ]
    }
    
    if( length(spanOverlap) != 0 ){
      overlap <- bind_rows(overlap, data.frame(overlap = "spans",
                                               query.index = from(spanOverlap),
                                               query = grObject[ from(spanOverlap) ],
                                               annotation.index = to(spanOverlap),
                                               annotation = exonAnnotation[ to(spanOverlap) ],
                                               stringsAsFactors = FALSE))
      
      # remove redundant spanning matches from anyOverlap
      anyOverlap <- anyOverlap[ -which(anyOverlap %in% spanOverlap) ]
    }
    
    # save matches that are neither exact, within or spanning
    if( length(anyOverlap) != 0 ){
      overlap <- bind_rows(overlap, data.frame(overlap = "other",
                                               query.index = from(anyOverlap),
                                               query = grObject[ from(anyOverlap) ],
                                               annotation.index = to(anyOverlap),
                                               annotation = exonAnnotation[ to(anyOverlap) ],
                                               stringsAsFactors = FALSE))
    }
  }
  
  if( length(noOverlap) != 0 ){
    overlap <- bind_rows(overlap, data.frame(overlap = "none",
                                             query.index = which(grObject %in% noOverlap),
                                             query = noOverlap,
                                             stringsAsFactors = FALSE))
  }
 
  # remove duplicate rows
  overlap <- unique(overlap)
  
  # get counts for input coordinates, all detected matches and all overlap types
  counts <- c( grObject %>% length(),
               overlap[which(overlap$overlap != "none"), ] %>% nrow(),
               overlap[which(overlap$overlap == "equal"), ] %>% nrow(),
               overlap[which(overlap$overlap == "within"), ] %>% nrow(),
               overlap[which(overlap$overlap == "spans"), ] %>% nrow(),
               overlap[which(overlap$overlap == "other"), ] %>% nrow(),
               overlap[which(overlap$overlap == "none"), ] %>% nrow())
  
  # plot counts for every type of overlap
  plot <- barplot(height = counts,
                  names = c("Input", "Matches", "Equal", "Within", "Spans", "Other", "None"),
                  xlab = "Type of Overlap",
                  ylab = "Count",
                  ylim = c(0, 1.1*max(counts)))
  
  text(x = plot, y = counts, label = counts, pos = 3, cex = 1)
  
  return(overlap)
}
