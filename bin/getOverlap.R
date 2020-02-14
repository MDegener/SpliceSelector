# ADD DESCRIPTION

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
      
      # remove redundant exact matches from withinOverlap, spanOverlap and anyOverlap
      withinOverlap <- withinOverlap[ -which(withinOverlap %in% exactOverlap) ]
      spanOverlap <- spanOverlap[ -which(spanOverlap %in% exactOverlap) ] 
      anyOverlap <- anyOverlap[ -which(anyOverlap %in% exactOverlap) ]
      
      # save matches in dataframe and annotate with its specific type
      overlap <- bind_rows(overlap, data.frame(overlap = "equal",
                                               query.index = from(exactOverlap),
                                               query = grObject[ from(exactOverlap) ],
                                               annotation.index = to(exactOverlap),
                                               annotation = exonAnnotation[ to(exactOverlap) ],
                                               stringsAsFactors = FALSE))
    }
    
    if( length(withinOverlap) != 0 ){
      
      # remove redundant within matches from anyOverlap
      anyOverlap <- anyOverlap[ -which(anyOverlap %in% withinOverlap) ]
      
      overlap <- bind_rows(overlap, data.frame(overlap = "within",
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
      spanOverlap <- data.frame(overlap = "",
                                query.index = from(spanOverlap),
                                query = grObject[ from(spanOverlap) ],
                                annotation.index = to(spanOverlap),
                                annotation = exonAnnotation[ to(spanOverlap) ],
                                stringsAsFactors = FALSE)
      
      # distinguish coordinates that span one or multiple exons
      spanMatches <- spanOverlap %>% group_by(query.index) %>% summarise(n_distinct(annotation.exon_id))
      spanOverlap$overlap <- map_chr(spanOverlap$query.index, 
                                     function( queryIndex ) {
                                       # check if number of span matches is equal to 1 for a given query
                                       if ( spanMatches[ queryIndex == spanMatches$query.index, 2] == 1){
                                         "spans one" 
                                       } else {
                                         "spans multiple"
                                       }
                                     })

      overlap <- bind_rows(spanOverlap)
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
  
  # get number of exon matches per input coordinates
  overlap %>%
    group_by(query.index) %>%
    summarise(n_distinct(annotation.exon_id)) -> exonMatches

  # get unique matches of input coordinates to an exon
  uniqueMatches <- sum( exonMatches[, 2] == 1)
  
  # ADD COMMENT
  matchCounts <- c(overlap$query.index %>% unique() %>% length(),
                   overlap[which(overlap$overlap != "none"), ] %>% nrow(),
                   uniqueMatches,
                   overlap[which(overlap$overlap == "none"), ] %>% nrow())
  
  # ADD COMMENT
  plot <- barplot(height = matchCounts,
                  names = c("Input", "All Matches", "Unique", "No Match"),
                  ylab = "Count",
                  ylim = c(0, 1.1*max(matchCounts)))
  
  text(x = plot, y = matchCounts, label = matchCounts, pos = 3, cex = 1)
  
  # get counts for all overlap types
  overlapCounts <- c(overlap[which(overlap$overlap == "equal"), ] %>% nrow(),
                     overlap[which(overlap$overlap == "within"), ] %>% nrow(),
                     overlap[which(overlap$overlap == "spans one"), ] %>% nrow(),
                     overlap[which(overlap$overlap == "spans multiple"), ] %>% nrow(),
                     overlap[which(overlap$overlap == "other"), ] %>% nrow())

  # plot counts for every type of overlap
  plot <- barplot(height = overlapCounts,
                  names = c("Equal", "Within", "Spans One", "Spans Multiple", "Other"),
                  xlab = "Type of Overlap",
                  ylab = "Count",
                  ylim = c(0, 1.1*max(overlapCounts)))
  
  text(x = plot, y = overlapCounts, label = overlapCounts, pos = 3, cex = 1)
  
  return(overlap)
}
