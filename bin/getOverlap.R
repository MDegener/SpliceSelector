getOverlap <- function(grObject, exonAnnotation){
  
  require(rtracklayer, dplyr)

  # find overlap
  exactOverlap <- findOverlaps(grObject, exonAnnotation, type = "equal")
  withinOverlap <- findOverlaps(grObject, exonAnnotation, type = "within")
  noOverlap <- grObject[ grObject %outside% exonAnnotation ]
  
  # remove redundant exact matches from withinOverlap
  withinOverlap <- withinOverlap[ -which(withinOverlap %in% exactOverlap) ]
  
  # write overlap to dataframe and save type of overlap
  overlap <- data.frame(overlap = "equal",
                        query = grObject[ from(exactOverlap)], 
                        annotation = exonAnnotation[ to(exactOverlap) ],
                        stringsAsFactors = FALSE)
  
  overlap <- rbind(overlap, data.frame(overlap = "within",
                                       query = grObject[ from(withinOverlap)], 
                                       annotation = exonAnnotation[ to(withinOverlap) ],
                                       stringsAsFactors = FALSE))
  
  overlap <- bind_rows(overlap, data.frame(overlap = "none",
                                           query = noOverlap,
                                           stringsAsFactors = FALSE))
  

  counts = c(exactOverlap %>% from() %>% unique() %>% length(),
                 withinOverlap %>% from() %>% unique() %>% length(), 
                 length(noOverlap))
  
  counts = c(counts, length(grObject) - sum(counts))

  # plot counts for every type of overlap
  plot <- barplot(height = counts,
                  names = c("Exact", "Within", "None", "Other"),
                  xlab = "Type of Overlap",
                  ylab = "Count",
                  ylim = c(0, 1.1*max(matchCount)))
  
  text(x = plot, y = counts, label = counts, pos = 3, cex = 1)
  
  return(overlap)
}
