# ADD DESCRIPTION

plotOverlapOverview <- function(overlap, plotPrefix, outDir){
  
  library(dplyr)
  
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
  
  # plot overview on how coordinates are matched to exons
  jpeg(paste0(outDir,"/", plotPrefix, "_matchCounts.jpeg"),
              width = 1000, height = 1000, pointsize = 24)
  
  plot <- barplot(height = matchCounts,
                  names = c("Input", "All Matches", "Unique", "No Match"),
                  ylab = "Count",
                  ylim = c(0, 1.1*max(matchCounts)))
  
  text(x = plot, y = matchCounts, label = matchCounts, pos = 3, cex = 1)
  dev.off()
  
  # get counts for all overlap types
  overlapCounts <- c(overlap[which(overlap$overlap == "equal"), ] %>% nrow(),
                     overlap[which(overlap$overlap == "within"), ] %>% nrow(),
                     overlap[which(overlap$overlap == "spans one"), ] %>% nrow(),
                     overlap[which(overlap$overlap == "spans multiple"), ] %>% nrow(),
                     overlap[which(overlap$overlap == "other"), ] %>% nrow())
  
  # plot counts for every type of overlap
  jpeg(paste0(outDir,"/", plotPrefix, "_overlapCounts.jpeg"),
       width = 1000, height = 1000, pointsize = 24)
  
  plot <- barplot(height = overlapCounts,
                  names = c("Equal", "Within", "Spans One", "Spans Multiple", "Other"),
                  xlab = "Type of Overlap",
                  ylab = "Count",
                  ylim = c(0, 1.1*max(overlapCounts)))
  
  text(x = plot, y = overlapCounts, label = overlapCounts, pos = 3, cex = 1)
  dev.off()
}