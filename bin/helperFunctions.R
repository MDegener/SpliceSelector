
# not written by me 
# ref: https://stackoverflow.com/questions/20396582/order-a-mixed-vector-numbers-with-letters
multi.mixedorder <- function(..., na.last = TRUE, decreasing = FALSE){
  do.call(order, c(
    lapply(list(...), function(l){
      if(is.character(l)){
        factor(l, levels=gtools::mixedsort(unique(l)))
      } else {
        factor(as.character(l), levels=gtools::mixedsort(levels(l)))
      }
    }),
    list(na.last = na.last, decreasing = decreasing)
  ))
}

### USE SEPARATE FUNCTION INSTEAD (tidyr)
expandColumn <- function(df, col, sep, names){
  library(stringr)
  
  col_expanded <- data.frame(str_split(df[, col], sep, simplify = TRUE))
  
  colnames(col_expanded) <- names
  
  df <- cbind(df[, -match(col,names(df))], col_expanded)
  
  return(df)
}

readLeafcutterData <- function(effect_sizes, cluster_significance){
  
  lc_eff <- read.csv(effect_sizes, sep = "\t")
  lc_sig <- read.csv(cluster_significance, sep = "\t")
  
  lc_eff <- expandColumn(lc_eff, "intron", ":", c("chr", "start", "end", "cluster_id"))
  
  lc_sig <- expandColumn(lc_sig, "cluster", ":", c("chr", "cluster_id"))
  lc_sig <- subset(lc_sig, select=-c(chr))
  
  lc_results <- merge(lc_eff, lc_sig, by = "cluster_id")
  
  return(lc_results)
}

liftOverCoordinates <- function(df, given_assembly, target_assembly){
  library(liftOver)
  
  target_assembly <- (paste0(toupper(substring(target_assembly,1,1)), substring(target_assembly,2)))
  
  chain <- paste0(given_assembly, "To", target_assembly, ".over.chain")
  chain_path <- paste0("lib/", chain)
  
  chain_url <- paste0("https://hgdownload.soe.ucsc.edu/goldenPath/", 
                      given_assembly,"/liftOver/", chain, ".gz")
  
  if (!chain %in% list.files("./lib")){
    download.file(chain_url, destfile = paste0(chain_path, ".gz"), method = "wget", quiet = TRUE)
    R.utils::gunzip(paste0(chain_path, ".gz"), remove = TRUE)
  }
  
  chain <- import.chain(chain_path)
  
  grObject <- GRanges(seqnames = df$Chr, ranges = IRanges(start = df$Start, end = df$End))
  
  return(as.data.frame(liftOver(grObject, chain)))
}

