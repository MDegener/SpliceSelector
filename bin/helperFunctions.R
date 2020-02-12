# Get gene annotation based on identifier (e.g. Ensembl gene id) using BioMart Ensembl database
getBiomaRtAnnotation <- function(identifierList, identifierType, species){
  library(biomaRt)
  
  ensembl <- useMart(biomart="ensembl",dataset=paste0(species,"_gene_ensembl"))
  
  geneSymbol <- if (species == "hsapiens") "hgnc_symbol" else if (species == "mmusculus") "mgi_symbol"
  
  att <- c("ensembl_gene_id", geneSymbol, "chromosome_name", 
           "start_position", "end_position", "strand", "description")
  
  ann <- getBM(attributes=att, filters=identifierType, values=identifierList, mart=ensembl)
  
  return(ann)
}

readLeafcutterResults <- function(effect_sizes, cluster_significance){
  library(tidyr)
  
  lc_eff <- read.csv(effect_sizes, sep = "\t")
  lc_eff <- separate(lc_eff, "intron", into = c("chr", "start", "end", "cluster_id"), sep = ":")
  
  lc_sig <- read.csv(cluster_significance, sep = "\t")
  lc_sig <- separate(lc_sig, "cluster", into = c("chr", "cluster_id"), sep = ":")
  lc_sig <- subset(lc_sig, select= -chr)
  
  lc_results <- merge(lc_eff, lc_sig, by = "cluster_id")
  
  return(lc_results)
}

liftOverCoordinates <- function(df, given_assembly, target_assembly){
  # TODO: add functionality to input multiple assemblies
  # TODO: use GRange object as input, not dataframe
  
  library(liftOver)
  
  target_assembly <- (paste0(toupper(substring(target_assembly,1,1)), substring(target_assembly,2)))
  
  chain <- paste0(given_assembly, "To", target_assembly, ".over.chain")
  chain_path <- paste0("lib/", chain)
  
  chain_url <- paste0("https://hgdownload.soe.ucsc.edu/goldenPath/", 
                      given_assembly,"/liftOver/", chain, ".gz")
  
  # TODO: Return error if chain file cannot be found
  if (!chain %in% list.files("./lib")){
    download.file(chain_url, destfile = paste0(chain_path, ".gz"), method = "wget", quiet = TRUE)
    R.utils::gunzip(paste0(chain_path, ".gz"), remove = TRUE)
  }
createUCSCtrack <- function(grObject, outFile, trackName, trackDescription){
  require(rtracklayer)

  # create output file
  file.create(outFile)
  
  # open connection to file
  fileConn<- file(outFile)
  
  # export UCSC track as .bed file and add trackLine
  export.ucsc(grObject, con = fileConn, subformat = "BED", 
              name = trackName,
              description = trackDescription,
              visibility = "2",
              priority = 1)
  
  # Probably not necessary because "export.ucsc" already closes file connection
  #close(fileConn)
}

