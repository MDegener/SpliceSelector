library(DEXSeq)

setwd("/home/maxd/mnt/xomics/SpliceSelector")

source("bin/helperFunctions.R")

species <- "hsapiens"

if (species == "hsapiens"){
  spliceEvents <- read.csv("test/hs_spliceEvents.csv", sep = "\t", stringsAsFactors = FALSE)
  targetAssembly <- "hg38"

} else if (species == "mmusculus"){
  spliceEvents <- read.csv("test/mm_spliceEvents.csv", sep = "\t", stringsAsFactors = FALSE)
  targetAssembly <- "mm10"
}

# remove all commas in start/end position and convert to numeric
spliceEvents$Start <- as.numeric(gsub(",", "", spliceEvents$Start))
spliceEvents$End <- as.numeric(gsub(",", "", spliceEvents$End))

# create GRange object
grEvents <- GRanges(seqnames = Rle(spliceEvents$Chr),
                    ranges = IRanges(start = spliceEvents$Start, 
                                     end = spliceEvents$End),
                    strand = Rle(spliceEvents$Strand),
                    gene_symbol = spliceEvents$Symbol,
                    species = species)

# liftover coordinates for all unique genome assemblies in "spliceEvents"
for (givenAssembly in unique(spliceEvents$Assembly)){
  
  # check whether liftover is necessary
  if (givenAssembly != targetAssembly) {
    source("bin/liftOverGRange.R")
    
    liftedEventsTemp <- liftOverGRange(grEvents[which(spliceEvents$Assembly == givenAssembly)],
                                        givenAssembly, targetAssembly)
    
    grEvents[which(spliceEvents$Assembly == givenAssembly)] <- liftedEventsTemp
    
    rm(liftedEventsTemp, givenAssembly)
  }
}

### create custom UCSC track
# createUCSCtrack(grEvents, 
#                 outFile = paste0("test/", target_assembly,"_UCSC_customTrack.bed"),
#                 trackName = "DM1 Mis-Splicing",
#                 trackDescription = "Putative DM1-related splice events")


######## DEXSEQ ##########
# # load dexseq results
# load(paste0("test/dexseq_", species, ".RData"))
# 
# # remove "chr" prefix
# if( !any(grep("chr", levels(seqnames(dxr$genomicData)))) ){
#   seqlevels(grEvents) <- gsub("chr", "", seqlevels(grEvents))
# }
# 
# # delete all NA rows
# dxr <- dxr[ -which(is.na(dxr$stat)), ]
# 
# # get DEXSeq GRange object
# grDEX = dxr$genomicData
# 
# dx_overlaps <- findOverlaps(query = grEvents,
#                             subject = grDEX,
#                             type = "within")
# 
# # TEMP
# dxr <- as.data.frame(dxr)
# rownames(dxr) <- make.names(rownames(dxr), unique = TRUE)
# 
# dx_results <- as.data.frame(dxr[ to(dx_overlaps), ])

######## LEAFCUTTER ##########
# # load leafcutter results
# lc_results <- readLeafcutterResults("test/leafcutter_ds_effect_sizes.txt",
#                                     "test/leafcutter_ds_cluster_significance.txt")
# 
# # create GRange object based on leafcutter results
# # TODO: add strand info to leafcutter GRange?
# grLC <- GRanges(seqnames = Rle(gsub("chr", "", lc_results$chr)),
#                      ranges = IRanges(start = as.numeric(lc_results$start),
#                                       end = as.numeric(lc_results$end)))
# 
# lc_overlaps <- findOverlaps(query = grEvents,
#                             subject = grLC, 
#                             type = "within")

# Features to add:
# - visualize overlaps (venn diagram)
# - add gene ontology analysis

# TODO: create GRange object when liftOver is not necessary
#spliceEvents <- read.csv("test/mm_spliceEvents.csv", sep = ";", stringsAsFactors = FALSE)

# liftover coordinates and return GRange object
#grEvents <- liftOverCoordinates(spliceEvents, "mm8", "mm10")
