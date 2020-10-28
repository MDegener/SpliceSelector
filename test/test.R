######################
### SpliceSelector ###
######################

suppressPackageStartupMessages( library(stringr) )
suppressPackageStartupMessages( library(data.table) )

setwd("/home/maxd/SpliceSelector")

source("bin/helperFunctions.R")

print(str_extract(list.files("test/spliceEvents/"), "(.*)(?=_spliceEvents)"))
study <- "Nakamori2013"

spliceEvents <- read.csv(paste0("test/spliceEvents/",study,"_spliceEvents.csv"),
                         sep = "\t", stringsAsFactors = FALSE)

### lift over coordinates to most recent genome build ##############

if( substring(spliceEvents$Assembly[1],1,2) == "hg"){
  species <- "hsapiens"
  targetAssembly <- "hg38"
} else if ( substring(spliceEvents$Assembly[1],1,2) == "mm" ){
  species <- "mmusculus"
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

### create custom UCSC track ############

createUCSCtrack(grEvents,
                outFile = paste0("test/UCSC_customTracks/", study,"_UCSC_customTrack.bed"),
                trackName = paste0("DM1 Mis-Splicing (", study, ")"),
                trackDescription = paste0("Putative DM1-related splice events as reported by ", study))

### validate exon coordinates ###########

source("bin/getOverlap.R")
source("bin/plotOverlapOverview.R")

if (species == "hsapiens"){
  gtf <- import("lib/gencode.v32.annotation.gtf")

} else if (species == "mmusculus"){
  gtf <- import("lib/gencode.vM24.annotation.gtf")
}
 
exons <- subset(gtf[ which(gtf$type == "exon") ],
                select = c("gene_id", "transcript_id", "exon_id", "exon_number", "level"))

overlap <- getOverlap(grEvents, exons)

fwrite(overlap, paste0("test/tables/", study, "_overlap.csv"))

plotOverlapOverview(overlap, study, paste0(getwd(),"/test/plots"))

######## DEXSEQ ##########
# suppressPackageStartupMessages( library(DEXSeq) )
#
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
