library(DEXSeq)
library(GenomicRanges)

setwd("/home/maxd/mnt/xomics/SpliceSelector")

source("bin/helperFunctions.R")

species <- "hsapiens"

if (species == "hsapiens"){
  splice_events <- read.csv("test/hs_spliceEvents.csv", sep = "\t", stringsAsFactors = FALSE)
  target_assembly <- "hg38"
  
} else if (species == "mmusculus"){
  splice_events <- read.csv("test/mm_spliceEvents.csv", sep = "\t", stringsAsFactors = FALSE)
  target_assembly <- "mm10"
}

# remove all commas in start/end position and convert to numeric
splice_events$Start <- as.numeric(gsub(",", "", splice_events$Start))
splice_events$End <- as.numeric(gsub(",", "", splice_events$End))

# create GRange object and liftover coordinates if necessary
if ( any(splice_events$Assembly != target_assembly)) {
  events_GRange <- liftOverCoordinates(splice_events, unique(splice_events$Assembly), target_assembly)
  
} else {
  events_GRange <- GRanges(seqnames = Rle(df$Chr),
                           ranges = IRanges(start = splice_events$Start, end = splice_events$End),
                           strand = Rle(splice_events$Strand),
                           gene_symbol = splice_events$Symbol)
  
  genome(events_GRange) <- unique(splice_events$Assembly)
}

# load dexseq results
load(paste0("test/dexseq_", species, ".RData"))
### create custom UCSC track
# createUCSCtrack(grEvents, 
#                 outFile = paste0("test/", target_assembly,"_UCSC_customTrack.bed"),
#                 trackName = "DM1 Mis-Splicing",
#                 trackDescription = "Putative DM1-related splice events")

# remove "chr" prefix
if( !any(grep("chr", levels(seqnames(dxr$genomicData)))) ){
  seqlevels(events_GRange) <- gsub("chr", "", seqlevels(events_GRange))
}

# delete all NA rows
dxr <- dxr[ -which(is.na(dxr$stat)), ]

# get DEXSeq GRange object
dx_GRange = dxr$genomicData

dx_overlaps <- findOverlaps(query = events_GRange,
                            subject = dx_GRange,
                            type = "within")

dxr <- as.data.frame(dxr)
rownames(dxr) <- make.names(rownames(dxr), unique = TRUE)
dx_results <- as.data.frame(dxr[ to(dx_overlaps), ])

# # load leafcutter results
# lc_results <- readLeafcutterResults("test/leafcutter_ds_effect_sizes.txt",
#                                     "test/leafcutter_ds_cluster_significance.txt")
# 
# # create GRange object based on leafcutter results
# # TODO: add strand info to leafcutter GRange?
# lc_GRange <- GRanges(seqnames = Rle(gsub("chr", "", lc_results$chr)),
#                      ranges = IRanges(start = as.numeric(lc_results$start),
#                                       end = as.numeric(lc_results$end)))
# 
# lc_overlaps <- findOverlaps(query = events_GRange,
#                             subject = lc_GRange, 
#                             type = "within")

# Features to add:
# - visualize overlaps (venn diagram)
# - add gene ontology analysis

# TODO: create GRange object when liftOver is not necessary
#splice_events <- read.csv("test/mm_spliceEvents.csv", sep = ";", stringsAsFactors = FALSE)

# liftover coordinates and return GRange object
#events_GRange <- liftOverCoordinates(splice_events, "mm8", "mm10")
