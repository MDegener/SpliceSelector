library(DEXSeq)

setwd("/home/maxd/mnt/xomics/SpliceSelector")

source("bin/helperFunctions.R")

splice_events <- read.csv("example/spliceEvents.csv", sep = ";", stringsAsFactors = FALSE)

hg38_coordinates <- liftOverCoordinates(splice_events, "hg18", "hg38")

splice_events$Start <- hg38_coordinates$start
splice_events$End <- hg38_coordinates$end
splice_events <- splice_events[, -6]
splice_events <- splice_events[ multi.mixedorder(splice_events$Chr, splice_events$Start), ]

load("/home/maxd/mnt/xomics/data/results/23012020035247/dexseq/dexseq.RData")
#load("example/Dexseq.RData")
dx_results <- as.data.frame(dxr) # Why only data of chromosome x?? (seqnames)

lc_results <- readLeafcutterData("example/leafcutter_ds_effect_sizes.txt",
                                 "example/leafcutter_ds_cluster_significance.txt")

rm(dxr, hg38_coordinates)

lc_found <- data.frame()
dx_found <- data.frame()

for( chr in unique(splice_events$Chr) ){

  chr_subset <- splice_events[ splice_events$Chr == chr, ]
  lc_chr <- lc_results[ lc_results$chr == chr, ]
  dx_chr <- dx_results[ dx_results$genomicData.seqnames == substring(chr, 4), ]
    
  for(i in 1:nrow(chr_subset)){

    min <- chr_subset$Start[i]
    max <- chr_subset$End[i]
  
    lc_start <- as.numeric(levels(lc_chr$start))[lc_chr$start]
    lc_end <- as.numeric(levels(lc_chr$end))[lc_chr$end]
    
    # start/end of exon of interest should lie within intron cluster
    lc_found <- rbind(lc_found, lc_chr[ which(min >= lc_start && max <= lc_end), ])

    dx_start <- dx_chr$genomicData.start
    dx_end <- dx_chr$genomicData.end
    
    # start/end of exon of interest should contain exon counting bin
    dx_found <- rbind(dx_found, dx_chr[ which(min <= dx_start && max >= dx_end ), ])
  }
}
