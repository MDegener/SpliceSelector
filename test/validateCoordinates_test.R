setwd("/home/maxd/mnt/xomics/SpliceSelector")

library(rtracklayer)

source("bin/helperFunctions.R")
source("bin/liftOverGRange.R")
source("bin/getOverlap.R")

# m24_gtf <- import("lib/gencode.vM24.annotation.gtf")
# m24_exons <- subset(m24_gtf[ which(m24_gtf$type == "exon") ],
#                      select = c("gene_id", "transcript_id", "exon_id", "exon_number", "level"))
# 
# splice_events <- read.delim("test/mm_spliceEvents.csv", sep = "\t", stringsAsFactors = F)
# 
# mm8_gr <- GRanges(seqnames = Rle(splice_events$Chr),
#                    ranges = IRanges(start = splice_events$Start, end = splice_events$End),
#                    strand = Rle(splice_events$Strand),
#                    gene_symbol = splice_events$Symbol
# )
# 
# mm10_gr <- liftOverGRange(mm8_gr, "mm8", "mm10")
# 
# mm_overlap <- getOverlap(mm10_gr, m24_exons)

### create custom UCSC track
# createUCSCtrack(mm10_gr,
#                 outFile = paste0("test/mm10_UCSC_customTrack.bed"),
#                 trackName = "DM1 Mis-Splicing",
#                 trackDescription = "Putative DM1-related splice events")

gc32_gtf <- import("lib/gencode.v32.annotation.gtf")
gc32_exons <- subset(gc32_gtf[ which(gc32_gtf$type == "exon") ],
                     select = c("gene_id", "transcript_id", "exon_id", "exon_number", "level"))

splice_events <- read.delim("test/hs_spliceEvents.csv", sep = "\t", stringsAsFactors = F)

hg18_gr <- GRanges(seqnames = Rle(splice_events$Chr),
                   ranges = IRanges(start = splice_events$Start, end = splice_events$End),
                   strand = Rle(splice_events$Strand),
                   gene_symbol = splice_events$Symbol
)

hg38_gr <- liftOverGRange(hg18_gr, "hg18", "hg38")

hs_overlap <- getOverlap(hg38_gr, gc32_exons)

# write.table(overlap, "test/exonOverlaps.txt", row.names = F, col.names = T, quote = F, sep = "\t")
