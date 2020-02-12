### written by @DaanJG98

library(dplyr)
library(stringr)
library(DEXSeq)

setwd("/home/maxd/mnt/xomics/SpliceSelector")

# Function to create and save plot for DEXSeq results
createDEXSeqplot <- function(x, output){
  dir.create(file.path(output), recursive=T)
  png(paste0(output,x$Gene_name,".png"), width = 1000, height = 958)
  plotDEXSeq(dxr, x$groupID, fitExpToVar = "tissue",
             splicing = T, expression = T, norCounts = T, legend = T, cex=1.3, cex.axis=1.2, lwd=1, FDR=0.05)
  dev.off()
}

# Gene annotation used to give a gene symbol to exon bins
gene_annot <- paste0(getwd(), "/lib/Homo_sapiens.GRCh38.95_nodupes_genid.txt")
gene_annot <- read.table(gene_annot, header=T, sep="\t")

# Load DEXSeq results
dexseq_result <- paste0(getwd(), "/test/dexseq_hsapiens.RData")
load(dexseq_result)

# convert to simple data frame for easier handling
dxr_slim <- as.data.frame(dxr)

# TEMP: Remove .n from transcript id
dxr_slim$geneID <- str_replace(dxr_slim$groupID, "(\\.[0-9]*)", "")

# combine with annotation
dxr_annot <- merge(dxr_slim, gene_annot, by.x="geneID", by.y="Gene_id", all.x=T)
# remove rows with NA values
dxr_red <- na.omit(dxr_annot)
# filter for abs log2 FC > 1.0 and adj p-value < 0.05
dxr_red <- subset(dxr_red, dxr_red$padj < 0.05 & abs(dxr_red %>% dplyr::select(starts_with("log2"))) > 1.0)

# create plots
if (dim(dxr_red)[1]!=0){
    apply(dxr_red, 1, FUN=createDEXSeqplot, output=paste0(getwd(),"/test/plots/"))
}
# save results
write.csv(dxr_red, file=paste0(getwd(),"/test/dexseq_significant.csv"))
