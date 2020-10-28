# SpliceSelector Tool

The SpliceSelector package faclitates a targeted analysis of alternatively spliced candidate exons in any RNAseq dataset. The required input will be a list of genomic coordinates (chr:start-end), an exon annotation file (e.g. [GENCODE](https://www.gencodegenes.org/human/release_32.html)) and the results from common splicing analysis tools. Results from [MISO](https://miso.readthedocs.io/en/fastmiso/), [DEXSeq](http://bioconductor.org/packages/release/bioc/html/DEXSeq.html), [leafcutter](https://github.com/davidaknowles/leafcutter) and [MAJIQ](https://majiq.biociphers.org/) will be supported.

## Development state

Currently, only the following helper functions have been implemented:
  
**1. liftOverGRange.R**
  
* wrapper function for [UCSC's liftOver tool](https://genome.ucsc.edu/cgi-bin/hgLiftOver) 
* **usage:** liftOverGRange(grObject, givenAssembly, targetAssembly, libPath)
* converts coordinates in genomic range object from the given assembly to the target assembly
* liftOver chain is stored under libPath

**2. getOverlap.R**

* compares query coordinates with coordinates from a .gtf/.gff exon annotation file (e.g. GENCODE)
* creates a table listing all overlapping coordinates with their corresponding overlap type 
* **usage:** getOverlap(grObject, exonAnnotation)
* grObject is a genomic range object that can be generated with [rtracklayer](https://www.bioconductor.org/packages/release/bioc/html/rtracklayer.html) and contains the query coordinates
* overlap types:
    * "equal" - perfect match of start and end positions
    * "within" - query coordinates that fall within annotation coordinates of an exon
    * "spans one" - query coordinates that span annotation coordinates of one exon
    * "spans multiple" - query coordinates that span annotation coordinates of multiple exons
    * "none" - query coordinates that do not match to any exon in annotation coordinates
    * "other" - any type of overlap that is not covered by the other categories

**3. plotOverlapOverview.R**

* plots results of getOverlap.R script in two barplots
* **usage:** plotOverlapOverview(overlap, plotPrefix, outDir)

**4. getBiomaRtAnnotation** 

* wrapper function for querying BioMart Ensembl database to extract gene annotation (e.g. coordinates, description)
* **usage:** getBiomaRtAnnotation(identifierList, identifierType, species)

**5. createUCSCtrack**

* creates a custom UCSC track from a genomic range object to view supplied coordinates in the [UCSC genome browser](https://genome.ucsc.edu/)
* **usage:** createUCSCtrack(grObject, outFile, trackName, trackDescription)