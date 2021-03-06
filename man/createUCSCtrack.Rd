% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/createUCSCtrack.R
\name{createUCSCtrack}
\alias{createUCSCtrack}
\title{UCSC custom track generator}
\usage{
createUCSCtrack(grObject, outFile, trackName, trackDescription)
}
\arguments{
\item{grObject}{GRanges object with genomic coordinates (see GenomicRanges documentation for more info)}

\item{outFile}{filepath and name of output UCSC track}

\item{trackName}{name of UCSC track}

\item{trackDescription}{description of UCSC track}
}
\value{
None
}
\description{
Creates a custom UCSC track from a genomic range object to view supplied coordinates in the UCSC genome browser
}
