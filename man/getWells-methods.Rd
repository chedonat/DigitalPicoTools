% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DigitalPicoTools.R
\docType{methods}
\name{getWells}
\alias{getWells}
\alias{getWells,LFRset-method}
\alias{getWells,getWells-method}
\title{Retrieving wells and fragments distribution}
\usage{
getWells(object)

\S4method{getWells}{LFRset}(object)
}
\arguments{
\item{x}{an object of class \code{\link{LFRset}}}
}
\value{
A vector having the different wells and their number of fragments
}
\description{
This function is a method of class \code{\link{LFRset}}. It retrieves the list of wells having fragments in the set of LFR as input and gives the numbe rof fragments for each well

Retrieving wells and fragments distribution
}
\examples{
#Loading an object of class (LFRset)
data(FiveWellsLFRset)
getWells(FiveWellsLFRset)
#A005  A016  L006  L024  P009
#43755 15192 23448 21892 24397
}
\seealso{
\code{\link{getWellLFRset}}, \code{\link{getLFRset}}, \code{\link{getLFRStats}}, \code{\link{getLFRStatsPerWell}}, \code{\link{plot.LFRset}}
}

