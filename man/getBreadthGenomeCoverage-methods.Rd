% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DigitalPicoTools.R
\docType{methods}
\name{getBreadthGenomeCoverage}
\alias{getBreadthGenomeCoverage}
\title{Retrieves the fraction of the whole genomes covered by the Long fragments Reads  making the LFRset}
\usage{
getBreadthGenomeCoverage(object)
}
\arguments{
\item{x}{an object of class \code{\link{LFRset}}}
}
\value{
If a single well, a numeric value representing the breadth coverage of its LFRs. If multiple wells, the methods output a \code{\link{data.frame}} containing the breadth coverage per wells.
}
\description{
Retrieves the fraction of the whole genomes covered by the Long fragments Reads  making the LFRset
}
\examples{
#Loading an object of class (LFRset)
data(FiveWellsLFRset)
P009_LFR=getWellLFRset(FiveWellsLFRset,'P009')
print(getBreadthGenomeCoverage(P009_LFR))
# 4.63
print(getBreadthGenomeCoverage(FiveWellsLFRset))
# Wells BreadthCoverage
# 1  A005            7.22
# 2  A016            1.91
# 3  L006            3.69
# 4  L024            2.76
# 5  P009            4.63
}
\seealso{
\code{\link{getN50}}, \code{\link{getLFRStats}}, \code{\link{getLFRStatsPerWell}}, \code{\link{summary.LFRset}}
}

