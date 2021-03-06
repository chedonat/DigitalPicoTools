% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DigitalPicoTools.R
\name{plot.VariantAlleleInfo}
\alias{plot.VariantAlleleInfo}
\title{Plotting the variant allele information}
\usage{
\method{plot}{VariantAlleleInfo}(object, region = NULL, samplingRatio = 1,
  value = "WellsFraction", main = "")
}
\arguments{
\item{region}{Region to be plotted.  If NULL then the whole region covered by the variants will be considered}

\item{samplingRatio}{Ratio of points to plot. Only 1 each samplingRatio points will be plotted.}

\item{main}{The title of the plot. default: None}

\item{x}{an  \code{\link{VariantAlleleInfo}} object containing the list of variant and their genomic location and allele count information}

\item{Value}{Information to plot. Can be one of WellsFraction, WellsCounts, ReadsFraction, ReadsCounts.}
}
\value{
''
}
\description{
This is a generic function to plot different information on the variant allele counts. The information to plot is set by the parameter \emph{Value} which can be eithe rof the following :
\describe{
  \item{WellsFraction}{Wells allele fraction of the variant}
  \item{WellsCounts}{Wells allele counts of the variants }
  \item{ReadsFraction}{Reads allele fraction of the variant}
  \item{ReadsCounts}{Reads allele counts of the variants }
}
}
\details{
In both types, the paramater \emph{Value} allow to specifie which variable will be represented on the y-axis. It can be  'Coverage' for the LFR Average Coverage and 'Length' for the LFR lengths.
}
\examples{

#Example 1:
vcffile<-system.file('extdata','LFR_11152_Tissue_chr22_1000variants.vcf.gz',package='DigitalPicoTools')
AlleleInfo_df = getVariantAlleleInfo(vcffile)
plot(AlleleInfo_df,'chr22')

}
\seealso{
\code{\link{getVariantAlleleInfo}}, \code{\link{getVariantCoverageTable}}, \code{\link{VariantAlleleInfo}}
}

