#########################################################
### DigitalPicoTools : files Organization data.R - file to generate the Roxygen help documents for the data DigitalPicoTools.R - this file is the
###############################################
############# main file of the package VCFtools.R - This file contains the function for vcf files.




# Roxigen help for the package (top help)

#' DigitalPicoTools : A set of tools  for Digital analysis of whole genome sequencing of picogram quantities of DNA
#'
#' The main functions provided by the tools are : \code{\link{getVariantAlleleInfo}}, \code{\link{getLFRset}} and \code{\link{getPhasingInfo}}
#'
#'
#' The package include experimental data retrieved from a study conducted in the Ovarian Cancer Lab of the University of Oxford.
#'
#' For more detailed information on usage, see the package vignette, by typing
#' \code{vignette('DigitalPicoTools')}. All support questions should be emailed to the authors.
#'
#' @references
#'
#' DigitalPicoTools references:
#'
#' OncoPhase: A package for computing Somatic Mutation cellular Prevalence in cancer using haplotype phasing. Bioinformatics 2016. Submitted
#'
#' DigitalPicoTools reference (Paper):
#'
#'  Ovarian cancer haplotype sequencing reveals ubiquitous SOX2 overexpression in the premalignant fallopian tube epithelium
#'
#' @author Donatien Chedom-Fotso, Ahmed Ahmed, Christopher Yau,
#'
#' @docType package
#' @name DigitalPicoTools
#' @aliases DigitalPicoTools-package
#' @keywords package
#' @import ggplot2 Rsamtools S4Vectors VariantAnnotation
NULL




#' @export
chrom_list <- paste("chr", c(1:22, "X", "Y", "M"), sep = "")
# Chromosomes sizes
#' @export
hg19_chrsize <- list(chr1 = 249250621, chr2 = 243199373, chr3 = 198022430, chr4 = 191154276, chr5 = 180915260, chr6 = 171115067, chr7 = 159138663,
    chr8 = 146364022, chr9 = 141213431, chr10 = 135534747, chr11 = 135006516, chr12 = 133851895, chr13 = 115169878, chr14 = 107349540, chr15 = 102531392,
    chr16 = 90354753, chr17 = 81195210, chr18 = 78077248, chr19 = 59128983, chr20 = 63025520, chr21 = 48129895, chr22 = 51304566, chrX = 155270560,
    chrY = 59373566, chrM = 16571)
#' @export
chrom_offset <- c(0)
# chrom_middle<-c(0 + hg19_chrsize[chrom_list[1]][[1]]/2)


for (chr in 2:length(hg19_chrsize)) {
    chrom_offset <- c(chrom_offset, chrom_offset[chr - 1] + as.numeric(hg19_chrsize[chr - 1]))
    # chrom_middle<-c(chrom_middle, chrom_offset[chr]+hg19_chrsize[chrom_list[chr]][[1]]/2)
}


# Welcome message

.onAttach <- function(libname, pkgname) {
    packageStartupMessage("\n************\nWelcome to DigitalPicoTools package\n************")
}

# Set up some custom options

.onLoad <- function(libname, pkgname) {
    op <- options()
    op.DigitalPicoTools <- list(DigitalPicoTools.path = "~/R", DigitalPicoTools.install.args = "", DigitalPicoTools.name = "Donatien Chedom-Fotso",
        DigitalPicoTools.desc.author = "\"Donatien Chedom-Fotso <donatien.chedom@gmail.com>\n    [aut, cre]\"", DigitalPicoTools.desc.license = "GPL-2",
        DigitalPicoTools.desc.suggests = NULL, DigitalPicoTools.desc = list(), DigitalPicoTools.PatientName = "NoName", DigitalPicoTools.WD = getwd())
    toset <- !(names(op.DigitalPicoTools) %in% names(op))
    if (any(toset))
        options(op.DigitalPicoTools[toset])

    invisible()
}


#' The LFRdataframe class
#'
#' An S4 base class. All other LFR class containing a data frame Table as main object inherit from this base class.
#' @slot Table Object of class \code{\link{data.frame}}, main container of the \code{\link{LFRdataframe}} class
#' @seealso \code{\link{LFRset}}, \code{\link{SampleCoverage}}, \code{\link{VariantAlleleInfo}}
#' @export LFRdataframe
#' @exportClass LFRdataframe
LFRdataframe <- setClass("LFRdataframe", slots = c(Table = "data.frame"))

#' @export
as.data.frame.LFRdataframe <- function(object) {
    object@Table
}

# head <- function(x, ...) UseMethod('LFRdataframe')
#' @export
head.LFRdataframe <- function(object) {
    LFRdataframe(Table = head(object@Table))
}
# setGeneric('heaad') setMethod('head',signature(object='LFRdataframe'),function(object) { head.LFRdataframe(object)})
#' @export
tail.LFRdataframe <- function(object) {
    LFRdataframe(Table = tail(object@Table))
}

#' @export
print.LFRdataframe <- function(object) {
    cat("\n\t An object of class LFRdataframe with ", nrow(object), " rows and ", ncol(object), " columns\n")
    if (nrow(as.data.frame(object)) < 20) {
        print(object@Table)
    } else {
        print(head(object@Table))
        cat("\n\t\t.........\n")
        print(tail(object@Table))
    }
}






#' The WellSample class
#'
#' An S4 class containing the main information (Well ID, well Number, Bam file) of a well sample.
#'
#' @slot Well_ID Object of class \code{\link{character}} representing the ID or the name assigned to the well
#' @slot Well_Number Object of class \code{\link{numeric}} representing the number of the well
#' @slot Bamfilename Object of class \code{\link{character}} representing the file name of the well bam file.
#' @examples
#'  #Creating an instance of the class
#'  myWell = WellSample(Well_ID='W384',Well_Number=384,BAMfilename='extdata/P009.bam')
#'  print(myWell)
#'  # An S4 object of class WellSample
#'  #
#'  # Well_ID :  W384
#'  # Well_Number :  384
#'  # BAMfilename :  extdata/P009.bam
#'  #
#' @seealso \code{\link{initsWellSamples}}
#' @rdname WellSample
#' @name WellSample
#' @rdname WellSample
#' @aliases WellSample-class
#' @export WellSample
#' @exportClass WellSample
WellSample <- setClass("WellSample", slots = c(Well_ID = "character", Well_Number = "numeric", BAMfilename = "character"), prototype = list(Well_ID = "",
    Well_Number = 0, BAMfilename = ""), validity = function(object) {
    if (!file.exists(object@BAMfilename))
        warnings(paste(" The BAM file of the created well : ", object@Well_ID, " do not exists"))
})


#' Initialising a set of Well samples
#'
#' This function take as input a vector containing a series of well ID and output of list of \code{\link{WellSample}} objects corresponding to each well.
#'
#' @param wellSamples_list a vector  of well ID.
#' @param Bamprefix a character string representing a prefix to add to the well ID to match the bam file name in the system file. Should contains the directories making path to the file. Putting the absolute path to the file is recommaeded.  default value : empty string.
#' @param Bamsuffix a character string representing a suffix to add after the well_ID in order to match the Bam file name on the file system. This suffux should not include the extemsion '.bam'. Defaults value : empty string
#' @param unexists.action Define the action to perform when a bam file do not exists. Can be one of the following:
#' \describe{
#'   \item{\emph{unexists.pass}}{No particular action, consider the provided bam file and continue.}
#'   \item{\emph{unexists.fail}}{Interupt the run and generate an error}
#'   \item{\emph{unexists.exclude}}{Do not process the sample, remove it from the list of well and continue}
#' }
#' The default value is set at \emph{unexists.fail}
#' @return A list of objects of class \code{\link{WellSample}}
#' @examples
#' wellsID_file=system.file('extdata','wells_id_tissue.txt',package='DigitalPicoTools')
#'  wellsID_list<-unlist(read.table(wellsID_file))
#'  WellSample_list=initsWellSamples(wellsID_list,
#'     BAMprefix=paste(path.package('DigitalPicoTools'),'extdata/',sep='/'))
#'  print(WellSample_list[3])
#'  # An S4 object of class WellSample
#'  #
#'  # Well_ID :  L006
#'  # Well_Number :  3
#'  # BAMfilename :  extdata/L006.bam
#'  #
#' @seealso \code{\link{WellSample}}
#' @export
initsWellSamples <- function(WellSamples_list, BAMprefix = NULL, BAMsuffix = NULL, unexists.action = "unexists.fail") {


    WellSampleList <- list()
    action.type = c("unexists.pass", "unexists.fail", "unexists.exclude")
    if (!(unexists.action %in% c("unexists.pass", "unexists.fail", "unexists.exclude")))
        stop("\nThe parameter unexists.action should be one of the following : unexists.pass unexists.fail unexists.exclude")

    for (iWellSample in 1:length(WellSamples_list)) {

        wellid = as.character(WellSamples_list[iWellSample])

        if (!is.null(BAMprefix) || !is.null(BAMsuffix))
            BAMfilename = paste(BAMprefix, wellid, BAMsuffix, ".bam", sep = "")

        if (!file.exists(BAMfilename)) {
            if (unexists.action == "unexists.fail") {
                stop(paste(" \n\n The BAM file ", BAMfilename, " for well ", wellid, " do not exists \n\t check your parameters or set unexists.action to unexists.exclude",
                  sep = ""))
            } else if (unexists.action == "unexists.exclude") {
                next
            }

        }

        newWellSample = WellSample(Well_ID = wellid, Well_Number = iWellSample, BAMfilename = BAMfilename)
        WellSampleList <- append(WellSampleList, newWellSample)
    }

    return(WellSampleList)
}


setGeneric("setProperty", def = function(object, xname, xvalue) {
    standardGeneric("setProperty")
})

setMethod(f = "setProperty", signature = "WellSample", definition = function(object, xname, xvalue) {
    slot(object, xname) = xvalue
    return(object)
})

#' @export
print.WellSample <- function(object) {
    cat("\n An S4 object of class WellSample\n")
    cat("\n\t Well_ID : ", object@Well_ID)
    cat("\n\t Well_Number : ", object@Well_Number)
    cat("\n\t BAMfilename : ", object@BAMfilename)
    cat("\n\n")

}






#' The LFRset class
#'
#' An S4 class containing the set of Long Fragments Reads (LFR) retrieved from an LFRSeq experiment.
#'
#' @slot Table  A data frame serving as a containers of the LFRs. Each LFR is represented by the following fields :
#' \describe{
#'   \item{LFRName}{The name of the LFR. By default it is a combination of the well's name the LFR is originating from and the genomic location of the LFR (Chromosome, atart position and end position)}
#'   \item{Chrom}{The chromosome of the LFR}
#'   \item{Start}{The starting position  of the LFR}
#'   \item{End}{The ending position of the LFR}
#'   \item{Length}{The length of the LFR}
#'   \item{AvgCoverage}{The average Coverage of the LFR}
#'   \item{NbReads}{The number of reads making up the LFR}
#'   \item{Well_ID}{The name or ID of the well the LFR is originating from}
#'   \item{Well_Number}{The number of the well}
#' }
#' @rdname LFRset
#' @note This class is one of the main class of the package
#' @name LFRset
#' @rdname LFRset
#' @aliases LFRset-class
#' @seealso \code{\link{getLFRset}}, \code{\link{getWells}}, \code{\link{getLFRStats}}, \code{\link{plotLFRDistribution}}, \code{\link{plot.LFRset}}
#' @export LFRset
#' @exportClass LFRset
LFRset <- setClass("LFRset", slots = character(0), prototype = list(Table = data.frame(nrow = 0, ncol = 10)), validity = function(object) {
    colname_table = colnames(object@Table)
    compulsory_columns = c("LFR_name", "Chrom", "Start", "End", "Length", "AvgCoverage", "Nbreads", "Well_ID")
    missing_columns = setdiff(compulsory_columns, colname_table)
    if (length(missing_columns) > 0)
        stop("\n\n\t Non valid object FragmengtsInfo, the following columns are missing from the Table : ", missing_columns)
}, contains = "LFRdataframe")


#' Retrieving wells and fragments distribution
#'
#' This function is a method of class \code{\link{LFRset}}. It retrieves the list of wells having fragments in the set of LFR as input and gives the numbe rof fragments for each well
#'
#' @name getWells
#' @rdname getWells-methods
#' @param x an object of class \code{\link{LFRset}}
#' @return A vector having the different wells and their number of fragments
#' @seealso \code{\link{getWellLFRset}}, \code{\link{getLFRset}}, \code{\link{getLFRStats}}, \code{\link{getLFRStatsPerWell}}, \code{\link{plot.LFRset}}
#' @examples
#' #Loading an object of class (LFRset)
#' data(FiveWellsLFRset)
#' getWells(FiveWellsLFRset)
#' #A005  A016  L006  L024  P009
#' #43755 15192 23448 21892 24397
#' @docType methods
#' @export
setGeneric("getWells", def = function(object) {
    standardGeneric("getWells")
})


#' Retrieving wells and fragments distribution
#'
#' @rdname getWells-methods
#' @aliases getWells,getWells-method
#' @docType methods
setMethod(f = "getWells", signature = "LFRset", def = function(object) {
    t = table(as.data.frame(object)$Well_ID)
    wells = as.vector(t)
    names(wells) = names(t)
    wells
})

#' Retrieves all the LFRs oiginating from a given well
#'
#' @name getWellLFRset
#' @rdname getWellLFRset-methods
#' @param x an object of class \code{\link{LFRset}}
#' @param well either an object of class \code{\link{numeric}} representing the well number or an object of class \code{\link{character}} repreenting the well nanme.
#' @return An object of class \code{\link{LFRset}} containing the Long Fragments Reads in x located in the well \code{well}
#' @seealso \code{\link{getLFRset}}, \code{\link{getLFRStats}}, \code{\link{getLFRStatsPerWell}}, \code{\link{summary.LFRset}}
#' @examples
#' #Loading an object of class (LFRset)
#' data(FiveWellsLFRset)
#' print(getWells(FiveWellsLFRset))
#' #A005  A016  L006  L024  P009
#' #43755 15192 23448 21892 24397
#' well1_LFR=getWellLFRset(FiveWellsLFRset,1)
#' print(getWells(well1_LFR))
#' #A005
#' #43755
#' P009_LFR=getWellLFRset(FiveWellsLFRset,'P009')
#' print(getWells(P009_LFR))
#' #P009
#' #24397
#' @docType methods
#' @export
setGeneric("getWellLFRset", def = function(object, well) {
    standardGeneric("getWellLFRset")
})

#' Retrieves all the LFRs oiginating from a given well
#'
#' @rdname getWellLFRset-methods
#' @aliases getWellLFRset,getWellLFRset-method
#' @docType methods
setMethod(f = "getWellLFRset", signature = "LFRset", def = function(object, well) {
    # if well is numeric then it represents the number of the well, if a string then it represent the name of the well
    fragments_df = object@Table
    if (as.character(class(well)) == "character") {
        subfragments_df = fragments_df[!is.na(fragments_df$Well_ID) & fragments_df$Well_ID == as.character(well), ]
    } else if (!is.na(as.numeric(well))) {
        subfragments_df = fragments_df[fragments_df$Well_Number == as.numeric(well), ]
    } else {
        stop("\n The paramete well should be either numeric (The well number) or a string (The well name)")
    }

    LFRset(Table = subfragments_df)
})



#' Retrieves the fraction of the whole genomes covered by the Long fragments Reads  making the LFRset
#'
#' @name getBreadthGenomeCoverage
#' @rdname getBreadthGenomeCoverage-methods
#' @param x an object of class \code{\link{LFRset}}
#' @return If a single well, a numeric value representing the breadth coverage of its LFRs. If multiple wells, the methods output a \code{\link{data.frame}} containing the breadth coverage per wells.
#' @seealso \code{\link{getN50}}, \code{\link{getLFRStats}}, \code{\link{getLFRStatsPerWell}}, \code{\link{summary.LFRset}}
#' @examples
#' #Loading an object of class (LFRset)
#' data(FiveWellsLFRset)
#' P009_LFR=getWellLFRset(FiveWellsLFRset,'P009')
#' print(getBreadthGenomeCoverage(P009_LFR))
#' # 4.63
#' print(getBreadthGenomeCoverage(FiveWellsLFRset))
#' # Wells BreadthCoverage
#' # 1  A005            7.22
#' # 2  A016            1.91
#' # 3  L006            3.69
#' # 4  L024            2.76
#' # 5  P009            4.63
#' @docType methods
#' @export
setGeneric("getBreadthGenomeCoverage", def = function(object) {
    standardGeneric("getBreadthGenomeCoverage")
})

#' Compute the percentage of the genome covered by the LFR set.
#'
#' @rdname getBreadthGenomeCoverage
#' @aliases getBreadthGenomeCoverage,getBreadthGenomeCoverage-method
#' @docType methods
setMethod(f = "getBreadthGenomeCoverage", signature = "LFRset", def = function(object) {
    LFR_df = object@Table
    wells_list = unique(object@Table$Well_ID)
    wellCoverage <- c()
    for (well in wells_list) {
        subLFR_df = LFR_df[LFR_df$Well_ID == well, ]

        totalgenomecovered = sum(as.numeric(subLFR_df[subLFR_df$Chrom %in% chrom_list, "Length"]), na.rm = T)
        genome_length = sum(unlist(hg19_chrsize[1:23]))
        ratio = totalgenomecovered/genome_length
        percentagecverage = as.numeric(format(100 * ratio, digit = 3))
        wellCoverage <- c(wellCoverage, percentagecverage)
    }
    if (length(wellCoverage) == 1) {
        wellCoverage
    } else {
        coverage_df = as.data.frame(cbind(wells_list, as.numeric(wellCoverage)))
        colnames(coverage_df) = c("Wells", "BreadthCoverage")
        # coverage_df$BreadthCoverage=as.numeric(coverage_df$BreadthCoverage)
        coverage_df
    }

})


#' Compute the N50 value of an LFRset
#'
#' @name getN50
#' @rdname getN50-methods
#' @param x an object of class \code{\link{LFRset}}
#' @return A numeric containing the N50 value  of the LFR set.
#' @seealso \code{\link{getBreadthGenomeCoverage}}, \code{\link{getLFRStats}}, \code{\link{getLFRStatsPerWell}}, \code{\link{LFRset}}, \code{\link{getWellLFRset}},
#' @examples
#' data(FiveWellsLFRset)
#' print(getN50(FiveWellsLFRset))
#' @docType methods
#' @export
setGeneric("getN50", def = function(object) {
    standardGeneric("getN50")
})

#' Compute the N50 of an LFR set.
#'
#' @rdname getN50
#' @aliases getN50,getN50-method
#' @docType methods
setMethod(f = "getN50", signature = "LFRset", definition = function(object) {
    frag_length = object@Table$Length  # Retrieve all the lengths
    frag_length = rev(sort(frag_length))
    n50 = frag_length[cumsum(frag_length) >= sum(frag_length)/2][1]
    n50
})



#' A function returning the following statistics from an LFR dataset
#'
#' \describe{
#'   \item{N}{The number of Long Fragments reads in the set}
#'   \item{N50_value}{The N50 of the LFR set}
#'   \item{N50_totallength}{The Sum of LFR lengths greater or equal to the N50}
#'   \item{L50}{The L50 of the set of LFR}
#'   \item{L50_percent}{Percentage of LFR having a length greater or equal to the N50}
#'   \item{Max}{Maximum LFR length in the set}
#'   \item{Min}{Minimum LFR length in the set }
#'   \item{Mean}{TThe mean of the LFR lengths}
#'   \item{Median}{The median of the LFR length}
#' }
#'
#' @name getLFRStats
#' @rdname getLFRStats-methods
#' @param x an object of class \code{\link{LFRset}}
#' @return A list containing the following statistic on the LFR set :
#' N, N50_value, N50_totallength, L50_value, L50_percent, Max,Min, Mean, Median .
#' @seealso \code{\link{getN50}}, \code{\link{getLFRStatsPerWell}}, \code{\link{getBreadthGenomeCoverage}}, \code{\link{plotLFRDistribution}}
#' @examples
#' data(FiveWellsLFRset)
#' P009_LFR=getWellLFRset(FiveWellsLFRset,'P009')
#' print(getLFRStats(P009_LFR))
#' print(getLFRStats(FiveWellsLFRset))
#' @docType methods
#' @export
setGeneric("getLFRStats", def = function(object) {
    standardGeneric("getLFRStats")
})

#'  A function returning the following statistics from an LFR dataset.
#'
#' @rdname getLFRStats
#' @aliases getLFRStats,getLFRStats-method
#' @docType methods
setMethod(f = "getLFRStats", signature = "LFRset", definition = function(object) {
    frag_length = object@Table$Length  # Retrieve all the lengths
    frag_length = rev(sort(frag_length))
    frag_cumlength = cumsum(frag_length)
    sumfrag_length = sum(frag_length)
    n50value = frag_length[frag_cumlength >= sumfrag_length/2][1]
    n50totallength = frag_cumlength[frag_cumlength >= sumfrag_length/2][1]
    L50value = which(frag_length == n50value)[1]
    L50percent = L50value * 100/length(frag_length)[1]
    stats = as.numeric(format(c(length(frag_length), n50value, n50totallength, L50value, L50percent, max(frag_length), min(frag_length), mean(frag_length), median(frag_length)), digits = 2,scientific=F))

    names(stats) = c("N", "N50_value", "N50_totallength", "L50_value", "L50_percent", "Max", "Min", "Mean", "Median")
    stats
})



#' A function returning the following statistics per wells from an LFR dataset
#'
#' \describe{
#'   \item{N}{The number of Long Fragments reads in the set}
#'   \item{N50_value}{The N50 of the LFR set}
#'   \item{N50_totallength}{The Sum of LFR lengths greater or equal to the N50}
#'   \item{L50}{The L50 of the set of LFR}
#'   \item{L50_percent}{Percentage of LFR having a length greater or equal to the N50}
#'   \item{Max}{Maximum LFR length in the set}
#'   \item{Min}{Minimum LFR length in the set }
#'   \item{Mean}{TThe mean of the LFR lengths}
#'   \item{Median}{The median of the LFR length}
#' }
#'
#' @name getLFRStatsPerWell
#' @rdname getLFRStatsPerWell-methods
#' @param x an object of class \code{\link{LFRset}}
#' @return A data frame  containing for each well having an LFR in the LFR set  the following statistics :N, N50_value,N50_totallength,L50_value,L50_percent,Max,Min,Mean,Median .
#' @seealso \code{\link{getN50}}, \code{\link{getLFRStats}}, \code{\link{getBreadthGenomeCoverage}}, \code{\link{plotLFRDistribution}}, \code{\link{summary.LFRset}}, \code{\link{getLFRset}}, \code{\link{plot.LFRset}}
#' @examples
#' data(FiveWellsLFRset)
#' P009_LFR=getWellLFRset(FiveWellsLFRset,'P009')
#' print(getLFRStatsPerWell(FiveWellsLFRset))
#' @docType methods
#' @export
setGeneric("getLFRStatsPerWell", def = function(object) {
    standardGeneric("getLFRStatsPerWell")
})

#' A function returning the following statistics per wells from an LFR dataset
#'
#' @rdname getLFRStatsPerWell
#' @aliases getLFRStatsPerWell,getLFRStatsPerWell-method
#' @docType methods
setMethod(f = "getLFRStatsPerWell", signature = "LFRset", definition = function(object) {
    wells_list = unique(object@Table$Well_ID)
    stats_df <- data.frame()
    for (iwell in wells_list) {
        df = object@Table
        frag_length = df[df$Well_ID == iwell, ]
        stats_df = rbind(stats_df, getLFRStats(LFRset(Table = frag_length)))
    }
    stats_df = cbind.data.frame(wells_list, stats_df)
    names(stats_df) = c("Well", "N", "N50_value", "N50_totallength", "L50_value", "L50_percent", "Max", "Min", "Mean", "Median")
    stats_df
})


#' Plot the LFR distribution
#'
#' Plot the distribution of the LFR across the genome and across all the wells. Suitable to be call on genomic region of less than 2MB length. Segments represent fragments, the y-axis is the wells and the x-axis the genomic location. The thickness of the segment is proportional to the averagecoverage. The blue segment is an illustration of a fagment of 100kbp (if Region <5Mb) or 1Mb if (region >5Mb)  with and average coverage of 1.
#'
#' @name plotLFRDistribution
#' @rdname plotLFRDistribution-methods
#' @param x an object of class \code{\link{LFRset}}
#' @param chrom the chromosome to be considered
#' @param region The genomic region to plot. For best results consider a region of less than 2MB.
#' @param minLength Conside ronly the LFR having a length greater or equal to minLength.
#' @return ''
#' @seealso \code{\link{plot}}, code{\link{summary.LFRset}}, \code{\link{getLFRset}}, \code{\link{getWells}}
#' @examples
#' cat('\n Loading an LFRset data frame ')
#' LFR_11152Tissue_10k_file<-system.file('extdata',
#'        'D1000_10k_allfragments.tsv.gz',package = 'DigitalPicoTools')
#' D1000_10k_11152Tissue_LFR=read.table(file=LFR_11152Tissue_10k_file)
#' D1000_10k_11152Tissue_LFRset=LFRset(Table=D1000_10k_11152Tissue_LFR)
#' cat('\n Plotting the Fragment Distribution ')
#' plotLFRDistribution(D1000_10k_11152Tissue_LFRset,'chr2',region='chr22:22000000-23000000')
#' @docType methods
#' @export
setGeneric("plotLFRDistribution", def = function(object, chrom, region = NULL, minLength = 10000) {
    standardGeneric("plotLFRDistribution")
})


#' A method for plotting the LFR distribution
#' @rdname plotLFRDistribution
#' @aliases plotLFRDistribution,plotLFRDistribution-method
#' @docType methods
setMethod(f = "plotLFRDistribution", signature = "LFRset", definition = function(object, region = NULL, minLength = 10000) {
    plot_LFRChromosomeDistribution(object@Table, region, minLength)
})


#' Summary statistics on an LFRset Object
#'
#' This is a generic function to compute the common summary statistic on the Length,
#' Average coverage and the numbe rof reads of an LFR set object
#'
#' @param x an  \code{\link{LFRset}} object containing the list of Long Fragment Reads
#' @return The five statistics summary (Min, Max, Mean, Median, 1st Quart, 3rd Quart) for the Length, the Average Coverage and the Number of Reads of the LFRs.
#' @seealso \code{\link{getBreadthGenomeCoverage}}, \code{\link{getLFRStats}}, \code{\link{getLFRset}},  \code{\link{getLFRStatsPerWell}}, \code{\link{getN50}}
#' @examples
#' LFR_11152Tissue_10k_file<-system.file('extdata',
#'       'D1000_10k_allfragments.tsv.gz',package = 'DigitalPicoTools')
#' D1000_10k_11152Tissue_LFR=read.table(file=LFR_11152Tissue_10k_file)
#' D1000_10k_11152Tissue_LFRset=LFRset(Table=D1000_10k_11152Tissue_LFR)
#' summary(D1000_10k_11152Tissue_LFRset)
#' #  Summary on the LFRset object
#' #
#' # Fragments Length
#' # Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#' # 10000   12340   15800   18700   21920  150200
#' #
#' # Fragments Average Coverage
#' # Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#' # 0.2747  0.9621  1.3730  1.7560  2.1050 28.9100
#' #
#' # Fragments Nbreads
#' # Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#' # 12.0    53.0    90.0   154.4   174.0  5240.0
#'
#' @export
summary.LFRset <- function(object) {
    cat(" \n\t Summary on the LFRset object ")

    cat("\n\n\t Fragments Length\n")
    print(summary(object@Table$Length))

    cat("\n\n\t Fragments Average Coverage\n")
    print(summary(object@Table$AvgCoverage))

    cat("\n\n\t Fragments Nbreads\n")
    print(summary(object@Table$Nbreads))



}


#' Plotting an LFR data set
#'
#' This is a generic function to plot  in different way the Length and the Average coverage of the long Frangment Reads within an LFR data set. Two type of plot can be performed.
#' \describe{
#'   \item{ScatterPlot}{this plot represent a scatter plot between the Length and the Average Coverage}
#'   \item{Landscape}{This plot type represent the distribution of the LFR across the whole genome. It adds a third dimension to the length and the coverage which is the genomic location. }
#' }
#'
#'  In both types, the paramater \emph{Value} allow to specifie which variable will be represented on the y-axis. It can be  'Coverage' for the LFR Average Coverage and 'Length' for the LFR lengths.
#'
#' @param x an  \code{\link{LFRset}} object containing the list of Long Fragment Reads
#' @param type the type of the plot. Can be either 'ScatterPlot' either 'Landscape'. The default value is 'ScatterPlot'.
#' @param Value The variable to represent on the y-Axis. Can be either 'Coverage' either 'length'. default is 'Coverage'.
#' @param minlength Minimum length of a fragment to be considered. Default 10000.
#' @param color The points color. default: black
#' @param main The title of the plot. default: None
#' @return ''
#' @seealso \code{\link{plotLFRDistribution}}, \code{\link{plot.LFRset}}, \code{\link{getWells}}, \code{\link{getWellLFRset}}, \code{\link{getLFRStats}}
#' @examples
#'
#'  #Example 1:
#'  wellsID_file=system.file('extdata','wells_id_tissue.txt',package='DigitalPicoTools')
#'  wellsID_list<-unlist(read.table(wellsID_file))
#'  WellSample_list=initsWellSamples(wellsID_list,
#'      BAMprefix=paste(path.package('DigitalPicoTools'),'extdata/',sep='/'))
#'  LFR_Info = getLFRset(WellSample_list,3000)
#'  well1=getWellLFRset(LFR_Info,1)
#'  plot(well1,value='Length')
#'  plot(LFR_Info,type='Landscape', color='magenta',main='test')
#'  plot(well1,value='Coverage',type='Landscape')
#' @export
plot.LFRset <- function(object, type = "ScatterPlot", value = "Coverage", minlength = 10000, color = "black", main = "", well = NULL) {

    plot_fragments(object@Table, type =type, value = value,  minlength = minlength, color = color, main = main, well = well)
}




#' Plotting histograms on an LFR data set
#'
#' This is a generic function to plot  the histogram of  the Length or the Average coverage of the long Frangment Reads within an LFR data set. given the paramter \emph{value}, two values can be plotted: The fragment length and the fragment average coverage
#'
#' @param x an  \code{\link{LFRset}} object containing the list of Long Fragment Reads
#' @param Value The variable to plot the histogram for . Can be either 'Coverage' either 'length'. default is 'Coverage'.
#' @param minlength Minimum length of a fragment to be considered. Default 10000.
#' @param color histogram color. default: black
#' @param main The title of the plot. default: None
#' @return ''
#' @seealso \code{\link{plotLFRDistribution}}, \code{\link{plot.LFRset}}, \code{\link{getWells}}, \code{\link{getWellLFRset}}, \code{\link{getLFRStats}}
#' @examples
#'
#'  #Example 1:
#'  wellsID_file=system.file('extdata','wells_id_tissue.txt',package='DigitalPicoTools')
#'  wellsID_list<-unlist(read.table(wellsID_file))
#'  WellSample_list=initsWellSamples(wellsID_list,
#'      BAMprefix=paste(path.package('DigitalPicoTools'),'extdata/',sep='/'))
#'  LFR_Info = getLFRset(WellSample_list,3000)
#'  well1=getWellLFRset(LFR_Info,1)
#'  hist(well1,value='Length')
#'  hist(LFR_Info,value='Coverage', color='magenta',main='test')
#' @export
hist.LFRset <- function(object, value = "Coverage", minlength = 10000, color = "black", main = "", well = NULL) {

  hist_fragments(object@Table,  value = value, minlength = minlength, color = color, main = main, well =well)
}







#' The SampleCoverage class
#'
#' An S4 class containing the genomic location and allele information per sample for each variant from a VCF file.
#'
#' @slot Table  A data frame containing the genomic location and allele information per sample for each variant. It gives for each called variants the following information :
#' \describe{
#'   \item{Chrom}{The chromosome of the variant)}
#'   \item{Pos}{The position of the variant}
#'   \item{REF}{The reference sequence }
#'   \item{ALT}{The altered sequence}
#'   \item{Qual}{The quality of the variant call}
#'   \item{Filter}{The VCF filter of the variant call}
#'   \item{TC}{The total reads count supporting the variant}
#'   \item{TR}{The total reads count supporting the variant}
#'   \item{AF}{The reads allele fraction}
#'   \item{}{For each  well samples}
#'   \describe{
#'      \item{NR_sampleID}{The total numbe rof reads covering the variant location at this smaple}
#'      \item{NV_sampleID}{Number of reads supporting the variant sequence at this smaple}
#'      \item{AF_sampleID}{Variant Reads allele fraction at this sample}
#'   }
#'   \item{GT}{colon separated list of unphased genotype at each sample}
#' }
#'
#' @seealso \code{\link{VariantAlleleInfo}}, \code{\link{getVariantAlleleInfo}}, \code{\link{getVariantCoverageTable}}, \code{\link{plot.VariantAlleleInfo}}
#' @name SampleCoverage
#' @rdname SampleCoverage
#' @aliases SampleCoverage-class
#' @export SampleCoverage
#' @exportClass SampleCoverage
SampleCoverage <- setClass("SampleCoverage", slots = character(0), prototype = list(Table = data.frame()), validity = function(object) {


}, contains = "LFRdataframe")





#' The VariantAlleleInfo class
#'
#' An S4 class containing the allele information for each variant from a VCF file.
#'
#' @slot Table  A data frame containing the Allele information per variant. It gives for each called variants the following information :
#' \describe{
#'   \item{Chrom}{The chromosome of the variant)}
#'   \item{Pos}{The position of the variant}
#'   \item{REF}{The reference sequence }
#'   \item{ALT}{The altered sequence}
#'   \item{Qual}{The quality of the variant call}
#'   \item{Filter}{The VCF filter of the variant call}
#'   \item{TC}{The total reads count supporting the variant}
#'   \item{TR}{The total reads count supporting the variant}
#'   \item{AF}{The reads allele fraction}
#'   \item{WR}{The total well count at the position of the variant}
#'   \item{WV}{The count of well supporting the variant}
#'   \item{Well_Number}{The wells allele fraction}
#' }
#'
#' @seealso \code{\link{SampleCoverage}}, \code{\link{VariantAlleleInfo}}, \code{\link{getVariantAlleleInfo}},
#'          \code{\link{getVariantCoverageTable}}, \code{\link{plot.VariantAlleleInfo}}
#' @note This class is one of the main class of the package
#' @name VariantAlleleInfo
#' @rdname VariantAlleleInfo
#' @aliases VariantAlleleInfo-class
#' @export VariantAlleleInfo
#' @exportClass VariantAlleleInfo
VariantAlleleInfo <- setClass("VariantAlleleInfo", slots = character(0), prototype = list(Table = data.frame()), validity = function(object) {


}, contains = "LFRdataframe")


#' Plotting the variant allele information
#'
#' This is a generic function to plot different information on the variant allele counts. The information to plot is set by the parameter \emph{Value} which can be eithe rof the following :
#' \describe{
#'   \item{WellsFraction}{Wells allele fraction of the variant}
#'   \item{WellsCounts}{Wells allele counts of the variants }
#'   \item{ReadsFraction}{Reads allele fraction of the variant}
#'   \item{ReadsCounts}{Reads allele counts of the variants }
#' }
#'
#'  In both types, the paramater \emph{Value} allow to specifie which variable will be represented on the y-axis. It can be  'Coverage' for the LFR Average Coverage and 'Length' for the LFR lengths.
#'
#' @param x an  \code{\link{VariantAlleleInfo}} object containing the list of variant and their genomic location and allele count information
#' @param region Region to be plotted.  If NULL then the whole region covered by the variants will be considered
#' @param samplingRatio Ratio of points to plot. Only 1 each samplingRatio points will be plotted.
#' @param Value Information to plot. Can be one of WellsFraction, WellsCounts, ReadsFraction, ReadsCounts.
#' @param main The title of the plot. default: None
#' @return ''
#' @seealso \code{\link{getVariantAlleleInfo}}, \code{\link{getVariantCoverageTable}}, \code{\link{VariantAlleleInfo}}
#' @examples
#'
#' #Example 1:
#' vcffile<-system.file('extdata','LFR_11152_Tissue_chr22_1000variants.vcf.gz',package='DigitalPicoTools')
#' AlleleInfo_df = getVariantAlleleInfo(vcffile)
#' plot(AlleleInfo_df,'chr22')
#'
#' @export
plot.VariantAlleleInfo <- function(object, region = NULL, samplingRatio = 1, value = "WellsFraction", main = "") {


    plot_AlleleInfo(as.data.frame(object), region, samplingRatio, value, main)


}

VariantAlleleInfo.summary <- function() {

}


#' Retrieveing variants genomic and alleles count information
#'
#'This function takes as input a VCF file obtained from an LFR experiment and output and returns of class \code{\link{VariantAlleleInfo}} containing the genomic and allele count information of each variants.
#' @note  The function can also be used for a multiple sample and non LFR experiment. In this case, the well allele count and fraction  information will assess the presence of the variant across the different samples.
#' @param VCFFilePath A character string representing the path to a VCF file obtained from an LFR experiment study. The file should contains one genotype information column per well sample.
#' @param wells_id   A vector of character string representing the well sample ID to assign to the columns name corresponding to each sample genotyping information. If not provided the column names in the VCF file are maintained.  length(well_id) should be equal to the number of column of the VCF file - 9 following to the VCF format specifications.
#' @param region A character string representing the region of the genome to consider given in the following format :  'chrom:start-end'. Default value :NULL and all the variants in the VCF file will be considered.
#' @param OnlyPassFilter A logical value. If set to true, only the variant having their calling filter status set at 'PASS'
#' @param FiltersToInclude A vector of character string representing a list of calling filter status. All the variant having any of them will be included in the analysis.
#' @param FiltersToExclude A vector of character string representing a list of calling filter status. All the variant having any of them will be excluded from  the analysis given that the filte ris not included in FiltersToInclude.
#' @param RetrieveWellsIDs A logical value specifying if the wells ID supporting the reference and  the SNP should be retrieved or not. Usefull for the Phasing analysis. (Default : FALSE)
#' @return An object of class  \code{\link{VariantAlleleInfo}}.
#' @examples
#'  vcffile<-system.file('extdata','LFR_11152_Tissue_chr22.vcf.gz',package='DigitalPicoTools')
#'  AlleleInfo = getVariantAlleleInfo(vcffile)
#'  print(head(AlleleInfo))
#'  #                     Chrom      Pos REF ALT Qual          Filter TC TR        AF WR WV        WF     WFAdj
#'  # chr22_16060515_CTT_C chr22 16060515 CTT   C   68  SC;badReads;MQ  4  3 0.7500000  4  3 0.7500000 0.7500000
#'  # chr22_16114149_C_A   chr22 16114149   C   A   36  SC;badReads;MQ 10  2 0.2000000 10  2 0.2000000 0.2000000
#'  # chr22_16114636_CT_C  chr22 16114636  CT   C   28  badReads;MQ;QD 13  4 0.3076923 11  4 0.3636364 0.3636364
#'  # chr22_16114662_G_A   chr22 16114662   G   A  147     badReads;MQ 13  5 0.3846154  8  5 0.6250000 0.6250000
#'  # chr22_16194000_A_AT  chr22 16194000   A  AT   39     badReads;MQ  2  2 1.0000000  2  2 1.0000000 0.9558641
#'  # chr22_16228052_T_C   chr22 16228052   T   C   18 Q20;badReads;MQ  6  2 0.3333333  5  2 0.4000000 0.4000000
#' @seealso \code{\link{getVariantCoverageTable}}, \code{\link{plot.VariantAlleleInfo}}, \code{\link{SampleCoverage}}, \code{\link{VariantAlleleInfo}}
#' @export
getVariantAlleleInfo <- function(VCFFilePath, wells_id = NULL, region = NULL, part=NULL, OnlyPassFilter = FALSE, FiltersToInclude = NULL, FiltersToExclude = NULL,RetrieveIDs=FALSE,mode="wells") {
    # cat('\nLoading the vcf ...')
    variant_df = loading_vcf(VCFFilePath, wells_id, region)

    # cat('\n Retrieving variants with filter PASS ...')
    if (OnlyPassFilter) {
      variant_df = select_passed(variant_df)
    } else {
      variant_df = variant_df

        if (!is.null(FiltersToInclude) || !is.null(FiltersToExclude))
          variant_df = select_filter(variant_df, FiltersToInclude, FiltersToExclude)
    }


    start_position=1
    end_position=nrow(variant_df)

    cat("/n number of rows is : ", nrow(variant_df))

    if (!is.null(part)){
      number_of_parts=as.numeric(unlist(strsplit(part,":"))[2])
      part_index=as.numeric(unlist(strsplit(part,":"))[1])

      part_size= ceiling(nrow(variant_df)%/% number_of_parts)
      start_position = (part_index-1) * part_size + 1
      end_position = (part_index-1) * part_size + part_size
      if(end_position>nrow(variant_df)){
        end_position=nrow(variant_df)
      }

    }


    variant_df = variant_df[start_position:end_position,]
    cat("\n Extracting table of allele count ...")
    cat("\n \t ",end_position-start_position +1 ," mutations to process : From mutation :", start_position, " to mutation ", end_position," ...")
    # variantWGA_df_tabular=get_tabular(variantWGA_df_passed,wells_id)

    variant_df = get_alleleinfotabular(variant_df,RetrieveIDs=RetrieveIDs, mode=mode)

    VariantAlleleInfo(Table = variant_df)

}


#' Retrieveing variants genomic and alleles count information per sample
#'
#'This function takes as input a VCF file obtained from an LFR experiment and output and returns of class \code{\link{SampleCoverage}} containing the genomic and allele count information per sample of each variants.
#' @note  The function can also be used on a multiple sample and non LFR experiment.
#' @inheritParams getVariantAlleleInfo
#' @return An object of class  \code{\link{SampleCoverage}}.
#' @examples
#'  vcffile<-system.file('extdata','LFR_11152_Tissue_chr22.vcf.gz',package='DigitalPicoTools')
#'  CoverageTable = getVariantCoverageTable(vcffile)
#'  print(head(CoverageTable))
#' @seealso \code{\link{getVariantAlleleInfo}}, \code{\link{plot.VariantAlleleInfo}}, \code{\link{SampleCoverage}}, \code{\link{VariantAlleleInfo}}
#' @export
getVariantCoverageTable <- function(VCFFilePath, wells_id = NULL, region = NULL, OnlyPassFilter = FALSE, FiltersToInclude = NULL, FiltersToExclude = NULL) {


    cat("\nLoading the vcf ...")
    variant_df = loading_vcf(VCFFilePath, wells_id)

    cat("\n Retrieving variants with filter PASS ...")
    if (OnlyPassFilter) {
        variant_df = select_passed(variant_df)
    } else {
        variant_df = variant_df

        if (!is.null(FiltersToInclude) || !is.null(FiltersToExclude))
            variant_df = select_filter(variant_df, FiltersToInclude, FiltersToExclude)
    }

    cat("\n Extracting table of allele count ...")
    # variantWGA_df_tabular=get_tabular(variantWGA_df_passed,wells_id)
    variant_df_tabular = get_coveragetabular(variant_df)

    SampleCoverage(Table = variant_df_tabular)

}


.unlist <- function(x) {
    ## do.call(c, ...) coerces factor to integer, which is undesired
    x1 <- x[[1L]]
    if (is.factor(x1)) {
        structure(unlist(x), class = "factor", levels = levels(x1))
    } else {
        do.call(c, x)
    }
}



#' Computing Long fragment Reads from the wells
#'
#'This function takes as input a list of wells and compute the list of Long Fragments Reads sequenced in each well. Each LFR is represented by its name, its chromosome, its start and end position, its length, its average coverage, the number of reads making up the fragment and the well ID it is originating from. The function return an objoct of class \code{\link{LFRSet}}.
#' @note  'The function creates temproary index of the bam files. So you need to have writing accesss to the location of the Bam files'.
#' @param well_list A list of objects of class \code{\link{WellSample}} representing the wells to compute the fragments from. well_list can also be a vector of character string representing the file names to each well BAM file.
#' @param mindistance The minimum distance separating two consecutive LFR. Two paired reads located within that distance are considered to be originating from the same fragments of DNA.
#' @param unexists.action Define the action to perform when a well bam file do not exists. Can be one of the following:
#' \describe{
#'   \item{\emph{unexists.pass}}{No particular action, consider the provided bam file and continue.}
#'   \item{\emph{unexists.fail}}{Interupt the run and generate an error}
#'   \item{\emph{unexists.exclude}}{Do not process the sample, skipt it or remove it from the list of well and continue}
#'  }
#'  The default value is set at 'unexists.fail'
#' @param minReads Minimum number of reads a fragment should contains to be considered as valid. Default value : 2
#' @param minFragmentsLength  Minimum length of a fragment for it to be considered as valid. default : 1000
#' @param region  A character string representing the region of the genome to consider given in the following format :  'chrom:start-end'. Default value :NULL and the Long Fragments will be retrieved across the whole BAM file.
#' @return An object of class  \code{\link{LFRset}}.
#' @examples
#'   wellsID_file=system.file('extdata','wells_id_tissue.txt',package='DigitalPicoTools')
#'   wellsID_list<-unlist(read.table(wellsID_file))
#'   WellSample_list=initsWellSamples(wellsID_list,
#'    BAMprefix=paste(path.package('DigitalPicoTools'),  'extdata/',sep='/'))
#'   FiveWellsLFRset = getLFRset(WellSample_list,3000)
#'   print(head(FiveWellsLFRset))
#'
#'   # An object of class LFRdataframe with   rows and   columns
#'   #                   LFR_name Chrom  Start    End Length AvgCoverage Nbreads Well_Number Well_ID
#'   # 1   A005_chr1_47518_51209  chr1  47518  51209   3691   0.2915199       4           1    A005
#'   # 2   A005_chr1_64725_72351  chr1  64725  72351   7626   0.2265932       7           1    A005
#'   # 3   A005_chr1_77103_78637  chr1  77103  78637   1534   0.6147327       4           1    A005
#'   # 4   A005_chr1_82636_87332  chr1  82636  87332   4696   0.3017462       6           1    A005
#'   # 5  A005_chr1_98437_104080  chr1  98437 104080   5643   0.3202197       7           1    A005
#'   # 6 A005_chr1_110942_113908  chr1 110942 113908   2966   0.2862441       3           1    A005

#' @seealso \code{\link{getWellLFRset}}, \code{\link{getLFRStats}, \code{\link{plotLFRDistribution}}}, \code{\link{getWells}},  \code{\link{plotLFRDistribution}}, \code{\link{plot.LFRset}}, \code{\link{summary.LFRset}}, \code{\link{getLFRStatsPerWell}}, \code{\link{getBreadthGenomeCoverage}},  \code{\link{getN50}}
#' @export
getLFRset <- function(wells_list, mindistance, unexists.action = "unexists.fail", minNbReads = 2, minFragmentsLength = 1000,  region = NULL)
  {
   chrom_list<-paste('chr',c(1:22,'X','Y','M'),sep='')
    LFR_wgs_df <- data.frame()
    for (iwell in 1:length(wells_list)) {
        if (class(wells_list[iwell]) == "character") {
            BAMfile_name <- wells_list[iwell]
            Well_Number = iwell
            filenameparts = unlist(strsplit(BAMfile_name, "/"))
            filenameparts_filename = filenameparts[length(filenameparts)]
            Well_ID = unlist(strsplit(filenameparts_filename, "\\."))[1]
            #We need to take the first string before any "_" since it is used as separators in the fragment name
            if(grepl("_",Well_ID  )){
              warning(paste(" The well ID ", Well_ID," contains the character _,  ", unlist(strsplit(Well_ID, "_"))[1]  , " will be considered"))
              Well_ID = unlist(strsplit(Well_ID, "_"))[1]
            }

            # WellID=paste('W', Well_Number,sep='')
        } else if (as.character(class(wells_list[[iwell]])) == "WellSample") {
            BAMfile_name <- wells_list[[iwell]]@BAMfilename
            Well_Number = wells_list[[iwell]]@Well_Number
            Well_ID = wells_list[[iwell]]@Well_ID
            # WellID=wells_list[[iwell]]@ID
        } else {
            stop(paste(" < ", as.character(class(wells_list[[iwell]])), ">  Incorrect class for well ", iwell, " Well class should be either WellSample (S4 object WellSample) either a string (the bam file name)",
                sep = ""))
        }



        if (!file.exists(BAMfile_name)) {
            if (unexists.action == "unexists.fail") {
                stop("\n The well bam file  ", BAMfile_name, " do not exists, Check the parameters or set unexists.action to unexists.exclude ")
            } else {
                warning("\n The well bam file  ", BAMfile_name, " do not exists, Well will be skipped")
            }
        }

        cat("\n\n Processing  ", Well_ID, "\n")


        BAMFile = BAMfile_name
        indexBam(BAMFile)
        p2 <- ScanBamParam(what = scanBamWhat())
        BAM <- scanBam(BAMFile, param = p2)
        # store names of BAM fields
        BAM_field <- names(BAM[[1]])
        # go through each BAM field and unlist
        list <- lapply(BAM_field, function(y) .unlist(lapply(BAM, "[[", y)))
        # store as data frame
        BAM_df <- do.call("DataFrame", list)
        names(BAM_df) <- BAM_field

        #Extract only the region concerned

        if (!is.null(region) && region!="" && region!="*") {
          GRegion = parseregion(region)
          if (!is.null(GRegion$Chrom)) {
            BAM_df = BAM_df[!is.na(BAM_df$rname) & as.character(BAM_df$rname) ==  GRegion$Chrom,]
            if (nrow(BAM_df)==0){
              warnings(" \n No sequencing reads in the Bam file corresponding to the region specified. Bam file will be skipped \n\n")
              next
            }
            if (!is.null(GRegion$Start) && !is.null(GRegion$End)) {
              startPosition = GRegion$Start
              endPosition = GRegion$End
              BAM_df = BAM_df[BAM_df$pos >= startPosition & BAM_df$pos <= endPosition, ]
              if (nrow(BAM_df)==0){
                warnings(" \n No sequencing reads in the Bam file corresponding to the region specified. Bam file will be skipped \n\n")
                next
              }
            }
          }
        }


        # Consider only the reads having both mate in the same chromosome
        BAM_df = BAM_df[!is.na(BAM_df$mrnm) & !is.na(BAM_df$rname) & BAM_df$mrnm == BAM_df$rname, ]

        # We need to extract some flags value
        isfirstmate <- function(flag) {
            as.integer(intToBits(as.integer(flag)))[7]
        }  # extract the fourth bits before the last. In flag bit 4 is other segment unmapped.
        issecondmate <- function(flag) {
            as.integer(intToBits(as.integer(flag)))[8]
        }  # extract the fourth bits before the last. In flag bit 4 is other segment unmapped.
        isconcordant <- function(flag) {
            as.integer(intToBits(as.integer(flag)))[2]
        }  # extract the fourth bits before the last. In flag bit 4 is other segment unmapped.
        BAM_df["IsFirstMate"] = unlist(lapply(unlist(BAM_df$flag), function(x) isfirstmate(as.numeric(x))))
        BAM_df["IsSecondMate"] = unlist(lapply(unlist(BAM_df$flag), function(x) issecondmate(as.numeric(x))))
        BAM_df["IsConcordant"] = unlist(lapply(unlist(BAM_df$flag), function(x) isconcordant(as.numeric(x))))



        ### }



        chrom_list=intersect(chrom_list,as.character(unique(BAM_df$rname)))

        for (mychrom in chrom_list) {


            cat(" ", mychrom, "....")
            BAM_df_chr = BAM_df[BAM_df$rname == mychrom, ]

            if(nrow(BAM_df_chr)==0)
              next

            BAM_mate1_df = BAM_df_chr[BAM_df_chr$IsFirstMate == 1, ]
            BAM_mate2_df = BAM_df_chr[BAM_df_chr$IsSecondMate == 1, ]
            rownames(BAM_mate1_df) = BAM_mate1_df$qname
            rownames(BAM_mate2_df) = BAM_mate2_df$qname

            # Retrieve proper paired reads
            pair_names = intersect(rownames(BAM_mate1_df[BAM_mate1_df$IsConcordant == 1, ]), rownames(BAM_mate2_df[BAM_mate2_df$IsConcordant ==
                1, ]))
            pair_mate1_df = BAM_mate1_df[pair_names, ]
            pair_mate2_df = BAM_mate2_df[pair_names, ]

            mate1_df = pair_mate1_df
            mate2_df = pair_mate2_df
            mate1_df["epos"] = as.numeric(as.character(unlist(mate1_df["pos"]))) + as.numeric(as.character(unlist(mate1_df["qwidth"])))
            mate2_df["epos"] = as.numeric(as.character(unlist(mate2_df["pos"]))) + as.numeric(as.character(unlist(mate2_df["qwidth"])))

            reads_df = mate1_df[c("rname", "pos", "qwidth", "mpos", "isize", "epos")]
            reads_df[c("qwidth2", "epos2")] = mate2_df[c("qwidth", "epos")]
            reads_df["insert"] = abs(as.numeric(as.character(unlist(reads_df["isize"]))))
            reads_df["rwidth"] = as.numeric(as.character(unlist(reads_df["qwidth"]))) + as.numeric(as.character(unlist(reads_df["qwidth2"])))
            reads_df["rstart"] = pmin.int(as.numeric(as.character(unlist(reads_df["pos"]))), as.numeric(as.character(unlist(reads_df["mpos"]))),
                na.rm = T)
            reads_df["rend"] = pmax.int(as.numeric(as.character(unlist(reads_df["epos"]))), as.numeric(as.character(unlist(reads_df["epos2"]))),
                na.rm = T)
            reads_df["length"] = pmin.int(as.numeric(as.character(unlist(reads_df["insert"]))), as.numeric(as.character(unlist(reads_df["rwidth"]))),
                na.rm = T)

            reads_df = reads_df[order(as.numeric(as.character(unlist(reads_df$rstart)))), ]


            reads_fragment = reads_df[c("rstart", "rend", "length")]
            # reads_fragment['gap']=c(0,diff(as.numeric(as.character(unlist(reads_fragment['rend'])))))
            N = nrow(reads_fragment)

            N = nrow(reads_fragment)
            if(N==0) next
            if(N==1){
              reads_fragment["gap"] = 1e+05
            }else{
              reads_fragment["gap"] = c(1e+05, as.numeric(as.character(unlist(reads_fragment[2:N, "rstart"]))) - as.numeric(as.character(unlist(reads_fragment[1:(N -
                                                                                                                                                                    1), "rend"]))))
            }

            reads_fragment["newfragment"] = as.numeric(reads_fragment$gap > mindistance)


            # plot(reads_fragment$gap, type='h')



            # Construct the fragments Indices of fragments
            fragments_indices <- which(reads_fragment$newfragment == 1)
            # Number of fragments
            Nfrag = length(fragments_indices)
            # Starts of fragments
            fragments_starts <- reads_fragment[fragments_indices, "rstart"]
            # Ends of fragments
            fragments_ends <- c(reads_fragment[fragments_indices - 1, "rend"], reads_fragment[N, "rend"])
            # Number fragments
            fragments_nbreads <- c(diff(fragments_indices), N - fragments_indices[Nfrag] + 1)
            # total bases in fragments
            fragments_sumbases = c()
            for (i in 1:length(fragments_indices)) fragments_sumbases = c(fragments_sumbases, sum(reads_fragment$length[fragments_indices[i]:(fragments_indices[i] +
                fragments_nbreads[i] - 1)]))

            # We constitute the fragments.
            LFR_df = matrix(ncol = 9, nrow = Nfrag)
            LFR_df <- as.data.frame(LFR_df)
            names(LFR_df) <- c("LFR_name", "Chrom", "Start", "End", "Length", "AvgCoverage", "Nbreads", "Well_Number", "Well_ID")
            LFR_df["LFR_name"] = paste(Well_ID, mychrom, fragments_starts, fragments_ends, sep = "_")
            LFR_df["Chrom"] = rep(mychrom, Nfrag)
            LFR_df["Start"] = fragments_starts
            LFR_df["End"] = fragments_ends
            LFR_df["Length"] = fragments_ends - fragments_starts
            LFR_df["AvgCoverage"] = fragments_sumbases/LFR_df$Length
            LFR_df["Nbreads"] = fragments_nbreads
            LFR_df["Well_ID"] = rep(Well_ID, Nfrag)
            LFR_df["Well_Number"] = rep(Well_Number, Nfrag)
            # LFR_df['WellID'] = rep(WellID,Nfrag)


            LFR_df2 = LFR_df[order(LFR_df$Length, decreasing = T), ]

            # cat('\n longuest fragments of the chromosome ') print(head(LFR_df2))

            # Select only fragment having not less than minNbReads reads
            LFR_df = LFR_df[as.numeric(LFR_df$Nbreads) >= minNbReads, ]

            # Select only fragment with a length not less than minFragmentsLength
            LFR_df = LFR_df[as.numeric(LFR_df$Length) >= minFragmentsLength, ]

            if (nrow(LFR_wgs_df) == 0) {
                LFR_wgs_df = LFR_df
            } else {
                LFR_wgs_df = rbind(LFR_wgs_df, LFR_df)
            }




        }

    }


    LFRset(Table = LFR_wgs_df)



}




#' Compute the Phaisng Code
#'
#'This function computes the phasing codes of each mutations. Mutations with the same phasing code are considered as phased (e.g PH_101_0 and PH_101_0) . Mutations with the same phasing fragment code but different allele index are considered as not phased(e.g PH_101_1 and PH_101_0) .
#'
#' @export
ComputePhasingCode_recursive<-function(LFR_withMutations_df,Mutations_set)
{








 cat("\n *******\n Phasing Code computation \n ********")
 cat("\n\t ", nrow(LFR_withMutations_df)," Long Fragment Reads to scan")
 cat("\n\t ", nrow(Mutations_set)," Mutations to phase")



 cat("\n\n Remove mutations with no wells count supporting the variant...")
 #This is odd, there are some mutations with no well count supporting the variant? We remove them from now.
 Mutations_set =Mutations_set[Mutations_set$WV!=0, ]
 cat("\n ", nrow(Mutations_set)," Mutations left..")

 cat("\n retrieve any potential homozygous mutations with zero Ref Allele count)")
 cat("\n\t ", nrow(Mutations_set[Mutations_set$WR==0, ]), " Potential homozygous mutations ")
 Mutations_set =Mutations_set[Mutations_set$WR!=0, ]


  #Matrix of Phasing_Code

  MutationPhasingCode_df<<-as.data.frame(matrix(nrow=nrow(Mutations_set),ncol=5))
  MutationPhasingCode_df[1:2]<<-Mutations_set[1:2]
  names(MutationPhasingCode_df) <<- c("Chrom","Pos",paste("PhasingCode",1:3,sep=""))
  rownames(MutationPhasingCode_df) <<- rownames(Mutations_set)

  list_of_mutations_to_phase=rownames(MutationPhasingCode_df)

  LFR_withMutations_df= LFR_withMutations_df[LFR_withMutations_df$MutationsOnSNP!="" | LFR_withMutations_df$MutationsOnREF!="" ,]
  cat("\n\t ", nrow(LFR_withMutations_df)," Long Fragment Reads with at least one mutations")

  cat("\n\t Counting the mutations per LFR")
  LFR_withMutations_df["NbOnSNP"] = unlist(lapply(as.character(unlist(LFR_withMutations_df$MutationsOnSNP)), function(x) if(x!="") length(unlist(strsplit(x,":"))) else 0))
  LFR_withMutations_df["NbOnREF"] = unlist(lapply(as.character(unlist(LFR_withMutations_df$MutationsOnREF)), function(x) if(x!="") length(unlist(strsplit(x,":"))) else 0))
  LFR_withMutations_df["NbMutations"] =  LFR_withMutations_df["NbOnSNP"] +  LFR_withMutations_df["NbOnREF"]


  LFR_withMutations_df =LFR_withMutations_df[LFR_withMutations_df$NbMutations>1,]
  cat("\n\t ", nrow(LFR_withMutations_df)," Long Fragment Reads with at least two mutations\n ")



  i <- sapply(LFR_withMutations_df, is.factor)
  LFR_withMutations_df[i] <- lapply(LFR_withMutations_df[i], as.character)


  #Vector of LFR Flag/ We will mark all the treated LFR

  nbLFR=nrow(LFR_withMutations_df)
  markedLFR<-vector("logical", nbLFR)
  markedLFR<- rep(FALSE,nbLFR)
  names(markedLFR)<-as.character(unlist(LFR_withMutations_df["LFR_name"]))


  cat("\n\nE  Extracting the phasing codes from ", nbLFR, " LFR ...")

  LFR_lst= rownames(LFR_withMutations_df)

  iphase=0

  stepmarker=max(1,length(LFR_lst)%/%100) #we take max to prevent the value to be 0.


  for (ilfr in 1:length(LFR_lst)){
    if((ilfr-1)%%stepmarker==0)
      cat("\n LFR No ", ilfr, "step ", (ilfr-1)%/%stepmarker +1,  " out of ", 100)

    myLFR=LFR_lst[ilfr]
    # cat("\n>",myLFR," value", markedLFR[myLFR])
    if(markedLFR[myLFR]){
      # cat(" Treated")
      next

    }


    iphase=iphase+1
    phasingcode<-paste("Phase_", iphase,sep="")
    SNPflag=1
    markedLFR =  process_fragments(LFR_withMutations_df,myLFR,phasingcode, 1,markedLFR,list_of_mutations_to_phase#,MutationPhasingCode_df
    )

    #  if(ilfr==3) stop()
  }

  MutationPhasingCode_df
}





#' Compute the Phaisng Code
#'
#'This function computes the phasing codes of each mutations. Mutations with the same phasing code are considered as phased (e.g PH_101_0 and PH_101_0) . Mutations with the same phasing fragment code but different allele index are considered as not phased(e.g PH_101_1 and PH_101_0) .
#'
#' @export
getPhaseInformation<-function(LFRset,MutationsAlleInfoset, calltype=NULL, OnlyPassFilter = FALSE  )
{

  LFR_df= LFRset
  Mutations_df=MutationsAlleInfoset

  cat("\n\nComputing the list of mutations called on each LFR as SNP...")
  LFR_df["MutationsOnSNP"]=  unlist(lapply(as.character(unlist(LFR_df$LFR_name)), function(x) paste(getMutationsOfLFR(LFR_df, Mutations_df,x,"SNP" ),collapse=":")))
  cat("\n\nComputing the list of mutations called on each LFR as REF...")
  LFR_df["MutationsOnREF"]=  unlist(lapply(as.character(unlist(LFR_df$LFR_name)), function(x) paste(getMutationsOfLFR(LFR_df, Mutations_df,x,"REF" ),collapse=":")))



  #getMutationsOfLFR<-function(LFR_df, Mutations_df, LFRname, calltype=NULL)
  #ComputePhasingCode<-function(LFR_withMutations_df,Mutations_set)
}



#' Retrieve the  mutations on each LFR of an LFRSet
#'
#'This function computes the list of mutations present on a  Long Fragment Reads either on their Variant Allele or their reference Allele
#'
#' @export
getMutationsOfLFRSet<-function(LFRset, MutationsAlleInfoset, progress=F)
{
  LFR_df= LFRset
  Mutations_df=MutationsAlleInfoset
if(progress){
  nbLFR=nrow(LFR_df)
  step=round(nbLFR/100)
}
  cat("\n\nComputing the list of mutations called on each LFR ...")
  for(irow in 1:nrow(LFR_df)){
    LFR_df[irow,"MutationsOnSNP"] = paste(getMutationsOfLFR(LFR_df[irow,], Mutations_df,"SNP" ),collapse=":")
    LFR_df[irow,"MutationsOnREF"] = paste(getMutationsOfLFR(LFR_df[irow,],  Mutations_df,"REF" ),collapse=":")
    if(progress)
      if(irow %% step ==0)
        cat("\n",round(irow/step),"- Processing LFR ", irow )
  }

  LFR_df


}




#' Compute the Phaisng Code
#'
#'This function computes the phasing codes of each mutations. Mutations with the same phasing code are considered as phased (e.g PH_101_0 and PH_101_0) . Mutations with the same phasing fragment code but different allele index are considered as not phased(e.g PH_101_1 and PH_101_0) .
#'
#' @export
ComputePhasingCode_iterative<-function(LFR_withMutations_df,Mutations_set,RetrieveHomo=FALSE)
{


  if(RetrieveHomo){
    cat("\n retrieve any potential homozygote mutations with zero Ref Allele count)")
    cat("\n\t ", nrow(Mutations_set[Mutations_set$WR==0, ]), " Potential homozygous mutations ")
    Mutations_set =Mutations_set[Mutations_set$WR!=0, ]

  }
  MutationsAlleleInfo_df=Mutations_set
  FullLFR_df = LFR_withMutations_df


  cat("\n\n Remove mutations with no wells count supporting the variant...")
  #This is odd, there are some mutations with no well count supporting the variant? We remove them from now.
  MutationsAlleleInfo_df =MutationsAlleleInfo_df[MutationsAlleleInfo_df$WV!=0, ]
  cat("\n ", nrow(MutationsAlleleInfo_df)," Mutations left..")



  cat("\n Select only Filter mutations ")

 # if(Filter=="PASS"){
    cat("\n\t  Select those with PASS filter ")
    MutationsPhase_df=MutationsAlleleInfo_df
    MutationsPhase_df =MutationsPhase_df[MutationsPhase_df$Filter=="PASS",]
    cat("\n ", nrow(MutationsPhase_df)," Mutations to phase..")
#   } else if(Filter=="NoBadreads"){
#     cat("\n\n Generating Badreads flag")
#     MutationsAlleleInfo_df$Badreads = unlist(lapply(as.character(MutationsAlleleInfo_df$Filter), function(x) as.numeric(grepl("badReads",x))))
#
#     cat("\n\t  Select those with No Bad reads  ")
#     MutationsPhase_df=MutationsAlleleInfo_df
#     MutationsPhase_df =MutationsPhase_df[MutationsPhase_df$Badreads==0,]
#     cat("\n ", nrow(MutationsPhase_df)," Mutations to phase..")
#   }else{
#     stop("\n\n Filter should be either PASS either No Badreads")
#   }




  MutationsPhase_df["PhasingCode"] = rep("",nrow(MutationsPhase_df))
  MutationsPhase_df["Homo"] = rep(0,nrow(MutationsPhase_df))

  MutationsPhase_df["LongPhasingCode"] = rep("",nrow(MutationsPhase_df))
  MutationsPhase_df["PhasingCodeTrace"] = rep("",nrow(MutationsPhase_df))
  MutationsPhase_df["ambiguous"] = rep(0,nrow(MutationsPhase_df))



  alternate_phasingcode<-function(phasingcode){
    index=substr(phasingcode,nchar(phasingcode),nchar(phasingcode))
    if(index=="1"){
      altindex="0"
    }else if (index=="0"){
      altindex="1"
    }else{
      cat("\n Code : ", phasingcode)
      stop("Error, last digit :  of the code should be 1 or 0")
    }

    altcode=paste(substr(phasingcode,1,nchar(phasingcode)-1),altindex,sep="")
    altcode
  }

  for(imut in 1 :nrow(MutationsPhase_df))
  {

    #######  print(MutationsPhase_df[1:20, 16:20])



    if(imut%%1000==0) cat ("\n\n - ", imut)

    code=paste("Phase_",imut,sep="")
    code_same=paste(code,"_1",sep="")
    code_other=paste(code,"_0",sep="")

    mutation=rownames(MutationsPhase_df[imut,])

    mutation_position = as.numeric(unlist(strsplit(mutation,"_"))[2])
    if(!is.na(MutationsPhase_df[imut,"Homo"]) && MutationsPhase_df[imut,"Homo"]>0)
      next



    if(MutationsPhase_df[imut,"PhasingCode"]==""){
      MutationsPhase_df[imut,"PhasingCode"] = code_same
    }else
      MutationsPhase_df[imut,"PhasingCode"] = paste(MutationsPhase_df[imut,"PhasingCode"] ,"-",code_same,"-",sep="")


    if ( MutationsPhase_df[imut,"LongPhasingCode"] ==""){
      MutationsPhase_df[imut,"LongPhasingCode"] = code_same
    }else
    {
      MutationsPhase_df[mutation,"PhasingCodeTrace"]  =  MutationsPhase_df[mutation,"LongPhasingCode"]
    }


    #List of LFR with the mutation on SNP
    LFR_withmutation_onSNP= FullLFR_df[grepl(mutation,FullLFR_df$MutationsOnSNP) , ]

    #list of LFR with the mutation on REF
    LFR_withmutation_onREF= FullLFR_df[ grepl(mutation,FullLFR_df$MutationsOnREF), ]


    coLFR_mutations_sameAllele_onSNP=c()
    coLFR_mutations_OtherAllele_onSNP = c()
    coLFR_mutations_sameAllele_onREF = c()
    coLFR_mutations_OtherAllele_onREF = c()




    if(nrow(LFR_withmutation_onSNP)>0){
      coLFR_mutations_sameAllele_onSNP= setdiff(intersect(unique(unlist(strsplit(paste(as.character(unlist( LFR_withmutation_onSNP$MutationsOnSNP)),collapse=":"),":"))), rownames(MutationsPhase_df)),mutation)
      coLFR_mutations_OtherAllele_onSNP=intersect(unique(unlist(strsplit(paste(as.character(unlist( LFR_withmutation_onSNP$MutationsOnREF)),collapse=":"),":"))), rownames(MutationsPhase_df))
    }

    if(nrow(LFR_withmutation_onREF)>0){
      coLFR_mutations_sameAllele_onREF= intersect(unique(unlist(strsplit(paste(as.character(unlist( LFR_withmutation_onREF$MutationsOnREF)),collapse=":"),":"))), rownames(MutationsPhase_df))
      coLFR_mutations_OtherAllele_onREF= intersect(unique(unlist(strsplit(paste(as.character(unlist( LFR_withmutation_onREF$MutationsOnSNP)),collapse=":"),":"))), rownames(MutationsPhase_df))
    }

    coLFR_mutations_sameAllele = unique(c(coLFR_mutations_sameAllele_onSNP,coLFR_mutations_sameAllele_onREF))
    coLFR_mutations_otherAllele = unique(c(coLFR_mutations_OtherAllele_onSNP,coLFR_mutations_OtherAllele_onREF ))

    homo=intersect(coLFR_mutations_sameAllele,coLFR_mutations_otherAllele)
    if(length(homo)>0){
      warning("\n Exhibits homozygocity")
      MutationsPhase_df[homo, "Homo"] =MutationsPhase_df[homo, "Homo"]  + 1

      coLFR_mutations_sameAllele = setdiff(coLFR_mutations_sameAllele,homo)
      coLFR_mutations_otherAllele=setdiff(coLFR_mutations_otherAllele, homo)
    }



    coLFR_mutations = c(coLFR_mutations_sameAllele ,  coLFR_mutations_otherAllele)

    if(length(coLFR_mutations)==0)
      next

    coLFR_mutations_df = as.data.frame(matrix(nrow=length(coLFR_mutations), ncol=2))
    colnames(coLFR_mutations_df) = c("Position","SameAllele")
    rownames(coLFR_mutations_df) = coLFR_mutations
    coLFR_mutations_df[coLFR_mutations_sameAllele,"SameAllele"] =  rep(1,length(coLFR_mutations_sameAllele))
    coLFR_mutations_df[coLFR_mutations_otherAllele,"SameAllele"] =  rep(0,length(coLFR_mutations_otherAllele))
    coLFR_mutations_df$Position = unlist(lapply(coLFR_mutations, function(x) as.numeric(unlist(strsplit(x,"_"))[2])))
    coLFR_mutations_df = coLFR_mutations_df[order(coLFR_mutations_df$Position),]

    #We keep only the mutations onward... coming after mutation_position



    coLFR_mutations_df=coLFR_mutations_df[coLFR_mutations_df$Position > mutation_position, ]

    if(nrow(coLFR_mutations_df)==0)
      next

    for (comut in 1:nrow(coLFR_mutations_df)){

      coLFRmut = rownames(coLFR_mutations_df[comut,])

      if(coLFR_mutations_df[comut,"SameAllele"]==1){
        if( MutationsPhase_df[coLFRmut,"PhasingCode"]==""){
          MutationsPhase_df[coLFRmut,"PhasingCode"] =code_same
        }else{
          MutationsPhase_df[coLFRmut,"PhasingCode"] = paste(MutationsPhase_df[coLFRmut,"PhasingCode"],code_same,sep=":" )
        }
        if(MutationsPhase_df[coLFRmut,"LongPhasingCode"] ==""){
          MutationsPhase_df[coLFRmut,"LongPhasingCode"] = MutationsPhase_df[mutation,"LongPhasingCode"]
        }else if(MutationsPhase_df[coLFRmut,"LongPhasingCode"] !=MutationsPhase_df[mutation,"LongPhasingCode"]){

          if(MutationsPhase_df[mutation,"LongPhasingCode"]==code_same)
            MutationsPhase_df[mutation,"LongPhasingCode"] =MutationsPhase_df[coLFRmut,"LongPhasingCode"]

          MutationsPhase_df[mutation,"PhasingCodeTrace"] = paste(MutationsPhase_df[mutation,"PhasingCodeTrace"], MutationsPhase_df[mutation,"LongPhasingCode"] ,sep=":")
          code_list= unlist(strsplit(MutationsPhase_df[mutation,"PhasingCodeTrace"],":"))
          if(length(code_list)>1)
            MutationsPhase_df[mutation,"ambiguous"] =1

          TopCode=sort(table(code_list),decreasing = T)[1]
          if(as.numeric(TopCode)>1)
          {
            MutationsPhase_df[mutation,"LongPhasingCode"] =names(TopCode)
          }
        }
      }else  if(coLFR_mutations_df[comut,"SameAllele"]==0){
        if( MutationsPhase_df[coLFRmut,"PhasingCode"]==""){
          MutationsPhase_df[coLFRmut,"PhasingCode"] =code_other
        }else{
          MutationsPhase_df[coLFRmut,"PhasingCode"] = paste(MutationsPhase_df[coLFRmut,"PhasingCode"],code_other,sep=":" )
        }
        if(MutationsPhase_df[coLFRmut,"LongPhasingCode"] ==""){
          MutationsPhase_df[coLFRmut,"LongPhasingCode"] = alternate_phasingcode(MutationsPhase_df[mutation,"LongPhasingCode"])
        }else if(MutationsPhase_df[coLFRmut,"LongPhasingCode"] !=alternate_phasingcode(MutationsPhase_df[mutation,"LongPhasingCode"])){
          if(MutationsPhase_df[mutation,"LongPhasingCode"]==code_same)
            MutationsPhase_df[mutation,"LongPhasingCode"] =alternate_phasingcode(MutationsPhase_df[coLFRmut,"LongPhasingCode"])

          MutationsPhase_df[mutation,"PhasingCodeTrace"] = paste(MutationsPhase_df[mutation,"PhasingCodeTrace"], alternate_phasingcode(MutationsPhase_df[mutation,"LongPhasingCode"]) ,sep=":")


          code_list= unlist(strsplit(MutationsPhase_df[mutation,"PhasingCodeTrace"],":"))
          if(length(code_list)>1)
            MutationsPhase_df[mutation,"ambiguous"] =1

          TopCode=sort(table(code_list),decreasing = T)[1]

          if(length(TopCode)>1)
            MutationsPhase_df[mutation,"ambiguous"] =1
          if(as.numeric(TopCode)>1)
          {
            MutationsPhase_df[mutation,"LongPhasingCode"] =names(TopCode)
          }

        }



      }else{
        stop(" \n How come I found myself here????")
      }


    }

  }






  MutationsPhase_df




}






#' Compute the Phaisng Code
#'
#'This function computes the phasing codes of each mutations. Mutations with the same phasing code are considered as phased (e.g PH_101_0 and PH_101_0) . Mutations with the same phasing fragment code but different allele index are considered as not phased(e.g PH_101_1 and PH_101_0) .
#'
#' @export
ComputePhasingCode<-function(LFR_withMutations_df,Mutations_set,mode="recursive")
{



  #select the level 1 mutations Germline mutations with Pass filer, shall with retrieve the homozygous?
#
#
#
#      if(mode=="iterative"){
#        cat("\n\n\t Computing the phasing using the iterative approach")
#        mutationphasing_df= ComputePhasingCode_iterative(LFR_withMutations_df, Mutations_set,RetrieveHomo=FALSE)
#      }else if(mode=="recursive"){
#        cat("\n\n\t Computing the phasing using the reciursive approach")
#        mutationphasing_df= ComputePhasingCode_recursive(LFR_withMutations_df, Mutations_set)
#      } else  if(mode=="both"){
#        cat("\n\n\t Computing the phasing using both methods")
#
#
#        cat("\n\n\t Computing the phasing using the iterative approach")
#        mutationphasing_df_iterative= ComputePhasingCode_iterative(LFR_withMutations_df, Mutations_set,RetrieveHomo=FALSE)
#
#        cat("\n\n\t Computing the phasing using the reciursive approach")
#        mutationphasing_df_recursive= ComputePhasingCode_recursive(LFR_withMutations_df, Mutations_set)
#
#      }





  AllMutations_df=Mutations_set
  #Level 1

  cat("\n\n\n Compute level 1 phasecode ...\n\n\t")

  #retrieve Germline Mutation
  if ("IsGermline" %in% colnames(Mutations_set)){
    Mutations_set=Mutations_set[Mutations_set$IsGermline==1,]
    cat("\n \n Germline mutations : ", nrow(Mutations_set)," mutations\n")
  }else if ("Germline" %in% colnames(Mutations_set)){
    Mutations_set=Mutations_set[Mutations_set$Germline==1,]
    cat("\n \n Germline mutations : ", nrow(Mutations_set)," mutations \n")
  }else{
    cat("\n\n No flag for the germline mutations in the data, all mutations will be considered as Germline")
  }

  #retrieve PASS Filter
  if("Filter" %in% colnames(Mutations_set)){
    Mutations_set =Mutations_set[Mutations_set$Filter=="PASS", ]
    cat("\n \n PASS mutations : ", nrow(Mutations_set)," mutations \n")
  }else{
    cat("\n No Filte rflag in the data, all the mutations will be considered as PASS mutations")
  }



  if(mode=="iterative"){
    cat("\n\n\t Computing the phasing using the iterative approach")
    mutationphasing_df_iterative= ComputePhasingCode_iterative(LFR_withMutations_df, Mutations_set,RetrieveHomo=FALSE)
  }  else  if(mode=="recursive"){
    cat("\n\n\t Computing the phasing using the reciursive approach")
    mutationphasing_df_recursive= ComputePhasingCode_recursive(LFR_withMutations_df, Mutations_set)
  } else  if(mode=="both"){
    cat("\n\n\t Computing the phasing using both methods")


    cat("\n\n\t Computing the phasing using the iterative approach")
    mutationphasing_df_iterative= ComputePhasingCode_iterative(LFR_withMutations_df, Mutations_set,RetrieveHomo=FALSE)

    cat("\n\n\t Computing the phasing using the reciursive approach")
    mutationphasing_df_recursive= ComputePhasingCode_recursive(LFR_withMutations_df, Mutations_set)

  }




  if(mode=="recursive" || mode=="both"){
    cat("\n\n\n\t Treating the recursive PhasingCode \n\t*********")

    mutationphasing_df_recursive = mutationphasing_df_recursive[!is.na(mutationphasing_df_recursive$PhasingCode1) & mutationphasing_df_recursive$PhasingCode1!="", ]
    AllMutations_df$PhasingCode1 = rep(NA,nrow(AllMutations_df))
    AllMutations_df$PhasingCode1_Level = rep(NA,nrow(AllMutations_df))
    AllMutations_df[rownames(mutationphasing_df_recursive),"PhasingCode1"] = mutationphasing_df_recursive$PhasingCode1
    AllMutations_df[rownames(mutationphasing_df_recursive),"PhasingCode1_Level"] = rep(1,nrow(mutationphasing_df_recursive))


    cat("\n\nRetrieve Phasing Code data ...")

    PhasingCode1= unique(as.character(unlist(mutationphasing_df_recursive$PhasingCode1)))
    PhasingCode1_wells <- vector("list", length(PhasingCode1))
    names(PhasingCode1_wells) <- PhasingCode1
    PhasingCode1_FirstPos=lapply(PhasingCode1_wells,function(x) Inf)
    PhasingCode1_LastPos=lapply(PhasingCode1_wells,function(x) 0)

    AllMutations_df$WV_IDs =  as.character(AllMutations_df$WV_IDs )
    for (imut in 1:nrow(mutationphasing_df_recursive) ){
      mut=rownames(mutationphasing_df_recursive)[imut]
      pos= mutationphasing_df_recursive[mut,"Pos"]
      phasecode= mutationphasing_df_recursive[mut,"PhasingCode1"]
      wells_ids=unlist(strsplit(AllMutations_df[mut,"WV_IDs"],":"))
      PhasingCode1_wells[phasecode] = list(unique(c(unlist(PhasingCode1_wells[phasecode]),wells_ids)))
      if(pos < unlist(PhasingCode1_FirstPos[phasecode]))
        PhasingCode1_FirstPos[phasecode] = pos
      if(pos > unlist(PhasingCode1_LastPos[phasecode]))
        PhasingCode1_LastPos[phasecode] = pos
    }


    ###Level 2

    cat("\n\nCompute level 2 phasecode ...\n\n\t")

    Level2_mutations=rownames(AllMutations_df[is.na(AllMutations_df$PhasingCode1),])
    for(imut in 1:length(Level2_mutations)){
      if(imut %%1000==0)
        cat(" ", imut )
      mut = Level2_mutations[imut]
      pos=AllMutations_df[mut,"Pos"]
      wells_ids=unlist(strsplit(AllMutations_df[mut,"WV_IDs"],":"))
      PhaseCandidate = unique(c(names(PhasingCode1_LastPos[PhasingCode1_LastPos<pos+1000000 & PhasingCode1_LastPos>pos-1000000 ]),
                                names(PhasingCode1_FirstPos[PhasingCode1_FirstPos<pos+1000000 & PhasingCode1_FirstPos>pos-1000000 ])))
      if(length(PhaseCandidate)>10){
        PhaseCandidate = unique(c(names(PhasingCode1_LastPos[PhasingCode1_LastPos<pos+500000 & PhasingCode1_LastPos>pos-500000 ]),
                                  names(PhasingCode1_FirstPos[PhasingCode1_FirstPos<pos+500000 & PhasingCode1_FirstPos>pos-500000 ])))
      }
      if(length(PhaseCandidate)>10){
        PhaseCandidate = unique(c(names(PhasingCode1_LastPos[PhasingCode1_LastPos<pos+100000 & PhasingCode1_LastPos>pos-100000 ]),
                                  names(PhasingCode1_FirstPos[PhasingCode1_FirstPos<pos+100000 & PhasingCode1_FirstPos>pos-100000 ])))
      }
      if(length(PhaseCandidate)>10){
        PhaseCandidate = unique(c(names(PhasingCode1_LastPos[PhasingCode1_LastPos<pos+50000 & PhasingCode1_LastPos>pos-50000 ]),
                                  names(PhasingCode1_FirstPos[PhasingCode1_FirstPos<pos+50000 & PhasingCode1_FirstPos>pos-50000 ])))
      }

      if(length(PhaseCandidate)==0)
        next

      number_intersecting_wells= unlist(lapply(PhaseCandidate, function(x) length(intersect(PhasingCode1_wells[x],wells_ids)) ))
      names(number_intersecting_wells) = PhaseCandidate
      number_intersecting_wells=sort(number_intersecting_wells,decreasing = TRUE)

      if(number_intersecting_wells[1] > 0){
        AllMutations_df[mut,"PhasingCode1"] = names(number_intersecting_wells[1])
        AllMutations_df[mut,"PhasingCode1_Level"] = 2

      }




    }

  }





  if(mode=="iterative" || mode=="both"){

    cat("\n\n\n\t the Treating the iterative PhasingCode \n\t*********")
    mutationphasing_df_iterative = mutationphasing_df_iterative[!is.na(mutationphasing_df_iterative$LongPhasingCode) & mutationphasing_df_iterative$LongPhasingCode!="" & mutationphasing_df_iterative$Homo==0, ]
    AllMutations_df$PhasingCode2 = rep(NA,nrow(AllMutations_df))
    AllMutations_df$PhasingCode2_Level = rep(NA,nrow(AllMutations_df))
    AllMutations_df[rownames(mutationphasing_df_iterative),"PhasingCode2"] = mutationphasing_df_iterative$LongPhasingCode
    AllMutations_df[rownames(mutationphasing_df_iterative),"PhasingCode2_Level"] = rep(1,nrow(mutationphasing_df_iterative))




    cat("\n\nRetrieve Phasing Code data ...")

    PhasingCode2= unique(as.character(unlist(mutationphasing_df_iterative$LongPhasingCode)))
    PhasingCode2_wells <- vector("list", length(PhasingCode2))
    names(PhasingCode2_wells) <- PhasingCode2
    PhasingCode2_FirstPos=lapply(PhasingCode2_wells,function(x) Inf)
    PhasingCode2_LastPos=lapply(PhasingCode2_wells,function(x) 0)

    AllMutations_df$WV_IDs =  as.character(AllMutations_df$WV_IDs )
    for (imut in 1:nrow(mutationphasing_df_iterative) ){
      mut=rownames(mutationphasing_df_iterative)[imut]
      pos= mutationphasing_df_iterative[mut,"Pos"]
      phasecode= mutationphasing_df_iterative[mut,"LongPhasingCode"]
      wells_ids=unlist(strsplit(AllMutations_df[mut,"WV_IDs"],":"))
      PhasingCode2_wells[phasecode] = list(unique(c(unlist(PhasingCode2_wells[phasecode]),wells_ids)))
      if(pos < unlist(PhasingCode2_FirstPos[phasecode]))
        PhasingCode2_FirstPos[phasecode] = pos
      if(pos > unlist(PhasingCode2_LastPos[phasecode]))
        PhasingCode2_LastPos[phasecode] = pos
    }


    ###Level 2

    cat("\n\nCompute level 2 phasecode ...\n\n\t")

    Level2_mutations=rownames(AllMutations_df[is.na(AllMutations_df$PhasingCode2),])
    for(imut in 1:length(Level2_mutations)){
      if(imut %%1000==0)
        cat(" ", imut )
      mut = Level2_mutations[imut]
      pos=AllMutations_df[mut,"Pos"]
      wells_ids=unlist(strsplit(AllMutations_df[mut,"WV_IDs"],":"))
      PhaseCandidate = unique(c(names(PhasingCode2_LastPos[PhasingCode2_LastPos<pos+1000000 & PhasingCode2_LastPos>pos-1000000 ]),
                                names(PhasingCode2_FirstPos[PhasingCode2_FirstPos<pos+1000000 & PhasingCode2_FirstPos>pos-1000000 ])))
      if(length(PhaseCandidate)>10){
        PhaseCandidate = unique(c(names(PhasingCode2_LastPos[PhasingCode2_LastPos<pos+500000 & PhasingCode2_LastPos>pos-500000 ]),
                                  names(PhasingCode2_FirstPos[PhasingCode2_FirstPos<pos+500000 & PhasingCode2_FirstPos>pos-500000 ])))
      }
      if(length(PhaseCandidate)>10){
        PhaseCandidate = unique(c(names(PhasingCode2_LastPos[PhasingCode2_LastPos<pos+100000 & PhasingCode2_LastPos>pos-100000 ]),
                                  names(PhasingCode2_FirstPos[PhasingCode2_FirstPos<pos+100000 & PhasingCode2_FirstPos>pos-100000 ])))
      }
      if(length(PhaseCandidate)>10){
        PhaseCandidate = unique(c(names(PhasingCode2_LastPos[PhasingCode2_LastPos<pos+50000 & PhasingCode2_LastPos>pos-50000 ]),
                                  names(PhasingCode2_FirstPos[PhasingCode2_FirstPos<pos+50000 & PhasingCode2_FirstPos>pos-50000 ])))
      }

      if(length(PhaseCandidate)==0)
        next

      number_intersecting_wells= unlist(lapply(PhaseCandidate, function(x) length(intersect(PhasingCode2_wells[x],wells_ids)) ))
      names(number_intersecting_wells) = PhaseCandidate
      number_intersecting_wells=sort(number_intersecting_wells,decreasing = TRUE)

      if(number_intersecting_wells[1] > 0){
        AllMutations_df[mut,"PhasingCode2"] = names(number_intersecting_wells[1])
        AllMutations_df[mut,"PhasingCode2_Level"] = 2

      }




    }

  }









  AllMutations_df



}



#' Extraction the copy number information
#'
#'This function extract at each mutation locus from the Copy Number table outputed by the software OncoSNP
#'
#' @export
getCopyNumber<-function(OncoSNP_copynumber,mutations_table){
  cat("\n\t \t Creating empty matrices\n\n")
  CopyNumber_df<-matrix(nrow=nrow(mutations_table), ncol=ncol(OncoSNP_copynumber) )
  CopyNumber_df<-as.data.frame(CopyNumber_df)
  colnames(CopyNumber_df) <-c("Chromosome","EndPosition","CopyNumber","LOH","Rank","Loglik","nProbes","NormalFraction","TumourState","PloidyNo","MajorCopyNumber","MinorCopyNumber","LowCopyNumber")
  rownames(CopyNumber_df) = rownames(mutations_table)

  list_column<-c("Chromosome","EndPosition","CopyNumber","LOH","Rank","Loglik","nProbes","NormalFraction","TumourState","PloidyNo","MajorCopyNumber","MinorCopyNumber")

  CopyNumber_df[c("Chromosome","EndPosition")]<-mutations_table[c(1,3)] # Chromosome and end position

  imutation=0


  for (imut in 1:nrow(mutations_table))
  {
    imutation=imutation+1
    mut <- rownames(mutations_table[imut ,])

    step<-(endmut-startmut)%/%100
    if(step==0) step=1
    if (imutation%%step==0) cat("step :  ", floor(imutation/step), " ")

    mutcn_df<-OncoSNP_copynumber[ OncoSNP_copynumber$StartPosition<as.numeric(mutations_table[mut,"Pos"]) & OncoSNP_copynumber$EndPosition>as.numeric(mutations_table[mut ,"Pos"]) , ]


    if(nrow(mutcn_df)==0)
      next

    mutcn <- mutcn_df[which.max(as.vector(unlist(mutcn_df["Rank"]))),list_column]
    mut_cn_df_plot<-mutcn_df[mutcn_df$Rank<3,]
    mutcn_plot <- mut_cn_df_plot[which.max(as.vector(unlist(mut_cn_df_plot["Rank"]))),"CopyNumber"]

    CopyNumber_df[mut,3:13]= c(mutcn[3:12],  mutcn_plot )


  }

  CopyNumber_df



}








#' @export
get_allbutlast<-function(x)
{
  x = as.character(unlist(x))
  substr(x, 1, nchar(x)-1)
}

#' @export
is_phased<-function(mutation1, mutation2, phased_submatrix,minsample)
{
  phased_str=""
  phased_1=0
  phased_0=0
  phased_NA=0


  for (isample in 1: length(colnames(phased_submatrix)))
  {
    sample=colnames(phased_submatrix)[isample]
    name_sample=sample

    if(is.na(phased_submatrix[mutation1,sample])) next
    if(is.na(phased_submatrix[mutation2,sample])) next

    if (phased_submatrix[mutation1,sample]==phased_submatrix[mutation2,sample])
    {
      phased_1=phased_1+1
      phased_str=paste(phased_str,paste(1,name_sample,sep=":"),sep="|")

      #  cat(" ",phased_submatrix[mutation1,sample])
    }else if((substr(phased_submatrix[mutation1,sample],1,6)=="Phased") && (substr(phased_submatrix[mutation2,sample],1,6)=="Phased"))
    {
      if (get_allbutlast(phased_submatrix[mutation1,sample])==get_allbutlast(phased_submatrix[mutation2,sample]))
      {
        phased_0=phased_0+1
        phased_str=paste(phased_str,paste(0,name_sample,sep=":"),sep="|")
      }
    }


  }

  is_phased=NA

  phased_occurence = paste(" 1:", phased_1, " 0:", phased_0,sep="")

  if( phased_1<minsample && phased_0<minsample)
  {
    is_phased=NA
  }else{
    if(phased_1>phased_0)
      is_phased=1
    else if(phased_0>phased_1)
      is_phased=0

    #more than 3 sample support every alternative, then there is an ambiguity.  report NA
    if((phased_1>3) && (phased_0>3))
      is_phased=NA
  }



  list(phased=is_phased, occurence= phased_occurence, str=phased_str)

}







#Similar like Is phased but in a more strictly manner, phased if there is absolutely no evidence supporting the not phased alternative
#' @export
is_strictlyphased<-function(mutation1, mutation2, phased_submatrix,minsample)
{
  phased_str=""
  phased_1=0
  phased_0=0
  phased_NA=0


  for (isample in 1: length(colnames(phased_submatrix)))
  {
    sample=colnames(phased_submatrix)[isample]
    name_sample=sample

    if(is.na(phased_submatrix[mutation1,sample])) next
    if(is.na(phased_submatrix[mutation2,sample])) next

    if (phased_submatrix[mutation1,sample]==phased_submatrix[mutation2,sample])
    {
      phased_1=phased_1+1
      phased_str=paste(phased_str,paste(1,name_sample,sep=":"),sep="|")

      #  cat(" ",phased_submatrix[mutation1,sample])
    }else if((substr(phased_submatrix[mutation1,sample],1,6)=="Phased") && (substr(phased_submatrix[mutation2,sample],1,6)=="Phased"))
    {
      if (get_allbutlast(phased_submatrix[mutation1,sample])==get_allbutlast(phased_submatrix[mutation2,sample]))
      {
        phased_0=phased_0+1
        phased_str=paste(phased_str,paste(0,name_sample,sep=":"),sep="|")
      }
    }


  }

  is_phased=NA

  phased_occurence = paste(" 1:", phased_1, " 0:", phased_0,sep="")

  if( phased_1<minsample && phased_0<minsample)
  {
    is_phased=NA
  }else if ( phased_1>0 && phased_0>0){
    is_phased=NA
  }else if(phased_1>0)
    is_phased=1
  else if(phased_0>0)
    is_phased=0

  #more than 3 sample support every alternative, then there is an ambiguity.  report NA
  #if((phased_1>3) && (phased_0>3))
  #  is_phased=NA
  #}



  list(phased=is_phased, occurence= phased_occurence, str=phased_str)

}

#' @export
list_phasing<-function(mutation, phased_submatrix,minsample,mode="strict")
{

  phased=vector()
  not_phased=vector()
  phasing_stroccurence=""
  phasing_strstr=""
  for (jmut in 1: nrow(phased_submatrix))
  {
    #  mutationj=  rownames(phased_submatrix[jmut,])
    mutationj=  rownames(phased_submatrix[jmut,colnames(phased_submatrix),drop=FALSE])

    if(mutationj==mutation) next
    # cat("\n mutationj is ",mutationj)
    if(!is.null(mode) && (mode=="strict")){
      phasing_result= is_strictlyphased(mutation,mutationj, phased_submatrix,minsample)
    }else{
      phasing_result= is_phased(mutation,mutationj, phased_submatrix,minsample)
    }

    # cat (" here it is ", phasing_result$phased)

    if (!is.na(phasing_result$phased))
    {

      #  cat (" here it is "); print(phasing_result)

      phasing_stroccurence=paste(phasing_stroccurence,phasing_result$occurence,sep="/")
      phasing_strstr=paste(phasing_strstr,phasing_result$str,sep="/")
      if (phasing_result$phased==1)
        phased<-c(phased,mutationj)
      else if (phasing_result$phased==0)
        not_phased<-c(not_phased,mutationj)

    }

  }

  list(phased_list=phased,notphased_list=not_phased, phased_occurence= phasing_stroccurence,phased_str=phasing_strstr)

}




#' Extraction phasing association
#'
#'This function take a matrix or submatrix and compute the phasing association between the somatic mutations and the germline mutations, for each somatic mutations present in the matrix its compute its list of phased germline mutations
#'
#' @export
GetPhasingAssociation<-function(phasedmatrix, FromType, ToType, maxdistance, number_node,totalparts_node,minsample, cifsm,mode="strict")
{



  #For each mutations of type "FromType" we compute the list of mutation of type "ToType" phased to it. Type cand be any of the following three : "somatic", "germline", 'any"

  cilsm=ncol(phasedmatrix)


  if (!(FromType %in% c("somatic","germline","any") && ToType %in%  c("somatic","germline","any") ))
    stop(" FromType and ToType must be one  of the following : 'somatic', 'germline' or  'any' ")

  if(FromType=="somatic") from_germline_status=0
  if(FromType=="germline") from_germline_status=1

  if(ToType=="somatic") to_germline_status=0
  if(ToType=="germline") to_germline_status=1

  #We restrict the rows
  cat("\n\n restrict the matrix and prepare the list of ")
  phasedmatrix =phasedmatrix[rowSums(!is.na(phasedmatrix[cifsm:cilsm]))>0,]
  fromtypephased_chr=phasedmatrix[1:(cifsm-1)]
  masterfromtypephased_chr=phasedmatrix
  #Important, we extract the submatrix only for the From Type mutations
  if(FromType!="any")
  {
    fromtypephased_chr  = fromtypephased_chr[fromtypephased_chr$IsGermline==from_germline_status, ]
    masterfromtypephased_chr=masterfromtypephased_chr[masterfromtypephased_chr$IsGermline==from_germline_status, ]
  }
  if(ToType!="any")
    phasedmatrix= phasedmatrix[phasedmatrix$IsGermline==to_germline_status ,]

  fromtypephased_chr["Phased_List"] = rep("", nrow(fromtypephased_chr))
  fromtypephased_chr["NotPhased_List"] = rep("", nrow(fromtypephased_chr))




  cat("\n\t*******************\n Extract position to run \n\t *********************")
  start_position_node=1
  end_position_node = nrow(fromtypephased_chr)
  if (number_node !=0)
  {
    cat("\n\n\t We extract the ", number_node,"th part out of ",totalparts_node,"  parts")
    full_length=nrow(fromtypephased_chr)
    part_length=ceiling(full_length/totalparts_node)
    start_position_node= (number_node-1) * part_length + 1
    end_position_node = start_position_node + part_length -1
    if (end_position_node >nrow(fromtypephased_chr) )
      end_position_node = nrow(fromtypephased_chr)

    cat("\n\n\t The run will go from mutation  ",start_position_node," to mutation ", end_position_node)
  }else
  {
    cat("\n\t No subpart extracted for this node, the whole chromosome considered...")
  }

  startmut=start_position_node
  endmut = end_position_node
  #fromtypephased_chr=fromtypephased_chr[startmut:endmut, ]
  #startmut=start_position_node + 600
  imutation=0
  for (imut in startmut:endmut)
  {
    imutation =imutation+1;
    mut <- rownames(fromtypephased_chr[imut,])
    mut_pos=as.numeric(fromtypephased_chr[imut,"Pos"])


    output=(((imutation-1) %%1)==0)
    if(output)
      cat (" \n mut ",imutation,"/ ",endmut-startmut," (", mut,") :  ",imut,"/",endmut,"/",nrow(fromtypephased_chr) )
    #cat("\n retrieval of ", ToType," mutations")

    minpos <- mut_pos -maxdistance;  maxpos<-  mut_pos +maxdistance
    subToType_phasing=phasedmatrix[as.numeric(phasedmatrix$Pos)>minpos & as.numeric(phasedmatrix$Pos)<maxpos & !is.na(phasedmatrix$Pos) ,cifsm:cilsm]
    nb_ToTypes = nrow(subToType_phasing)
    if(output)   cat ("  ", nb_ToTypes, " ", ToType, " mutations located within ", maxdistance, " bp  "  )

    # submatrix_phasing = rbind(phasedmatrix[mut,cifsm:cilsm], subToType_phasing)



    #If the number of found mutations is to high, then  we consider mutations located within half the previous distance
    #If a still  really high number of mutations around (>1000), then we choose only the mutations not close by less than 50 bp form each other.

    Reduce=T
    if(Reduce){
      #if (length(unlist(phasing_info$phased_list))>100)
      maxdistance_reduce = maxdistance
      while (nb_ToTypes>2000){
        maxdistance_reduce =  maxdistance_reduce/2
        minpos <- mut_pos -maxdistance_reduce;  maxpos<-  mut_pos +maxdistance_reduce
        subToType_phasing=phasedmatrix[as.numeric(phasedmatrix$Pos)>minpos & as.numeric(phasedmatrix$Pos)<maxpos & !is.na(phasedmatrix$Pos)  ,cifsm:cilsm]
        nb_ToTypes = nrow(subToType_phasing)

        if(output)   cat ("\n\t Too high we go with   ", nb_ToTypes, " ", ToType, "  located within ", maxdistance_reduce, " bp" )

      }



    }



    cat("\n level 1")

    # submatrix_phasing = rbind(phasedmatrix[mut,cifsm:cilsm], subToType_phasing)
    submatrix_phasing = rbind(masterfromtypephased_chr[mut,cifsm:cilsm], subToType_phasing)

    #Only sample where the somatic has a phasing information
    somatic_phasing=submatrix_phasing[mut,]
    samples_toconsider=names( somatic_phasing[,!is.na( somatic_phasing),drop=F])
    submatrix_phasing=submatrix_phasing[samples_toconsider]

    #Only germline  with phasing information
    not_naphasing=rowSums(!is.na(submatrix_phasing))
    germline_withphase=names(not_naphasing[not_naphasing>0])
    submatrix_phasing= submatrix_phasing[germline_withphase,samples_toconsider, drop=F]
    # cat("\n")
    # print(head(submatrix_phasing))


    if(output)   cat ("\n\t   ", nrow(submatrix_phasing), " germline mutations "  )



    if (nb_ToTypes>0)
      phasing_info = list_phasing(mut,submatrix_phasing, minsample,mode=mode )
    else
      next



    if(output)  cat ("\n",length(unlist(phasing_info$phased_list)), "Phased, ", length(unlist(phasing_info$notphased_list)), " Not phased \n Phasing string : ",phasing_info$phased_str )

    if (length(unlist(phasing_info$phased_list))==0)
      next

    cat("\n level 2")

    phased_list<-unlist(phasing_info$phased_list);  notphased_list<-unlist(phasing_info$notphased_list)
    phased_list = setdiff(phased_list,notphased_list); notphased_list=setdiff(notphased_list,phased_list)

    if(output)  cat (length( phased_list), " non conflicting Phased, ", length(notphased_list), "non conflicting  Not phased")

    fromtypephased_chr[mut , "Phased_List"] = paste( phased_list, collapse=":")
    fromtypephased_chr[mut, "NotPhased_List"] = paste(notphased_list, collapse=":")
  }

  if (totalparts_node==1)
    write.csv(fromtypephased_chr_nothead, file =paste(Patient_Folder,"Phasing/Chromosomes/phasingassociation_",FromType,"_",ToType,"_",Patient_Name,"_",chrom, ".csv",sep=""),quote=F)

  fromtypephased_chr[startmut:endmut,]
}














