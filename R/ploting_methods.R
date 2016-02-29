

#' @export
parseregion <- function(region) {

    chrom = NULL
    startPosition = NULL
    endPosition = NULL

    region_parts = unlist(strsplit(region, ":"))
    if (length(region_parts > 1)) {
        chrom = region_parts[1]


        if (length(region_parts) > 1) {
            coordinates = unlist(strsplit(region_parts[2], "-"))
            startPosition = as.numeric(coordinates[1])
            endPosition = as.numeric(coordinates[2])
        }

    }
    list(Chrom = chrom, Start = startPosition, End = endPosition)
}


#' @export
number_ticks <- function(n) {
    function(limits) pretty(limits, n)
}

#' @export
plot_fragments <- function(fragment_df, type = "ScatterPlot", Value = "Coverage", minlength = 10000, color = "black", main = "", well = NULL) {

    if (is.null(well)) {
        well = unique(fragment_df$well_name)
        if (length(well) > 0)
            warning("\n More than one well provided. All  will be considered. Restrict to a single well by setting the well parameter ")


    } else {
        fragment_df = fragment_df[fragment_df$well_name == well, ]
    }




    fragment_df = fragment_df[fragment_df$Length > minlength, ]

    if (nrow(fragment_df) == 0) {
        stop(" \n\n No fragments with the given parameters")
    }

    if (type == "ScatterPlot") {
        if (Value == "Coverage") {
            myplot <- ggplot(fragment_df, aes(Length, AvgCoverage)) + geom_point(colour = color) + ggtitle(main) + ylab("Average Coverage") + xlab("Fragment length") +
                theme(axis.text = element_text(size = 12, face = "bold"), axis.title = element_text(size = 12, face = "bold"), plot.title = element_text(size = 14,
                  face = "bold"), legend.title = element_text(size = 10, face = "bold"))

        } else if (Value == "Length") {

            myplot <- ggplot(fragment_df, aes(AvgCoverage, Length)) + geom_point(colour = color) + ggtitle(main) + ylab("Fragment length") + xlab(" Average Coverage") +
                theme(axis.text = element_text(size = 12, face = "bold"), axis.title = element_text(size = 12, face = "bold"), plot.title = element_text(size = 14,
                  face = "bold"), legend.title = element_text(size = 10, face = "bold"))


        } else {
            stop("\nValue should be either Coverage either Length")
        }

    } else if (type == "Landscape") {

        frag_df = fragment_df[fragment_df$Length > minlength, ]



        frag_df["Chrnum"] <- unlist(lapply(as.character(unlist(frag_df["Chrom"])), function(x) if (x == "chrX") {
            23
        } else if (x == "chrY") {
            24
        } else if (x == "chrM") {
            25
        } else as.numeric(substr(x, 4, nchar(x)))))

        frag_df["Pos_offset"] = as.numeric(unlist(frag_df["Start"])) + as.numeric(unlist(lapply(as.numeric(unlist(frag_df["Chrnum"])), function(x) if (!is.na(x)) as.numeric(chrom_offset[x]) else NA)))



        frag_df <- na.omit(frag_df)
        frag_df$Length = as.numeric(frag_df$Length)
        frag_df$Pos_offset = as.numeric(frag_df$Pos_offset)
        frag_df$AvgCoverage = as.numeric(frag_df$AvgCoverage)
        frag_df["Point_size1"] = floor(as.numeric(unlist(frag_df["Length"]))/10000)
        frag_df["Point_size2"] = floor(as.numeric(unlist(frag_df["AvgCoverage"])))


        size_list <- unique(sort(as.numeric(unlist(frag_df["Point_size1"]))))
        indices = seq(1, length(size_list), 2)

        if (Value == "Coverage") {
            myplot <- ggplot(frag_df, aes(x = Pos_offset, AvgCoverage, size = factor(Point_size1))) + geom_point(colour = color) + ggtitle(main) +
                geom_vline(xintercept = chrom_offset[1:25]) + theme(axis.text = element_text(size = 14, face = "bold", color = "white"), axis.title = element_text(size = 14,
                face = "bold"), plot.title = element_text(size = 16, face = "bold"), legend.title = element_text(size = 12, face = "bold")) + scale_size_discrete(name = "Fragment \n Length",
                breaks = size_list[indices], labels = paste(size_list[indices], "0K", sep = "")) + ylab("Average Coverage") + xlab("Whole Genome") +
                ylim(0, 1.1 * max(as.numeric(frag_df$AvgCoverage))) + annotate("text", x = chrom_offset[1:25] + (unlist(hg19_chrsize[1:25]) + c(rep(0,
                24), 5e+07))/2, y = 0, label = c(1:22, "X", "Y", "M"))

        } else if (Value == "Length") {
            myplot <- ggplot(frag_df, aes(x = Pos_offset, Length, size = factor(Point_size2))) + geom_point(colour = color) + ggtitle(main) + geom_vline(xintercept = chrom_offset[1:25]) +
                theme(axis.text = element_text(size = 14, face = "bold", color = "white"), axis.title = element_text(size = 14, face = "bold"),
                  plot.title = element_text(size = 16, face = "bold"), legend.title = element_text(size = 12, face = "bold")) + scale_size_discrete(name = "Average \n Coverage",
                labels = paste(unique(sort(as.numeric(unlist(frag_df["Point_size2"])))), "x", sep = "")) + ylab("Fragment Length") + xlab("Whole Genome") +
                annotate("text", x = chrom_offset[1:25] + (unlist(hg19_chrsize[1:25]) + c(rep(0, 24), 5e+07))/2, y = -2, label = c(1:22, "X", "Y",
                  "M"))

        }




    } else {
        stop(" Parameter type shuld be either ScatterPlot either Landscape")
    }


    print(myplot)



}

#' @export
plot_LFRChromosomeDistribution <- function(fragment_df, region = NULL, minLength = 10000) {


    chrom = NULL
    startPosition = NULL
    endPosition = NULL

    if (!is.null(region)) {
        GRegion = parseregion(region)
        if (!is.null(GRegion$Chrom)) {
            fragment_df = fragment_df[fragment_df$Chrom == GRegion$Chrom, ]
            chrom = GRegion$Chrom
            startPosition = 1
            endPosition = hg19_chrsize[chrom]
            if (!is.null(GRegion$Start) && !is.null(GRegion$End)) {
                startPosition = GRegion$Start
                endPosition = GRegion$End
                fragment_df = fragment_df[fragment_df$Start >= startPosition & fragment_df$End <= endPosition, ]
            }
        }
    }

    cat("region is ", region)

    if (nrow(fragment_df) == 0)
        stop("\n No Long Fragment Reads in the specified region ")

    # if(!is.null()) unlist(strsplit(region,':')) chrom = region_parts[1] startPosition = 1 endPosition = hg19_chrsize[chrom] fragment_df =
    # fragment_df[fragment_df$Chrom == chrom , ] if(length(region_parts)>1){ coordinates = unlist(strsplit(region_parts[2],'-')) startPosition =
    # as.numeric(coordinates[1]) endPosition = as.numeric(coordinates[2]) fragment_df = fragment_df[fragment_df$Chrom == chrom & fragment_df$Start
    # >= startPosition & fragment_df$End <= endPosition,] } if(chrom!=mychrom){ warning('Chromosome in region do not match chromosome parameter. The
    # chromosome specified in the region will be considered') } mychrom=chrom }else{ fragment_df = fragment_df[fragment_df$Chrom == mychrom , ] }


    regionsize = endPosition - startPosition

    # fragment_df= fragment_wgs_df[fragment_wgs_df$Chrom==mychrom,]
    fragment_df = fragment_df[!is.na(fragment_df$Length) & as.numeric(unlist(as.character(fragment_df$Length))) > minLength, ]


    fragment_df$AvgCoverage = as.numeric(as.character(unlist(fragment_df$AvgCoverage)))
    fragment_df$End = as.numeric(as.character(unlist(fragment_df$End)))
    fragment_df$Start = as.numeric(as.character(unlist(fragment_df$Start)))
    fragment_df$Well_Number = as.numeric(as.character(unlist(fragment_df$Well_Number)))

    maxavgcoverage = max(fragment_df$AvgCoverage, na.rm = T)
    maxwidth = 3

    y_blue = floor(median(unique(fragment_df$Well_Number)))
    x_blue = startPosition + regionsize/3

    if (regionsize > 5e+06) {
        length_blue = 1e+06
    } else {
        length_blue = 1e+05
    }

    fragment_df["Width"] = fragment_df$AvgCoverage * maxwidth/maxavgcoverage


    p1 <- ggplot(fragment_df, aes(Start, Well_Number)) + geom_segment(aes(x = Start, y = Well_Number, xend = End, yend = Well_Number), size = fragment_df$Width) +
        geom_segment(x = x_blue, y = y_blue, xend = x_blue + length_blue, yend = y_blue, color = "blue", size = 1) + theme(axis.text = element_text(size = 12,
        face = "bold"), axis.title = element_text(size = 12, face = "bold"), plot.title = element_text(size = 14, face = "bold"), legend.title = element_text(size = 10,
        face = "bold"))

    print(p1)

    list(region = paste(chrom, ":", startPosition, "-", endPosition, sep = ""), blueSegment_start = x_blue, blueSegment_length = length_blue, blueSegment_well = y_blue)
}







#' @export
plot_allfraction <- function(variant_df, region = NULL, samplingRatio = 1, Value = "WellsFraction", main = "") {


    if (!is.null(region)) {
        GRegion = parseregion(region)
        if (!is.null(GRegion$Chrom)) {
            variant_df = variant_df[variant_df$Chrom == GRegion$Chrom, ]
            mychrom = GRegion$Chrom
            startPosition = 1
            endPosition = hg19_chrsize[mychrom]
            if (!is.null(GRegion$Start) && !is.null(GRegion$End)) {
                startPosition = GRegion$Start
                endPosition = GRegion$End
                variant_df = variant_df[variant_df$Pos >= startPosition & variant_df$Pos <= endPosition, ]
            }
        }
    }






    # if(!is.null(region)){ region_parts= unlist(strsplit(region,':')) chrom = region_parts[1] startPosition = 1 endPosition = hg19_chrsize[chrom]
    # fragment_df = variant_df[variant_df$Chrom == chrom , ] if(length(region_parts)>1){ coordinates = unlist(strsplit(region_parts[2],'-'))
    # startPosition = as.numeric(coordinates[1]) endPosition = as.numeric(coordinates[2]) variant_df = variant_df[variant_df$Chrom == chrom &
    # variant_df$Start >= startPosition & variant_df$End <= endPosition,] } if(chrom!=mychrom){ warning('Chromosome in region do not match
    # chromosome parameter. The chromosome specified in the region will be considered') } mychrom=chrom }else{ variant_df =
    # variant_df[variant_df$Chrom == mychrom , ] }

    if (nrow(variant_df) == 0)
        stop("\n No variants in the specified region")

    step = floor(1/samplingRatio)

    indices = seq(1, nrow(variant_df), step)
    variant_df = variant_df[indices, ]



    master_df = variant_df
    master_df["Pos_offset"] = as.numeric(unlist(master_df["Pos"])) + as.numeric(unlist(lapply(as.character(unlist(master_df["Chrom"])), function(x) if (!is.na(x)) as.numeric(chrom_offset[x]) else NA)))



    if (Value == "WellsFraction") {
        column = "WF"
    } else if (Value == "WellsCount") {
        column = "WR"
    } else if (Value == "ReadsFraction") {
        column = "AF"
    } else if (Value == "ReadsCount") {
        column = "TC"
    } else {
        stop("\n Invalid Value parameter")
    }
    plot(unlist(master_df["Pos"]), unlist(master_df[column]), cex = 0.2, main = main, ylab = Value, cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.5,
        col = "black", xlab = NA, xlim = c(0, as.numeric(hg19_chrsize[mychrom])))
    v0 = as.numeric(chrom_offset)
    l = length(v0)
    v0[l + 1] = master_df$Pos_offset[nrow(master_df)]
    l = length(v0)
    v1 = (v0[1:(l - 1)] + v0[2:l])/2
    abline(v = v0, lwd = 2, col = "black")
    mtext(chrom_list, side = 1, line = 1, at = v1, las = 2)


}
