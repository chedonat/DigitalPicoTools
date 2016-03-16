

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
plot_fragments <- function(fragment_df, type = "ScatterPlot", value = "Coverage", minlength = 10000, color = "black", main = "", well = NULL) {

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
        if (value == "Coverage") {
            myplot <- ggplot(fragment_df, aes(Length, AvgCoverage)) + geom_point(colour = color) + ggtitle(main) + ylab("Average Coverage") + xlab("Fragment length") +
                theme(axis.text = element_text(size = 12, face = "bold"), axis.title = element_text(size = 12, face = "bold"), plot.title = element_text(size = 14,
                  face = "bold"), legend.title = element_text(size = 10, face = "bold"))

        } else if (value == "Length") {

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

        if (value == "Coverage") {
            myplot <- ggplot(frag_df, aes(x = Pos_offset, AvgCoverage, size = factor(Point_size1))) + geom_point(colour = color) + ggtitle(main) +
                geom_vline(xintercept = chrom_offset[1:25]) + theme(axis.text = element_text(size = 14, face = "bold", color = "white"), axis.title = element_text(size = 14,
                face = "bold"), plot.title = element_text(size = 16, face = "bold"), legend.title = element_text(size = 12, face = "bold")) + scale_size_discrete(name = "Fragment \n Length",
                breaks = size_list[indices], labels = paste(size_list[indices], "0K", sep = "")) + ylab("Average Coverage") + xlab("Whole Genome") +
                ylim(0, 1.1 * max(as.numeric(frag_df$AvgCoverage))) + annotate("text", x = chrom_offset[1:25] + (unlist(hg19_chrsize[1:25]) + c(rep(0,
                24), 5e+07))/2, y = 0, label = c(1:22, "X", "Y", "M"))

        } else if (value == "Length") {
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

    myplot


}




#' @export
hist_fragments <- function(fragment_df, type = "ScatterPlot", value = "Coverage", minlength = 10000, color = "black", main = "", well = NULL) {

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

    if (value == "Coverage") {

      myplot <- ggplot(data =fragment_df, aes(x = AvgCoverage)) +
        geom_histogram(
          fill = rep(color,33),
          color = "dodgerblue2")  + ggtitle(main) +
        ylab("Frequency") + xlab("Average Coverage")+ theme(axis.text=element_text(size=12,face="bold"),  axis.title=element_text(size=12,face="bold"),  plot.title = element_text(size=16,face="bold"),  legend.title=element_text(size=10,face="bold"))+
        scale_x_continuous(breaks=number_ticks(10)) +
        scale_y_continuous(breaks=number_ticks(10))

    } else if (value == "Length") {

      myplot <- ggplot(data =fragment_df, aes(x = Length)) +
        geom_histogram(
          fill = rep(color,33),
          color = "dodgerblue2")  + ggtitle(main) +
        ylab("Frequency") + xlab("Fragments Length")+ theme(axis.text=element_text(size=12,face="bold"),  axis.title=element_text(size=12,face="bold"),  plot.title = element_text(size=16,face="bold"),  legend.title=element_text(size=10,face="bold"))+
        scale_x_continuous(breaks=number_ticks(10)) +
        scale_y_continuous(breaks=number_ticks(10))


    } else {
      stop("\nValue should be either Coverage either Length")
    }




  print(myplot)

  myplot



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
plot_AlleleInfo <- function(variant_df, region = NULL, samplingRatio = 1, value = "WellsFraction", main = "") {

if(samplingRatio>1)
  stop("\n the argument sampling ration should be < 1")


  hg19_chrsize <- list(chr1 = 249250621, chr2 = 243199373, chr3 = 198022430, chr4 = 191154276, chr5 = 180915260, chr6 = 171115067, chr7 = 159138663,
                       chr8 = 146364022, chr9 = 141213431, chr10 = 135534747, chr11 = 135006516, chr12 = 133851895, chr13 = 115169878, chr14 = 107349540, chr15 = 102531392,
                       chr16 = 90354753, chr17 = 81195210, chr18 = 78077248, chr19 = 59128983, chr20 = 63025520, chr21 = 48129895, chr22 = 51304566, chrX = 155270560,
                       chrY = 59373566, chrM = 16571)


  if (!is.null(region)) {
    GRegion = parseregion(region)
    if (!is.null(GRegion$Chrom)) {
      variant_df = variant_df[variant_df$Chrom == GRegion$Chrom, ]
      mychrom = GRegion$Chrom
      startPosition = 1
      PosPosition = hg19_chrsize[mychrom]
      if (!is.null(GRegion$Start) && !is.null(GRegion$Pos)) {
        startPosition = GRegion$Start
        PosPosition = GRegion$Pos
        variant_df = variant_df[variant_df$Pos >= startPosition & variant_df$Pos <= PosPosition, ]
      }
    }
  }else{
    warnings("\n No region provided, The first chromosome is taken as the region \n ")
    mychrom<-unlist(variant_df["Chrom"])[1]
    cat("\n :  ", mychrom)
  }




  if (nrow(variant_df) == 0)
    stop("\n No variants in the specified region")

  step = floor(1/samplingRatio)

  indices = seq(1, nrow(variant_df), step)
  variant_df = variant_df[indices, ]



  master_df = variant_df
  master_df["Pos_offset"] = as.numeric(unlist(master_df["Pos"])) + as.numeric(unlist(lapply(as.character(unlist(master_df["Chrom"])), function(x) if (!is.na(x)) as.numeric(chrom_offset[x]) else NA)))



  if (value == "WellsFraction") {
    column_to_plot= "WF"
  } else if (Value == "WellsCount") {
    column_to_plot= "WR"
  } else if (value == "ReadsFraction") {
    column_to_plot= "AF"
  } else if (value == "ReadsCount") {
    column_to_plot= "TC"
  } else {
    stop("\n Invalid value parameter")
  }


  present_chromosome<-as.character(unlist(unique(variant_df["Chrom"])))
  present_chromosome=present_chromosome[!is.na(present_chromosome)]
  offset_chromosome<-list()
  offset_chromosome[present_chromosome[1]]=0

  if(length(present_chromosome)>1)
    for (ichr in 2: length(present_chromosome))
      offset_chromosome[present_chromosome[ichr]]=as.numeric(offset_chromosome[present_chromosome[ichr-1]])+as.numeric(hg19_chrsize[present_chromosome[ichr-1]])

  variant_df["Pos_offset"]=variant_df["Pos"] + as.numeric(unlist(lapply(as.character(unlist(variant_df["Chrom"])), function(x) if (!is.na(x)) as.numeric(offset_chromosome[x] ) else NA )))

  plot(unlist(variant_df["Pos_offset"]),unlist(variant_df[ column_to_plot]), cex=0.2,main=main ,ylab=value,xlab="",cex.main=1.5, col="black",ylim=c(0,1) )
  v0=as.numeric(offset_chromosome); l=length(v0); v0[l+1]= variant_df$Pos_offset[nrow(variant_df)];l=length(v0);
  v1=(v0[1:(l-1)] + v0[2:l])/2
  #chrom_indices=which(chrom_list==present_chromosome)
  abline(v =v0 ,lwd=2,col = "black"); mtext(present_chromosome,side=1,line=2, at=v1,las=2)





}
