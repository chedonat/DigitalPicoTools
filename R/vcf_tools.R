
# Loading a VCF files
#' @export
loading_vcf <- function(vcf_file, samples_id = NULL, region = NULL) {

  cat("\n\n\t Loading the variants at   : ", vcf_file)
  variant_df <- read.table(file = vcf_file, header = F)


  # If samples_id not provided, following the VCF specifications, the columns names from field 9 to the last field are considered.
  #If samples_id is provided as a single integer n, the first n character of the samples_id column are considered.

  if (is.null(samples_id)) {
    samples_id=samples(scanVcfHeader(vcf_file))
  }else if(class(samples_id)=="numeric")
  {
    samples_names=samples(scanVcfHeader(vcf_file))
    samples_id = substr(samples_names,1,samples_id)
  }else if (9 + length(samples_id) != ncol(variant_df))
    stop(" \n Number of columns do not match number of IDs provided. Was expecting ", ncol(variant_df) - 9, " IDs, but got ", length(samples_id),
         " IDs")



  names(variant_df) = c("Chrom", "Pos", "ID", "REF", "ALT", "Qual", "Filter", "Info", "Format", samples_id)
  rownames(variant_df) = make.unique(paste(variant_df$Chrom, "_", variant_df$Pos, "_", gsub(",", "-", variant_df$REF), "_", gsub(",", "-", variant_df$ALT),
                                           sep = ""))

  if (!is.null(region)) {
    GRegion = parseregion(region)
    if (!is.null(GRegion$Chrom)) {
      variant_df = variant_df[variant_df$Chrom == GRegion$Chrom, ]
      if (!is.null(GRegion$Start) && !is.null(GRegion$End)) {
        startPosition = GRegion$Start
        endPosition = GRegion$End
        variant_df = variant_df[variant_df$Pos >= startPosition & variant_df$Pos <= endPosition, ]
      }
    }
  }




  variant_df
}




# Extract only the variant which passed all filters
#' @export
select_passed <- function(variant_df) {
    variant_df[variant_df$Filter == "PASS", ]
}

# Extract only the variant which has any of the passe dlist of filter
select_filter <- function(variant_df, filter_to_include, filter_to_exclude) {
    # if(length(filter_to_include)>0)
    {
        reg_expression_filter_to_include = paste(filter_to_include, collapse = "|")
        selection_list_to_include = unlist(lapply(unlist(as.character(variant_df$Filter)), function(x) grepl(reg_expression_filter_to_include, x)))

    }
    # if(length(filter_to_exclude)>0)
    {
        reg_expression_filter_to_exclude = paste(filter_to_exclude, collapse = "|")
        selection_list_to_exclude = unlist(lapply(unlist(as.character(variant_df$Filter)), function(x) grepl(reg_expression_filter_to_exclude, x)))
    }

    if (length(filter_to_include) == 0)
        selection_list_to_include = !selection_list_to_exclude
    if (length(filter_to_exclude) == 0)
        selection_list_to_exclude = !selection_list_to_include

    variant_df[selection_list_to_include & !selection_list_to_exclude, ]
}

# Generate a tabular form of the VCF files
#' @export
get_coveragetabular <- function(variant_df) {

    samples_id = colnames(variant_df[10:ncol(variant_df)])

    nb_samples = length(samples_id)

    # Select only well formed info line

    variant_df = variant_df[!is.na(variant_df$Info) & variant_df$Info != "" & grepl("TC=", variant_df$Info), ]

    variant_passed = variant_df
    variant_passed = variant_df[c("Chrom", "Pos", "REF", "ALT", "Qual", "Filter")]

    cat("\nTC...")
    variant_passed["TC"] = unlist(lapply(unlist(variant_df["Info"]), function(x) as.numeric(unlist(strsplit(unlist(strsplit(as.character(x), ";"))[unlist(lapply(unlist(strsplit(as.character(x),
        ";")), function(y) grepl("TC=", y)))], "="))[2])))
    cat("\nTR...")
    variant_passed["TR"] = unlist(lapply(unlist(variant_df["Info"]), function(x) as.numeric(unlist(strsplit(unlist(strsplit(as.character(x), ";"))[unlist(lapply(unlist(strsplit(as.character(x),
        ";")), function(y) grepl("TR=", y)))], "="))[2])))
    cat("\nAF...")
    variant_passed["AF"] = variant_passed["TR"]/variant_passed["TC"]

    cat("\nGT/NR/NV \n Processing sample/well:  : ")

    for (isample in 1:nb_samples) {
        cat(" ", isample)



      values=unlist(lapply(as.character(unlist(variant_df[samples_id[isample]])),
                           function(x) unlist(strsplit(x, ":"))[c(1, 5, 6)]))

      variant_passed[paste(c("GT_", "NR_", "NV_"), samples_id[isample], sep = "")] = matrix(values, ncol=3,byrow=T)



      #  variant_passed[paste(c("GT_", "NR_", "NV_"), samples_id[isample], sep = "")] = unlist(lapply(as.character(unlist(variant_df[samples_id[isample]])),
      #      function(x) unlist(strsplit(x, ":"))[c(1, 5, 6)]))

        # variant_passed[paste('NR_',samples_id[isample],sep='')] = unlist(lapply(as.character(unlist(variant_df[samples_id[isample]])), function(x)
        # unlist(strsplit(x,':'))[5] )) variant_passed[paste('NV_',samples_id[isample],sep='')] =
        # unlist(lapply(as.character(unlist(variant_df[samples_id[isample]])), function(x) unlist(strsplit(x,':'))[6] ))
       # variant_passed[paste("AF_", samples_id[isample], sep = "")] = as.numeric(unlist(variant_passed[paste("NV_", samples_id[isample], sep = "")]))/as.numeric(unlist(variant_passed[paste("NR_",
       #     samples_id[isample], sep = "")]))

        # if(isample==1) variant_passed['GT']=unlist(lapply(as.character(unlist(variant_df[samples_id[isample]])), function(x)
        # unlist(strsplit(x,':'))[1] )) else variant_passed['GT']=paste(unlist(variant_passed['GT']),
        # unlist(lapply(as.character(unlist(variant_df[samples_id[isample]])), function(x) unlist(strsplit(x,':'))[1] )),sep=':')
    }

    cat("\n\n All information retrieved \n\n")
    variant_passed

}


# For 384 wells
#' @export
get_alleleinfotabular <- function(variant_df,RetrieveWellsIDs=FALSE) {

  ##FORMAT=<ID=GT,Number=1,Type=String,Description="Unphased genotypes">
  ##FORMAT=<ID=NR,Number=.,Type=Integer,Description="Number of reads covering variant location in this sample">
  ##FORMAT=<ID=NV,Number=.,Type=Integer,Description="Number of reads containing variant in this sample">
    samples_id = colnames(variant_df[10:ncol(variant_df)])


    variant_passed = get_coveragetabular(variant_df)

    newvariant_df = variant_passed[1:9]

    newvariant_df["WC"] = as.numeric(unlist(apply(variant_passed, 1, function(x) sum(x[paste("NR_", samples_id, sep = "")] > 0 & (x[paste("NR_",
        samples_id, sep = "")] == x[paste("NV_", samples_id, sep = "")] | x[paste("NV_", samples_id, sep = "")] == 0)))))
    newvariant_df["WV"] = as.numeric(unlist(apply(variant_passed, 1, function(x) sum(x[paste("NV_", samples_id, sep = "")] > 0 & (x[paste("NR_",
        samples_id, sep = "")] == x[paste("NV_", samples_id, sep = "")] | x[paste("NV_", samples_id, sep = "")] == 0)))))
    newvariant_df["WR"] = newvariant_df["WC"] -  newvariant_df["WV"]
    wellfraction = as.numeric(unlist(newvariant_df["WV"]/newvariant_df["WC"]))
    newvariant_df["WF"] = wellfraction

    if(RetrieveWellsIDs){

      newvariant_df["WV_IDs"] =apply(variant_passed, 1, function(x)  paste(lapply(names(which( x[paste("NV_", samples_id, sep = "")] > 0 & (x[paste("NR_",  samples_id, sep = "")] == x[paste("NV_", samples_id, sep = "")] | x[paste("NV_", samples_id, sep = "")] == 0))),function(x) substr(x,4,nchar(x))),collapse=":"))

      newvariant_df["WC_IDs"] = apply(variant_passed, 1,        function(x)   paste(lapply(names(which(x[paste("NR_", samples_id, sep = "")] > 0 &                                                 (x[paste("NR_",samples_id, sep = "")] == x[paste("NV_", samples_id, sep = "")] | x[paste("NV_", samples_id, sep = "")] == 0))),function(x) substr(x,4,nchar(x))),collapse=":"))

      newvariant_df["WR_IDs"] = apply(newvariant_df, 1,        function(x)   paste(setdiff(unlist(strsplit(x["WC_IDs"],":")), unlist(strsplit(x["WV_IDs"],":"))),collapse=":"))

      newvariant_df["WC_IDs"] = NULL


    }
  #  newvariant_df["WFAdj"] = unlist(lapply(wellfraction, function(x) if (!is.na(x) && (x == 1))
  #      runif(1, 0.95, 1) else x))

    newvariant_df

}


# plot_allfraction<-function(variant_df,sample, ratio,title=sample) { step = floor(1/ratio) indices=seq(1,nrow(variant_df),step)
# variant_df=variant_df[indices,] #Chromsomes list chrom_list<-paste('chr',c(1:22,'X','Y','M'),sep='') #Chromosomes sizes
# hg19_chrsize<-list(chr1=249250621,chr2=243199373,chr3=198022430, chr4=191154276, chr5=180915260, chr6=171115067, chr7=159138663,
# chr8=146364022, chr9=141213431, chr10=135534747, chr11=135006516, chr12=133851895, chr13=115169878, chr14=107349540,chr15=102531392,
# chr16=90354753, chr17=81195210, chr18=78077248, chr19=59128983, chr20=63025520,chr21=48129895,chr22=51304566, chrX=155270560,
# chrY=59373566,chrM=16571) #Chromosome offset chrom_offset<-c(0) #chrom_middle<-c(0 + hg19_chrsize[chrom_list[1]][[1]]/2) for (chr in
# 2:length(hg19_chrsize)) { chrom_offset<-c(chrom_offset,chrom_offset[chr-1]+as.numeric(hg19_chrsize[chr-1])) # chrom_middle<-c(chrom_middle,
# chrom_offset[chr]+hg19_chrsize[chrom_list[chr]][[1]]/2) } names(chrom_offset) = chrom_list master_df=variant_df #
# master_df=as.numeric(unlist(master_df['Pos'])) + as.numeric(unlist(lapply(as.character(master_df['Chrom'])), function(x) if (!is.na(x) &
# x!='chrX' & x!='chrM') as.numeric(master_dfchrom_offset[x] ) else NA ))) master_df['Pos_offset']=as.numeric(unlist(master_df['Pos'])) +
# as.numeric(unlist(lapply(as.character(unlist(master_df['Chrom'])), function(x) if (!is.na(x) ) as.numeric(chrom_offset[x] ) else NA )))
# cat('\n we plot', sample) plot(unlist(master_df['Pos_offset']),unlist(master_df[paste('AF_',sample,sep='')]), cex=0.2,main=paste('Allele
# Fraction plot : ',title),ylab='Allelefraction',xlab='n',xaxt='n',cex.main=1, col='black' ) v0=as.numeric(chrom_offset); l=length(v0); v0[l+1]=
# master_df$Pos_offset[nrow(master_df)];l=length(v0); v1=(v0[1:(l-1)] + v0[2:l])/2 abline(v =v0 ,lwd=2,col = 'black');
# mtext(chrom_list,side=1,line=1, at=v1,las=2) } plot_allfraction_v2<-function(variant_df, ratio, chrom) { step = floor(1/ratio)
# indices=seq(1,nrow(variant_df),step) variant_df=variant_df[indices,] #Chromsomes list chrom_list<-paste('chr',c(1:22,'X','Y','M'),sep='')
# #Chromosomes sizes hg19_chrsize<-list(chr1=249250621,chr2=243199373,chr3=198022430, chr4=191154276, chr5=180915260, chr6=171115067,
# chr7=159138663, chr8=146364022, chr9=141213431, chr10=135534747, chr11=135006516, chr12=133851895, chr13=115169878,
# chr14=107349540,chr15=102531392, chr16=90354753, chr17=81195210, chr18=78077248, chr19=59128983, chr20=63025520,chr21=48129895,chr22=51304566,
# chrX=155270560, chrY=59373566,chrM=16571) #Chromosome offset chrom_offset<-c(0) #chrom_middle<-c(0 + hg19_chrsize[chrom_list[1]][[1]]/2) for
# (chr in 2:length(hg19_chrsize)) { chrom_offset<-c(chrom_offset,chrom_offset[chr-1]+as.numeric(hg19_chrsize[chr-1])) #
# chrom_middle<-c(chrom_middle, chrom_offset[chr]+hg19_chrsize[chrom_list[chr]][[1]]/2) } names(chrom_offset) = chrom_list master_df=variant_df
# # master_df=as.numeric(unlist(master_df['Pos'])) + as.numeric(unlist(lapply(as.character(master_df['Chrom'])), function(x) if (!is.na(x) &
# x!='chrX' & x!='chrM') as.numeric(master_dfchrom_offset[x] ) else NA ))) master_df['Pos_offset']=as.numeric(unlist(master_df['Pos'])) +
# as.numeric(unlist(lapply(as.character(unlist(master_df['Chrom'])), function(x) if (!is.na(x) ) as.numeric(chrom_offset[x] ) else NA )))
# par(mfrow=c(5,1)) cat('\n we plot the well fraction Adjusted') plot(unlist(master_df['Pos_offset']),unlist(master_df['WFAdj']),
# cex=0.2,main=paste(chrom, 'Well Fraction Adjusted'),ylab=' Well fraction',cex.main=2,cex.label=1.5, col='black',xlab=NA ) cat('\n we plot the
# well fraction') plot(unlist(master_df['Pos_offset']),unlist(master_df['WF']), cex=0.2,main=paste(chrom, 'Well Fraction'),ylab=' Well
# fraction',cex.main=2,cex.label=1.5, col='black',xlab=NA ) # v0=as.numeric(chrom_offset); l=length(v0); v0[l+1]=
# master_df$Pos_offset[nrow(master_df)];l=length(v0); # v1=(v0[1:(l-1)] + v0[2:l])/2 # abline(v =v0 ,lwd=2,col = 'black');
# mtext(chrom_list,side=1,line=1, at=v1,las=2) cat('\n we plot the well count') plot(unlist(master_df['Pos_offset']),unlist(master_df['WR']),
# cex=0.2,main=paste(chrom, ' Well count'),ylab='Total Well Count',cex.main=2,cex.label=1.5, col='black',xlab=NA ) #
# v0=as.numeric(chrom_offset); l=length(v0); v0[l+1]= master_df$Pos_offset[nrow(master_df)];l=length(v0); #v1=(v0[1:(l-1)] + v0[2:l])/2 #
# abline(v =v0 ,lwd=2,col = 'black'); mtext(chrom_list,side=1,line=1, at=v1,las=2) cat('\n we plot the Read fraction')
# plot(unlist(master_df['Pos_offset']),unlist(master_df['AF']), cex=0.2,main=paste(chrom, 'Read Fraction '),ylab='Reads
# fraction',cex.main=2,cex.label=1.5, col='black',xlab=NA ) # v0=as.numeric(chrom_offset); l=length(v0); v0[l+1]=
# master_df$Pos_offset[nrow(master_df)];l=length(v0); # v1=(v0[1:(l-1)] + v0[2:l])/2 # abline(v =v0 ,lwd=2,col = 'black');
# mtext(chrom_list,side=1,line=1, at=v1,las=2) cat('\n we plot the Read count') plot(unlist(master_df['Pos_offset']),unlist(master_df['TC']),
# cex=0.2,main=paste(chrom, ' Reads count'),ylab='Total Reads Count',cex.main=2,cex.label=1.5, col='black',xlab=NA ) #
# v0=as.numeric(chrom_offset); l=length(v0); v0[l+1]= master_df$Pos_offset[nrow(master_df)];l=length(v0); # v1=(v0[1:(l-1)] + v0[2:l])/2 #
# abline(v =v0 ,lwd=2,col = 'black'); mtext(chrom_list,side=1,line=1, at=v1,las=2) } plot_allfraction_chrom<-function(variant_df, ratio, chrom)
# { variant_df =variant_df[variant_df$Chrom==chrom,] step = floor(1/ratio) indices=seq(1,nrow(variant_df),step) variant_df=variant_df[indices,]
# #Chromsomes list chrom_list<-paste('chr',c(1:22,'X','Y','M'),sep='') #Chromosomes sizes
# hg19_chrsize<-list(chr1=249250621,chr2=243199373,chr3=198022430, chr4=191154276, chr5=180915260, chr6=171115067, chr7=159138663,
# chr8=146364022, chr9=141213431, chr10=135534747, chr11=135006516, chr12=133851895, chr13=115169878, chr14=107349540,chr15=102531392,
# chr16=90354753, chr17=81195210, chr18=78077248, chr19=59128983, chr20=63025520,chr21=48129895,chr22=51304566, chrX=155270560,
# chrY=59373566,chrM=16571) #Chromosome offset chrom_offset<-c(0) #chrom_middle<-c(0 + hg19_chrsize[chrom_list[1]][[1]]/2) for (chr in
# 2:length(hg19_chrsize)) { chrom_offset<-c(chrom_offset,chrom_offset[chr-1]+as.numeric(hg19_chrsize[chr-1])) # chrom_middle<-c(chrom_middle,
# chrom_offset[chr]+hg19_chrsize[chrom_list[chr]][[1]]/2) } names(chrom_offset) = chrom_list master_df=variant_df par(mfrow=c(5,1)) cat('\n we
# plot the well fraction') plot(unlist(master_df['Pos']),unlist(master_df['WF']), cex=0.2,main=paste(chrom, 'Well Fraction'),ylab=' Well
# fraction',cex.main=1.5,cex.lab=1.5,cex.axis=1.5, col='black',xlab=NA,xlim=c(0, as.numeric(hg19_chrsize[chrom])) ) cat('\n we plot the well
# fraction Adjusted') plot(unlist(master_df['Pos']),unlist(master_df['WFAdj']), cex=0.2,main=paste(chrom, 'Well Fraction Adjusted'),ylab=' Well
# fraction',cex.main=1.5,cex.lab=1.5,cex.axis=1.5, col='black',xlab=NA,xlim=c(0, as.numeric(hg19_chrsize[chrom])) ) #
# v0=as.numeric(chrom_offset); l=length(v0); v0[l+1]= master_df$Pos_offset[nrow(master_df)];l=length(v0); # v1=(v0[1:(l-1)] + v0[2:l])/2 #
# abline(v =v0 ,lwd=2,col = 'black'); mtext(chrom_list,side=1,line=1, at=v1,las=2) cat('\n we plot the well count')
# plot(unlist(master_df['Pos']),unlist(master_df['WR']), cex=0.2,main=paste(chrom, ' Well count'),ylab='Total Well
# Count',cex.main=1.5,cex.lab=1.5,cex.axis=1.5, col='black',xlab=NA,xlim=c(0, as.numeric(hg19_chrsize[chrom]) )) # v0=as.numeric(chrom_offset);
# l=length(v0); v0[l+1]= master_df$Pos_offset[nrow(master_df)];l=length(v0); #v1=(v0[1:(l-1)] + v0[2:l])/2 # abline(v =v0 ,lwd=2,col = 'black');
# mtext(chrom_list,side=1,line=1, at=v1,las=2) cat('\n we plot the Read fraction') plot(unlist(master_df['Pos']),unlist(master_df['AF']),
# cex=0.2,main=paste(chrom, 'Read Fraction '),ylab='Reads fraction',cex.main=1.5,cex.lab=1.5,cex.axis=1.5, col='black',xlab=NA,xlim=c(0,
# as.numeric(hg19_chrsize[chrom]) )) # v0=as.numeric(chrom_offset); l=length(v0); v0[l+1]= master_df$Pos_offset[nrow(master_df)];l=length(v0); #
# v1=(v0[1:(l-1)] + v0[2:l])/2 # abline(v =v0 ,lwd=2,col = 'black'); mtext(chrom_list,side=1,line=1, at=v1,las=2) cat('\n we plot the Read
# count') plot(unlist(master_df['Pos']),unlist(master_df['TC']), cex=0.2,main=paste(chrom, ' Reads count'),ylab='Total Reads
# Count',cex.main=1.5,cex.lab=1.5,cex.axis=1.5, col='black',xlab=NA,xlim=c(0, as.numeric(hg19_chrsize[chrom]) )) # v0=as.numeric(chrom_offset);
# l=length(v0); v0[l+1]= master_df$Pos_offset[nrow(master_df)];l=length(v0); # v1=(v0[1:(l-1)] + v0[2:l])/2 # abline(v =v0 ,lwd=2,col =
# 'black'); mtext(chrom_list,side=1,line=1, at=v1,las=2) }
