
#' Retrieve the  mutations on an LFR
#'
#'This function computes the list of mutations present on a  Long Fragment Reads either on their Variant Allele or their reference Allele
#'
#' @export
getMutationsOfLFR<-function(LFR_info, Mutations_df, calltype=NULL)
{
  if(is.null(LFR_info))
    if(length(LFR_info) ==4) {
      names(LFR_info)=c("LFR_name","Chrom","Start","End","Well_ID")
    }else{
      cat("\n LFR_info should contains the following", c("LFR_name ","Chrom ","Start ","End ","Well_ID "))
    }

  startpos = as.numeric(unlist(LFR_info["Start"]))
  endpos = as.numeric(unlist(LFR_info["End"]))
  chrom = as.character(unlist(LFR_info["Chrom"]))
  wellID= as.character(unlist(LFR_info["Well_ID"]))

  subMutations_df=Mutations_df[Mutations_df$Pos>=startpos  & Mutations_df$Pos<=endpos & Mutations_df$Chrom==chrom,]

  if(is.null(calltype))
  {
    mutations_list=subMutations_df[!is.na(subMutations_df$WV_IDs) & (grepl(wellID,subMutations_df$WV_IDs) |grepl(wellID,subMutations_df$WR_IDs)),  ]
  }else if(calltype=="SNP"){
    mutations_list=subMutations_df[!is.na(subMutations_df$WV_IDs) & grepl(wellID,subMutations_df$WV_IDs),  ]
    #mutations_list=subMutations_df[!is.na(subMutations_df$WV_IDs) & wellID %in%subMutations_df$WV_IDs),  ]
  }else if(calltype=="REF"){
    mutations_list=subMutations_df[grepl(wellID,subMutations_df$WR_IDs),  ]
  }else{
    stop("\n\n Calltype shuld be either SNP either REF \n\n ")
  }

  rownames(mutations_list)


}




#After the computation done with getMutationsOfLFR, simply retrieve the list of mutations in the LFR under REF or under variant.
#' @export
getMutationsOfLFR2<-function(MutationsLFR_df,LFRname, calltype=NULL){

  #cat("\n\n",LFRname,"\n\n")
  if(calltype=="SNP")
    mutations_list=MutationsLFR_df[LFRname,"MutationsOnSNP"]
  if(calltype=="REF")
    mutations_list=MutationsLFR_df[LFRname,"MutationsOnREF"]
  mutations=""


  if(mutations_list !="")
    mutations=unlist(strsplit(mutations_list,":"))
  mutations

}


#' @export
AssignPhasingCode<-function(Mutation_list, phasing_code, SNPflag#,MutationPhasingCode_df
)
{
  MutationPhasingCode_df[Mutation_list, "PhasingCode1"] <<- paste(phasing_code,SNPflag,sep="_")
}

#' @export
getOverlapLFR<-function(MutationsLFR_df,LFRname)
{
  #rownames(MutationsLFR_df) = MutationsLFR_df$LFR_name
  LFRname_elements=unlist(strsplit(LFRname,"_"))
  #   startpos = MutationsLFR_df[LFRname,"Start"]
  #   endpos = MutationsLFR_df[LFRname,"End"]
  #   chrom = as.character(MutationsLFR_df[LFRname,"Chrom"])
  #   wellID= as.character(MutationsLFR_df[LFRname,"Well_ID"])
  startpos = as.numeric(LFRname_elements[3])
  endpos = as.numeric(LFRname_elements[4])
  chrom = LFRname_elements[2]
  wellID= LFRname_elements[1]

  Overlap_LFR_lst = MutationsLFR_df[(MutationsLFR_df$Start >= startpos  & MutationsLFR_df$Start <= endpos) &
                             (MutationsLFR_df$End >= startpos  & MutationsLFR_df$End <= endpos) &
                             MutationsLFR_df$LFR_name != LFRname ,]

  as.character(Overlap_LFR_lst$LFR_name)
}

#' @export
getOverlapMutations<-function(markedLFR,MutationsLFR_df,LFRname1, LFRname2)
{

  mutation_list1_snp= getMutationsOfLFR2(MutationsLFR_df,LFRname1,"SNP")
  mutation_list2_snp= getMutationsOfLFR2(MutationsLFR_df,LFRname2,"SNP")
  mutation_list1_ref= getMutationsOfLFR2(MutationsLFR_df,LFRname1,"REF")
  mutation_list2_ref= getMutationsOfLFR2(MutationsLFR_df,LFRname2,"REF")

  if(length(c(mutation_list2_snp,mutation_list2_ref)==0))
   # markedLFR[LFRname2]<<-TRUE
    eval.parent(substitute(markedLFR[LFRname2]<-TRUE))

  # cat("here")
  #cat("\n LFR1 ", getMutationsOfLFR(LFRname1))
  # cat("\n LFR2", getMutationsOfLFR(LFRname2))

  overlap11=intersect(mutation_list1_snp,mutation_list2_snp)
  overlap10=intersect(mutation_list1_snp,mutation_list2_ref)
  overlap00=intersect(mutation_list1_ref,mutation_list2_ref)
  overlap01=intersect(mutation_list1_ref,mutation_list2_snp)

  list(overlap11=overlap11,overlap10=overlap10,overlap00=overlap00, overlap01=overlap01)
}

#' @export
process_fragments<-function(MutationsLFR_df,LFRname, phasingCode, SNPflag,markedLFR#,MutationPhasingCode_df
)
{
  # cat("\n lfr is ", LFRname)

  if(markedLFR[LFRname]==FALSE){

    markedLFR[LFRname]<-TRUE

    # cat("\n markedLFR ",LFRname, " value:",markedLFR[LFRname] )

    mutations_snp=getMutationsOfLFR2(MutationsLFR_df,LFRname, "SNP")
    mutations_ref=getMutationsOfLFR2(MutationsLFR_df,LFRname, "REF")

    if(length(c(mutations_snp,mutations_ref))!=0){
      AssignPhasingCode(mutations_snp,phasingCode, SNPflag#,MutationPhasingCode_df
      )
      AssignPhasingCode(mutations_ref,phasingCode, (1-SNPflag)#, MutationPhasingCode_df
      )

      overlapLFR_lst<-getOverlapLFR(MutationsLFR_df,LFRname)
      # cat("\n number overlap ", length(overlapLFR_lst))

      SubMutationsLFR_df =MutationsLFR_df[unique(c(LFRname,overlapLFR_lst)),]
      for (lfr in overlapLFR_lst){
        if(markedLFR[lfr]){
          # cat(" Ttreated")
          next
        }

        overlapmutations<-getOverlapMutations(markedLFR,SubMutationsLFR_df,LFRname, lfr)
        lengthverlap=as.numeric(unlist(lapply(overlapmutations,length)))
        names(lengthverlap) = names(overlapmutations)


        if(max(lengthverlap)==0)
          next
        check=F
        if(check){
          cat("\n")
          print(lengthverlap)
        }


        #   cat("\n max", max(lengthverlap),"\n")

        #   print(lengthverlap)

        maxoverlap=names(lengthverlap[which.max(lengthverlap)])

        #  cat("\n treating recursively  ", lfr )

        if(maxoverlap %in% c("overlap11","overlap00")){
          markedLFR =  process_fragments(MutationsLFR_df,lfr, phasingCode,SNPflag,markedLFR#,MutationPhasingCode_df
          )
        }else if (maxoverlap %in% c("overlap10","overlap01")){
          markedLFR =   process_fragments(MutationsLFR_df,lfr, phasingCode,(1-SNPflag),markedLFR#,MutationPhasingCode_df
          )
        }else{
          cat(" \n unknown maxoverlap : ",maxoverlap )
        }
      }

    }else{
      # cat("\n No mutations called in  the Fragment")

    }
  }else
  {
    # cat("\n>> markedLFR ",LFRname, " value:",markedLFR[LFRname] )
    #cat(" treated")
  }


  markedLFR

}










