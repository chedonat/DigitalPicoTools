
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
getMutationsOfLFR2<-function(LFR_withMutations_df,LFRname, list_of_mutations_to_phase, calltype=NULL){

  #cat("\n\n",LFRname,"\n\n")
  if(calltype=="SNP")
    mutations_list=LFR_withMutations_df[LFRname,"MutationsOnSNP"]
  if(calltype=="REF")
    mutations_list=LFR_withMutations_df[LFRname,"MutationsOnREF"]
  mutations=""


  if(mutations_list !="")
    mutations=unlist(strsplit(mutations_list,":"))
  #We need to consider only  mutations which belong to list_of_mutations_to_phase
  mutations=intersect(mutations,list_of_mutations_to_phase)

  mutations

}


#' @export
AssignPhasingCode<-function(Mutation_list, phasing_code, SNPflag#,MutationPhasingCode_df
)
{
  MutationPhasingCode_df[Mutation_list, "PhasingCode1"] <<- paste(phasing_code,SNPflag,sep="_")
}

#' @export
getOverlapLFR<-function(LFR_withMutations_df,LFRname)
{
  #rownames(LFR_withMutations_df) = LFR_withMutations_df$LFR_name
  LFRname_elements=unlist(strsplit(LFRname,"_"))
  #   startpos = LFR_withMutations_df[LFRname,"Start"]
  #   endpos = LFR_withMutations_df[LFRname,"End"]
  #   chrom = as.character(LFR_withMutations_df[LFRname,"Chrom"])
  #   wellID= as.character(LFR_withMutations_df[LFRname,"Well_ID"])
  startpos = as.numeric(LFRname_elements[3])
  endpos = as.numeric(LFRname_elements[4])
  chrom = LFRname_elements[2]
  wellID= LFRname_elements[1]

  Overlap_LFR_lst = LFR_withMutations_df[(LFR_withMutations_df$Start >= startpos  & LFR_withMutations_df$Start <= endpos) &
                             (LFR_withMutations_df$End >= startpos  & LFR_withMutations_df$End <= endpos) &
                             LFR_withMutations_df$LFR_name != LFRname ,]

  as.character(Overlap_LFR_lst$LFR_name)
}

#' @export
getOverlapMutations<-function(markedLFR,LFR_withMutations_df,LFRname1, LFRname2,list_of_mutations_to_phase)
{

  mutation_list1_snp= getMutationsOfLFR2(LFR_withMutations_df,LFRname1,list_of_mutations_to_phase,"SNP")
  mutation_list2_snp= getMutationsOfLFR2(LFR_withMutations_df,LFRname2,list_of_mutations_to_phase,"SNP")
  mutation_list1_ref= getMutationsOfLFR2(LFR_withMutations_df,LFRname1,list_of_mutations_to_phase,"REF")
  mutation_list2_ref= getMutationsOfLFR2(LFR_withMutations_df,LFRname2,list_of_mutations_to_phase,"REF")

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
process_fragments<-function(LFR_withMutations_df,LFRname, phasingCode, SNPflag,markedLFR,list_of_mutations_to_phase#,MutationPhasingCode_df
)
{
  # cat("\n lfr is ", LFRname)

  if(markedLFR[LFRname]==FALSE){

    markedLFR[LFRname]<-TRUE

    # cat("\n markedLFR ",LFRname, " value:",markedLFR[LFRname] )

    mutations_snp=getMutationsOfLFR2(LFR_withMutations_df,LFRname, list_of_mutations_to_phase,"SNP")
    mutations_ref=getMutationsOfLFR2(LFR_withMutations_df,LFRname,list_of_mutations_to_phase, "REF")

    if(length(c(mutations_snp,mutations_ref))!=0){
      AssignPhasingCode(mutations_snp,phasingCode, SNPflag#,MutationPhasingCode_df
      )
      AssignPhasingCode(mutations_ref,phasingCode, (1-SNPflag)#, MutationPhasingCode_df
      )

      #Get the list of Overlaping LFR
      overlapLFR_lst<-getOverlapLFR(LFR_withMutations_df,LFRname)
      # cat("\n number overlap ", length(overlapLFR_lst))
      #create a submutation of LFR and the overlaping LFR
      SubLFR_withMutations_df =LFR_withMutations_df[unique(c(LFRname,overlapLFR_lst)),]

      #For each non processed LFR here goes the recurrrence
      for (lfr in overlapLFR_lst){
        if(markedLFR[lfr]){
          # cat(" Ttreated")
          next
        }

        #Get the list of Overlaping Mutations
        overlapmutations<-getOverlapMutations(markedLFR,SubLFR_withMutations_df,LFRname, lfr,list_of_mutations_to_phase)
        lengthOverlap=as.numeric(unlist(lapply(overlapmutations,length)))
        names(lengthOverlap) = names(overlapmutations)


        if(max(lengthOverlap)==0)
          next
        check=F
        if(check){
          cat("\n")
          print(lengthOverlap)
        }


        #   cat("\n max", max(lengthOverlap),"\n")

        #   print(lengthOverlap)

        maxoverlap=names(lengthOverlap[which.max(lengthOverlap)])

        notnulloverlap=names(lengthOverlap[lengthOverlap>0])
        #If both same allele and different allele present skip
        if(length(intersect(c("overlap11","overlap00"),notnulloverlap ))>0 && length(intersect(c("overlap01","overlap10"),notnulloverlap ))>0 ){
          next
        }


        #  cat("\n treating recursively  ", lfr )

        if(maxoverlap %in% c("overlap11","overlap00")){
          markedLFR =  process_fragments(LFR_withMutations_df,lfr, phasingCode,SNPflag,markedLFR, list_of_mutations_to_phase#,MutationPhasingCode_df
          )
        }else if (maxoverlap %in% c("overlap10","overlap01")){
          markedLFR =   process_fragments(LFR_withMutations_df,lfr, phasingCode,(1-SNPflag),markedLFR, list_of_mutations_to_phase#,MutationPhasingCode_df
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










