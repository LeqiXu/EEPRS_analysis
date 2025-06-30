ukbcase <- function(icd10=NULL,icd9=NULL,oper4=NULL,histology=NULL,
                    hesin=NULL,hesin_diag10=NULL,hesin_diag9=NULL,hesin_oper4=NULL,
                    deathdata=NULL,cancerdata=NULL,
                    icd10main=T,icd10sec=T,icd9main=T,icd9sec=T,oper4main=T,oper4sec=T,
                    death=T,cancer=T) {
  
  case <- NULL # initialization
  if(is.null(icd10) & is.null(icd9) & is.null(oper4) & is.null(histology)) {stop("check the definition of disease")}
  
  if(((icd10main | icd10sec | icd9main | icd9sec) & (!is.null(icd10))) == T){
    if(is.null(hesin)){stop('main hesin dataset is required')}
    case_diag <- ukbcase_diag(icd10,icd9,hesin,hesin_diag10,hesin_diag9,icd10main,icd10sec,icd9main,icd9sec)
    case <- case_diag
  }
  
  if(((oper4main | oper4sec) & (!is.null(oper4))) == T){
    if(is.null(hesin)){stop('main hesin dataset is required')}
    case_oper <- ukbcase_oper(oper4,hesin,hesin_oper4,oper4main,oper4sec)
    if(is.null(case)){case <- case_oper}
    else{case <- merge(case_diag,case_oper,all=T,by=c("eid","ins_index","epistart","epiend"))}
  }
  
  if(death==T & (!is.null(deathdata))){
    deathcase <- deathdata[which(deathdata$cause_of_death %in% icd10),]
    deathcase$death_description <- 'death'
    deathcase$date_of_death <- as.Date(deathcase$date_of_death,format='%Y-%m-%d')
    if(is.null(case)){
      case <- deathcase
    }else{case <- merge(case, deathcase, all=T,by.x=c("eid","diag_icd10","epistart"),
                          by.y=c("eid","cause_of_death","date_of_death"))}
  }
  
  if(cancer==T & (!is.null(cancerdata))){
    if(is.null(icd10) & is.null(icd9) & is.null(histology)) {stop("check the definition of disease")}
    cancercase <- cancerdata[unique(c(which(cancerdata$icd10 %in% icd10),which(cancerdata$icd9 %in% icd9),
                                      which(cancerdata$histology %in% histology))),c('eid','date','icd10','icd9','histology')]
    cancercase$date <- as.Date(cancercase$date)
    names(cancercase) <- c('eid','epistart','diag_icd10','diag_icd9','cancer_histology')
    if(is.null(case)){case <- cancercase}
    else{case <- merge(case, cancercase, all=T,by=c("eid","diag_icd10","diag_icd9","epistart"))}
  }
  
  # get unique case with the earlist onset date
  case_noepi=case[is.na(case$epistart),]
  case_noepi=case_noepi[!duplicated(case_noepi$eid),]
  #if( nrow( case_noepi)!=0){
  #  stop("Check if there are some update of the format of date.")
  #}
  case <- case[!is.na(case$epistart),]
  colid<-sapply(unique(case$eid),function(x){
    col<-which(case$eid==x & case$epistart==min(case$epistart[which(case$eid==x)],na.rm=T))[1]
    return(col)
  })
  case <- case[colid,]
  case <- case[!is.na(case$eid),]
  case_noepi=case_noepi[!case_noepi$eid %in% case$eid,]
  case=rbind(case, case_noepi)
  if('epiend' %in% colnames(case)){case$epiend <- as.Date(case$epiend)}
  return(case)
}



###################### diag
ukbcase_diag <- function(icd10=NULL,icd9=NULL,hesin=NULL,hesin_diag10=NULL,hesin_diag9=NULL,
                         icd10main=T,icd10sec=T,icd9main=T,icd9sec=T) {
  hesin_diag <- NULL # initialization
  if(icd10main==F & icd10sec==F & icd9main==F) {stop('at least one icd code should be given when using function ukbcase_diag')}

  # all icd10 from hesin_diag10 subdataset
  if(icd10sec==T & (!is.null(hesin_diag10))) {
    hesin_diag_2 <- hesin_diag10[which(hesin_diag10$diag_icd10 %in% icd10),c('eid','ins_index','diag_icd10')]
    # find the date of each record from the main hesin dateset
    hesin_diag_2 <- merge(hesin_diag_2,hesin, by=c('eid','ins_index'),all=F)
    hesin_diag_2$epistart <- as.Date(hesin_diag_2$epistart)

    hesin_diag <- hesin_diag_2
  }

  # all icd9 from hesin_diag9 dataset
  if(icd9sec==T & (!is.null(hesin_diag9))) {
    hesin_diag_3 <- hesin_diag9[which(hesin_diag9$diag_icd9 %in% icd9),c('eid','ins_index','diag_icd9')]
    # find the date of each record from the main hesin dateset
    hesin_diag_3 <- merge(hesin_diag_3,hesin, by=c('eid','ins_index'),all=F)
    hesin_diag_3$epistart <- as.Date(hesin_diag_3$epistart)

    if(is.null(hesin_diag)){
      hesin_diag <- hesin_diag_3
    }else{
      hesin_diag <- merge(hesin_diag,hesin_diag_3,all=T,by = c("eid","ins_index","epistart","epiend"))
    }
  }

  # select the epi with earlist date
  hesin_diag <- hesin_diag[which(!is.na(hesin_diag$eid)),]
  colid<-sapply(unique(hesin_diag$eid),function(x){
    col<-which(hesin_diag$eid==x & hesin_diag$epistart==min(hesin_diag$epistart[which(hesin_diag$eid==x)]))[1]
    return(col)
  })
  hesin_diag <- hesin_diag[colid,]
  hesin_diag <- hesin_diag[which(!is.na(hesin_diag$eid)),]

  return(hesin_diag)
}

###################### oper
ukbcase_oper <- function(oper4=NULL,hesin=NULL,hesin_oper4=NULL,
                         oper4main=T,oper4sec=T) {
  hesin_oper <- NULL # initialization
  if(oper4main==F & oper4sec==F) {stop('at least one oper4 code should be given when using function ukbcase_oper')}

  # secondary oper4 from hesin_oper4 subdataset
  if(oper4sec==T & (!is.null(hesin_oper4))) {
    oper_col <- which(hesin_oper4$oper4 %in% oper4)
    hesin_oper_2 <- hesin_oper4[oper_col,c('eid','ins_index','oper4')]
# # 6007659         6  K403     SMR     10     <NA>   <NA>
    # find the date of each record from the main hesin dateset
    hesin_oper_2 <- merge(hesin_oper_2,hesin, by=c('eid','ins_index'),all=F)
    hesin_oper_2$epistart <- as.Date(hesin_oper_2$epistart)

    hesin_oper <- hesin_oper_2
  }

  hesin_oper <- hesin_oper[which(!is.na(hesin_oper$eid)),]

  return(hesin_oper)
}
