# fucntion to return age 
  age_calculation <- function(data){
    data$onset_year <- as.integer(str_sub(data$epistart,1,4))
    data$birth_year <- unlist(cd1[which(cd1$eid %in% data$eid), 'year_birth'])
    data$onset_age <- data$onset_year - data$birth_year
    return(data)
  }

# read the baseline self-report information
  idex_case_report_base <- function(selfreport1,selfreport2,disease_age1,disease_age2,code1,code2=NULL){
    index <- which(unlist(selfreport1[,2]) %in% code1)
    age <- disease_age1[,2]
    age[-index] <- NA
    for (i in 3:nrow(selfreport1)){
      index_new <- which(unlist(selfreport1[,..i]) %in% code1)
      if (length(index_new)!=0){
        age_new <- disease_age1[,..i]
        age_new[-index_new] <- NA
        age <- apply(cbind(age,age_new),1,my.min)
        index <- c(index,index_new)}
    }
    if (!is.null(code2)){
      for (i in 2:nrow(selfreport2)){
        index_new <- which(unlist(selfreport2[,..i]) %in% code2)
        if (length(index_new)!=0){
          age_new <- disease_age2[,..i]
          age_new[-index_new] <- NA
          age <- apply(cbind(age,age_new),1,my.min)
          index <- c(index, index_new)}
      }
    }
    res <- list(index,age)
    return(res)
  }

  idex_case_report <- function(selfreport1,selfreport2,disease_age1,disease_age2,code1,code2=NULL){
    index <- which(unlist(selfreport1[,2]) %in% code1)
    age <- disease_age1[,2]
    age[-index] <- NA
    for (i in 3:ncol(selfreport1)){
      index_new <- which(unlist(selfreport1[,..i]) %in% code1)
      if (length(index_new)!=0){
        age_new <- disease_age1[,..i]
        age_new[-index_new] <- NA
        age <- apply(cbind(age,age_new),1,my.min)
        index <- c(index,index_new)}
    }
    if (!is.null(code2)){
      for (i in 2:ncol(selfreport2)){
        index_new <- which(unlist(selfreport2[,..i]) %in% code2)
        if (length(index_new)!=0){
          age_new <- disease_age2[,..i]
          age_new[-index_new] <- NA
          age <- apply(cbind(age,age_new),1,my.min)
          index <- c(index, index_new)}
      }
    }
    res <- list(index,age)
    return(res)
  }
  
  idex_cancer_report_base <- function(cancer_report, cancer_age, code1){
    index <- which(unlist(cancer_report[,2]) %in% code1)
    age <- cancer_age[,2]
    age[-index] <- NA
    #for (i in 3:19){
    for (i in 3:ncol(cancer_report)){
      index_new <- which(unlist(cancer_report[,..i]) %in% code1)
      if (length(index_new)!=0){
      age_new <- cancer_age[,..i]
      age_new[-index_new] <- NA
      age <- apply(cbind(age,age_new),1,my.min)
      index <- c(index,index_new)}
    }
    res <- list(index,age)
    return(res)
  }

  idex_cancer_report <- function(cancer_report, cancer_age, code1){
    index <- which(unlist(cancer_report[,2]) %in% code1)
    age <- cancer_age[,2]
    age[-index] <- NA
    for (i in 3:ncol(cancer_report)){
      index_new <- which(unlist(cancer_report[,..i]) %in% code1)
      if (length(index_new)!=0){
      age_new <- cancer_age[,..i]
      age_new[-index_new] <- NA
      age <- apply(cbind(age,age_new),1,my.min)
      index <- c(index,index_new)}
    }
    res <- list(index,age)
    return(res)
  }


# function to select the ealiest time
  my.min <- function(x) ifelse( !all(is.na(x)), min(x, na.rm=T), NA)
