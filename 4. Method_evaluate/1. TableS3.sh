# PRS evaluation
library(data.table)
library(stringr)
library(dplyr)
library(pROC)
library(writexl)

cov_choice = c("age_recruit","sex",paste0("PC",1:20))
total_covariates = fread("/gpfs/gibbs/pi/zhao/lx94/JointPRS/data/ukbb_data/cov_data/adjust.csv")
total_covariates = total_covariates %>% select(all_of(c("eid",cov_choice)))

pop = "EUR"

trait_list = c("HDL","LDL","TC","logTG","BMI","WHR","AgeSmk","SmkI","SmkC","CigDay","DrnkWk","Glu2h","HbA1c","eGFR","logCp",
                "AD","Angina","Asthma","AF","ADHD","ASD","BIP","BrC","CKD","CAD","HF","IBD","IS","LuC",
                "MDD","Osteoporosis","OvC","PAD","PBC","PSC","PrC","PAH","RA","SCZ","SLE","T2D")

continuous_trait = c("HDL","LDL","TC","logTG","BMI","WHR","AgeSmk","CigDay","DrnkWk","Glu2h","HbA1c","eGFR","logCp")

prs_final_table = data.table()

for (trait in trait_list){
## test phenotype
Trait_pheno <- fread(paste0("/gpfs/gibbs/pi/zhao/lx94/EEPRS/data/ukbb_pheno_data/ukbb_",trait,".tsv"))
Trait_pheno_name = c("eid",trait)
Trait_pheno = Trait_pheno[,..Trait_pheno_name]
colnames(Trait_pheno) = c("eid","pheno")
Trait_pheno = na.omit(Trait_pheno)

## depend on the trait type and pop type, we read the score accordingly
PRScs = fread(paste0("/gpfs/gibbs/pi/zhao/lx94/EEPRS/result/Trait_PRS/UKB_test_1over3_",trait,"_PRScsx_prs_",pop,".sscore"))
EEPRS_Word2Vec100 = fread(paste0("/gpfs/gibbs/pi/zhao/lx94/EEPRS/result/MIXPRS/UKB_test_1over3_",trait,"_EEPRS_word2vec100_prs_",pop,".sscore"))
EEPRS_Word2Vec100_PCA30 = fread(paste0("/gpfs/gibbs/pi/zhao/lx94/EEPRS/result/MIXPRS/UKB_test_1over3_",trait,"_EEPRS_word2vec100_PCA30_prs_",pop,".sscore"))
EEPRS_Word2Vec100_ICA30 = fread(paste0("/gpfs/gibbs/pi/zhao/lx94/EEPRS/result/MIXPRS/UKB_test_1over3_",trait,"_EEPRS_word2vec100_ICA30_prs_",pop,".sscore"))
EEPRS_gpticd3072_PCA67 = fread(paste0("/gpfs/gibbs/pi/zhao/lx94/EEPRS/result/MIXPRS/UKB_test_1over3_",trait,"_EEPRS_gpticd3072_PCA67_prs_",pop,".sscore"))
EEPRS_gpticd3072_ICA67 = fread(paste0("/gpfs/gibbs/pi/zhao/lx94/EEPRS/result/MIXPRS/UKB_test_1over3_",trait,"_EEPRS_gpticd3072_ICA67_prs_",pop,".sscore"))

## match the IID for pheno and PRS
Trait_pheno = Trait_pheno[which(Trait_pheno$eid %in% PRScs$IID),]
Trait_pheno_id = Trait_pheno$eid

PRScs = PRScs[match(Trait_pheno_id,PRScs$IID),]
EEPRS_Word2Vec100 = EEPRS_Word2Vec100[match(Trait_pheno_id,EEPRS_Word2Vec100$IID),]
EEPRS_Word2Vec100_PCA30 = EEPRS_Word2Vec100_PCA30[match(Trait_pheno_id,EEPRS_Word2Vec100_PCA30$IID),]
EEPRS_Word2Vec100_ICA30 = EEPRS_Word2Vec100_ICA30[match(Trait_pheno_id,EEPRS_Word2Vec100_ICA30$IID),]
EEPRS_gpticd3072_PCA67 = EEPRS_gpticd3072_PCA67[match(Trait_pheno_id,EEPRS_gpticd3072_PCA67$IID),]
EEPRS_gpticd3072_ICA67 = EEPRS_gpticd3072_ICA67[match(Trait_pheno_id,EEPRS_gpticd3072_ICA67$IID),]

test_data = data.table(PRScs = scale(PRScs$SCORE1_AVG),
                       EEPRS_Word2Vec100 = scale(EEPRS_Word2Vec100$SCORE1_AVG),
                       EEPRS_Word2Vec100_PCA30 = scale(EEPRS_Word2Vec100_PCA30$SCORE1_AVG),
                       EEPRS_Word2Vec100_ICA30 = scale(EEPRS_Word2Vec100_ICA30$SCORE1_AVG),
                       EEPRS_gpticd3072_PCA67 = scale(EEPRS_gpticd3072_PCA67$SCORE1_AVG),
                       EEPRS_gpticd3072_ICA67 = scale(EEPRS_gpticd3072_ICA67$SCORE1_AVG))
colnames(test_data) = c("PRScs","EEPRS_Word2Vec","EEPRS_Word2Vec_PCA","EEPRS_Word2Vec_ICA","EEPRS_GPT_PCA","EEPRS_GPT_ICA")


if (trait %in% continuous_trait){
## obtain the prs_table
prs_table = data.table(pop = c(pop), trait = c(trait), type = c("R2"),
                       PRScs = c(0), EEPRS_Word2Vec = c(0), EEPRS_Word2Vec_PCA = c(0), EEPRS_Word2Vec_ICA = c(0), EEPRS_GPT_PCA = c(0), EEPRS_GPT_ICA = c(0))

Trait_pheno$pheno = scale(Trait_pheno$pheno)
covariates = total_covariates[match(Trait_pheno_id,total_covariates$eid),]
pheno_covariates = cbind(Trait_pheno[,-1],covariates[,-1])
        
# null model in all individuals in UKBB dataset
linear_null = lm(pheno ~ . , data = pheno_covariates)
linear_null_summary = summary(linear_null)
linear_null_res2 = sum(linear_null_summary$residuals^2)

for (j in 1:6){
    data = pheno_covariates
    data$prs <- unlist(test_data[,..j])
    linear = lm(pheno ~ ., data=data)
    linear_summary=summary(linear)
    linear_summary_res2 = sum(linear_summary$residuals^2)

    prs_table[1,j+3] = 1 - linear_summary_res2/linear_null_res2
}

} else {
## obtain the prs_table
prs_table = data.table(pop = c(pop), trait = c(trait), type = c("AUC"),
                       PRScs = c(0), EEPRS_Word2Vec = c(0), EEPRS_Word2Vec_PCA = c(0), EEPRS_Word2Vec_ICA = c(0), EEPRS_GPT_PCA = c(0), EEPRS_GPT_ICA = c(0))

pheno = Trait_pheno[,-1]

for (j in 1:6){
    data = pheno
    data$prs <- unlist(test_data[,..j])
    glmfit = glm(pheno~prs, data=data,family=binomial(link="logit"))
    glmfit_prob = predict(glmfit, type="response")
    glmfit_auc = roc(data$pheno, glmfit_prob, quiet=T, plot=F)$auc
          
    prs_table[1,j+3] = glmfit_auc
}

}

prs_final_table = rbind(prs_final_table,prs_table)

}

write.table(prs_final_table,paste0("/gpfs/gibbs/pi/zhao/lx94/EEPRS/result/evaluation/EEPRS_benchmarking_R2_AUC.csv"),quote=F,sep='\t',row.names=F,col.names=T)
write_xlsx(prs_final_table,paste0("/gpfs/gibbs/pi/zhao/lx94/EEPRS/result/table/TableS3.xlsx"))