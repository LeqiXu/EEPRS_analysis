# PRS evaluation
library(data.table)
library(stringr)
library(dplyr)
library(tidyr)
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

cv_results_all_traits = list()

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
if (trait %in% continuous_trait) {
Trait_pheno$pheno = scale(Trait_pheno$pheno)
}

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

## Merge covariates
cov_dt = total_covariates[match(Trait_pheno_id, total_covariates$eid), ]
merged_dt = cbind(test_data, cov_dt[ , -1]) # drop 'eid' col from cov_dt
merged_dt = cbind(Trait_pheno[, -1], merged_dt) # drop 'eid' col from Trait_pheno

#=== Create folds ===#
set.seed(42)  # for reproducible folds
n_ind = nrow(merged_dt)
n_fold = 4
folds = sample(rep(1:n_fold, length.out = n_ind))
merged_dt[, fold := folds]

# We'll store results for each fold here:
fold_results = data.table()

#=== Fold Loop ===#
for (f in 1:n_fold) {
#---- (a) Split train vs. test ----#
train_dt = merged_dt[fold != f]
test_dt  = merged_dt[fold == f]

train_performance = c()

# Because you have continuous OR binary traits, let's handle that:
if (trait %in% continuous_trait) {
#  -- continuous trait: compare R^2 against a baseline model with covariates

cov_only_formula = as.formula(paste("pheno ~", paste(cov_choice, collapse = " + ")))
null_mod = lm(cov_only_formula, data = train_dt)
null_ssr = sum(resid(null_mod)^2)

all_prs_cols = c("PRScs", "EEPRS_Word2Vec","EEPRS_Word2Vec_PCA","EEPRS_Word2Vec_ICA","EEPRS_GPT_PCA","EEPRS_GPT_ICA")

for (prs_col in all_prs_cols) {
    form_str = paste("pheno ~", paste(cov_choice, collapse = " + "), "+", prs_col)
    fit = lm(as.formula(form_str), data=train_dt)
    fit_ssr = sum(resid(fit)^2)
    r2 = 1 - (fit_ssr / null_ssr)
    train_performance[prs_col] = r2
}
      
} else {
all_prs_cols = c("PRScs", "EEPRS_Word2Vec","EEPRS_Word2Vec_PCA","EEPRS_Word2Vec_ICA","EEPRS_GPT_PCA","EEPRS_GPT_ICA")

for (prs_col in all_prs_cols) {
    form_str = paste("pheno ~", prs_col)
    fit = glm(as.formula(form_str), data=train_dt, family=binomial(link="logit"))
    pred_probs = predict(fit, newdata=train_dt, type="response")
    auc_val = pROC::roc(train_dt$pheno, pred_probs, quiet=TRUE)$auc
    train_performance[prs_col] = auc_val
}
}

eeprs_cols = c("EEPRS_Word2Vec","EEPRS_Word2Vec_PCA","EEPRS_Word2Vec_ICA","EEPRS_GPT_PCA","EEPRS_GPT_ICA")
eeprs_best = eeprs_cols[ which.max(train_performance[eeprs_cols]) ]

test_perf = data.table(fold = f,trait = trait,method = c("PRScs", "EEPRS_optimal"),performance = NA_real_)

if (trait %in% continuous_trait) {
    null_mod_test = lm(cov_only_formula, data=test_dt)
    null_ssr_test = sum(resid(null_mod_test)^2)

    prs_fit_test = lm(as.formula(paste("pheno ~", paste(cov_choice, collapse = " + "), "+ PRScs")), data=test_dt)
    prs_ssr = sum(resid(prs_fit_test)^2)
    test_perf[method=="PRScs", performance := 1 - prs_ssr/null_ssr_test]

} else {
    prs_fit_test = glm(as.formula(paste("pheno ~", " PRScs")), data=test_dt, family=binomial(link="logit"))
    prs_probs = predict(prs_fit_test, newdata=test_dt, type="response")
    auc_val = pROC::roc(test_dt$pheno, prs_probs, quiet=TRUE)$auc
    test_perf[method=="PRScs", performance := auc_val]
}

if (trait %in% continuous_trait) {
    null_ssr_test = sum(resid(null_mod_test)^2)
    eeprs_fit_test = lm(as.formula(paste("pheno ~", paste(cov_choice, collapse = " + "), "+", eeprs_best)), data=test_dt)
    eeprs_ssr = sum(resid(eeprs_fit_test)^2)
    eeprs_r2 = 1 - eeprs_ssr/null_ssr_test
    test_perf[method=="EEPRS_optimal", performance := eeprs_r2]
} else {
    eeprs_fit_test = glm(as.formula(paste("pheno ~", eeprs_best)), data=test_dt, family=binomial(link="logit"))
    eeprs_probs = predict(eeprs_fit_test, newdata=test_dt, type="response")
    eeprs_auc = pROC::roc(test_dt$pheno, eeprs_probs, quiet=TRUE)$auc
    test_perf[method=="EEPRS_optimal", performance := eeprs_auc]
}

fold_results = rbind(fold_results, test_perf)
}

trait_cv_summary = fold_results[ , .(mean_performance = mean(performance)), by=.(trait, method)]
cv_results_all_traits[[trait]] = trait_cv_summary
}

#=== Combine everything & look at results ===#
final_cv_results = rbindlist(cv_results_all_traits)
print(final_cv_results)

write.table(final_cv_results,paste0("/gpfs/gibbs/pi/zhao/lx94/EEPRS/result/evaluation/EEPRS_optimal_R2_AUC.csv"),quote=F,sep='\t',row.names=F,col.names=T)

wide_table <- pivot_wider(final_cv_results, 
                          names_from = method, 
                          values_from = mean_performance)
write_xlsx(wide_table,paste0("/gpfs/gibbs/pi/zhao/lx94/EEPRS/result/table/TableS4.xlsx"))