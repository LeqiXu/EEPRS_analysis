library(data.table)
library(PheWAS)
library(dplyr)
library(ggplot2)
library(stringr)
library(vroom)

## ------------------------------------------------------
## 1. Pre-load and prepare all standard datasets
## ------------------------------------------------------

## 1a. Phecode reference maps: https://wei-lab.app.vumc.org/phecode
new.vocabulary.map <- vroom::vroom("/gpfs/gibbs/pi/zhao/lx94/EEPRS/data/phecode/orig_phecode_map.csv", .name = janitor::make_clean_names, delim = ",", col_types = c(vocabulary_id = "c", code = "c", phecode = "c"))
new.annotate.phenotype.description <- vroom::vroom("/gpfs/gibbs/pi/zhao/lx94/EEPRS/data/phecode/original_phecodes_pheinfo.csv", .name = janitor::make_clean_names, delim = ",", col_types = c(phecode = "c", description = "c", groupnum = "double", group = "c", color = "c"))
new.exclusion.map <- vroom::vroom("/gpfs/gibbs/pi/zhao/lx94/EEPRS/data/phecode/original_phecodes_exclusion.csv", .name = janitor::make_clean_names, delim = ",", col_types = c(code = "c", exclusion_criteria = "c"))
new.sex.restriction <- vroom::vroom("/gpfs/gibbs/pi/zhao/lx94/EEPRS/data/phecode/original_phecodes_gender_restriction.csv", .name = janitor::make_clean_names, delim = ",", col_types = c(phecode = "c", male_only = "logical", female_only = "logical"))
new.rollup.map <- vroom::vroom("/gpfs/gibbs/pi/zhao/lx94/EEPRS/data/phecode/phecode_rollup_map_1.2024.csv", .name = janitor::make_clean_names, delim = ",", col_types = c(code = "c", phecode_unrolled = "c"))

colnames(new.annotate.phenotype.description)[1] <- "phenotype"

## 1b. Load ICD10 data
ICD10_data = fread("/gpfs/gibbs/pi/zhao/cz354/zhao-data/cz354/UKB/table/hesin_diag.txt")

## 1c. Subset to test set EIDs
test_data = fread("/gpfs/gibbs/pi/zhao/lx94/JointPRS/data/ukbb_data/ancestry_info/EUR_test_1over3_doubleid.tsv", header = FALSE)
ICD10_data_test = ICD10_data[which(ICD10_data$eid %in% test_data$V1),c("eid","ins_index","diag_icd10")]
ICD10_data_test$vocabulary_id = "ICD10CM"

## 1d. Create phenotypes
id.vocab.code.index <- ICD10_data_test[,c("eid","vocabulary_id","diag_icd10","ins_index")]
colnames(id.vocab.code.index) <- c("id","vocabulary_id","code","index")
phenotypes <- createPhenotypes(id.vocab.code.index, vocabulary.map = as.data.frame(new.vocabulary.map), rollup.map = as.data.frame(new.rollup.map), exclusion.map = as.data.frame(new.exclusion.map), sex.restriction = as.data.frame(new.sex.restriction))

## 1e. Load covariates (once), subset to test set
cov_choice = c("age_recruit","sex",paste0("PC",1:20))
total_covariates = fread("/gpfs/gibbs/pi/zhao/lx94/JointPRS/data/ukbb_data/cov_data/adjust.csv")
total_covariates = total_covariates %>% select(all_of(c("eid",cov_choice)))
test_total_covariates = total_covariates[which(total_covariates$eid %in% test_data$V1),]
colnames(test_total_covariates) = c("id",cov_choice)

## ------------------------------------------------------
## 2. Function to run PheWAS for a specific PRS file
## ------------------------------------------------------
run_phewas_for_prs <- function(
  prs_file,            # The .sscore path
  prs_label,           # A string label, e.g. "PRS_V25" or "gpticd3072_V1"
  phenotypes,          # from createPhenotypes()
  covariates           # test_total_covariates
) {
  # 2a. Read the PRS file
  # Adjust these column selections to match your actual .sscore files
  dat_prs <- fread(prs_file)
  
  # e.g., if your .sscore has columns: IID, SCORE1_AVG
  # or if it has multiple SCORES, adapt as needed
  # Here we assume "IID" is ID, "SCORE1_AVG" is the main PRS, etc.
  dat_prs[, PRS := scale(SCORE1_AVG)] 
  # rename or keep as is; just ensure the final name we merge on is "id"
  setnames(dat_prs, "IID", "id")

  # optionally keep only two columns for clarity
  dat_prs <- dat_prs[, c("id", "PRS")]
  colnames(dat_prs) = c("id",prs_label)

  # 2b. Run PheWAS
  results <- phewas(
    outcomes      = phenotypes,
    predictors    = dat_prs,
    covariates    = covariates,
    additive.genotypes = FALSE,
    significance.threshold = c("bonferroni")
  )
  
  return(results)
}

## ------------------------------------------------------
## 3. Loop over multiple PRS embeddings
## ------------------------------------------------------
for(type in c("word2vec100_ICA30","gpticd3072_ICA67")){

if (type == "word2vec100_ICA30"){
sig_embedding <- c(25,19,12,2,15,11,22,27,4,14,24,5,8,1,10,20,16,13,9,6,29,7)
title_type <- "Word2Vec_ICA"
}

if (type == "gpticd3072_ICA67"){
sig_embedding <- c(54,1,32,15,22,9,56,67,43,25,18,38,12,2,50,20,57,28,37,52,30,49,29,23,21,11,3,51,40,55,42,19,48,5,63,24,34,45,44,35,31,4,7,46)
title_type <- "GPT_ICA"
}

for(i in sig_embedding){
component=paste0("V",i)
this_label <- paste0(title_type,"_", component)
print(paste0(this_label, " Start!"))

# Build file path
prs_path <- paste0("/gpfs/gibbs/pi/zhao/lx94/EEPRS/result/Embedding_PRS/UKB_test_1over3_",type,"_train_EUR_UKB_", component, "_PRScsx_prs_EUR.sscore")

# Run PheWAS
tmp_res <- run_phewas_for_prs(
    prs_file   = prs_path,
    prs_label  = this_label,
    phenotypes = phenotypes,
    covariates = test_total_covariates
  )

write.table(tmp_res, file=paste0("/gpfs/gibbs/pi/zhao/lx94/EEPRS/result/Embedding_PRS/PheWAS/UKB_",type,"_",component,"_PheWAS.csv"), row.names=F, col.names=T, quote=F, append=F, sep = "\t")

print(paste0(this_label, " Complete!"))
}
}
