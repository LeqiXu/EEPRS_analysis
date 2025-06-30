# Step1: Data preparation
library(data.table)
library(stringr)
library(dplyr)
library(pROC)

args <- commandArgs(trailingOnly = TRUE)

pop = "EUR"
trait = args[1]

continuous_trait = c("HDL","LDL","TC","logTG","BMI","WHR","AgeSmk","CigDay","DrnkWk","Glu2h","HbA1c","eGFR","logCp")

## test phenotype
Trait_pheno <- fread(paste0("/gpfs/gibbs/pi/zhao/lx94/EEPRS/data/ukbb_pheno_data/ukbb_",trait,".tsv"))
Trait_pheno_name = c("eid",trait)
Trait_pheno = Trait_pheno[,..Trait_pheno_name]
colnames(Trait_pheno) = c("eid","pheno")
Trait_pheno = na.omit(Trait_pheno)

## depend on the trait type and pop type, we read the score accordingly
single_trait_PRS = fread(paste0("/gpfs/gibbs/pi/zhao/lx94/EEPRS/result/Trait_PRS/UKB_all_",trait,"_PRScsx_prs_",pop,".sscore"))
single_trait_PRS = single_trait_PRS[,c("IID","SCORE1_AVG")]
colnames(single_trait_PRS) = c("IID",trait)

multi_trait_PRS_file = list.files(path = paste0("/gpfs/gibbs/pi/zhao/lx94/EEPRS/result/MTAG_meta_PRS/",trait,"/"), pattern = "UKB_all_mtag.*.sscore$", full.names = TRUE)
multi_trait_PRS = single_trait_PRS[,c("IID")]
for (file in multi_trait_PRS_file){
  file_base_name = sub("_EUR_PRScsx_EUR.sscore$", "", basename(file))
  
  sub_multi_trait_PRS = fread(file)
  sub_multi_trait_PRS = sub_multi_trait_PRS[,c("IID","SCORE1_AVG")]
  colnames(sub_multi_trait_PRS) = c("IID",file_base_name)

  multi_trait_PRS = merge(multi_trait_PRS,sub_multi_trait_PRS,by = c("IID"))
}

## match the IID for pheno and PRS for both training and evaluation data
Trait_pheno = Trait_pheno[which(Trait_pheno$eid %in% single_trait_PRS$IID),]
Train_id = fread(paste0("/gpfs/gibbs/pi/zhao/lx94/JointPRS/data/ukbb_data/ancestry_info/EUR_train_2over3_doubleid.tsv"))

Trait_pheno_train = Trait_pheno[which(Trait_pheno$eid %in% Train_id$V1),]
Trait_pheno_test = Trait_pheno[which(!(Trait_pheno$eid %in% Train_id$V1)),]

if (trait %in% continuous_trait) {
Trait_pheno_train$pheno = scale(Trait_pheno_train$pheno)
Trait_pheno_test$pheno = scale(Trait_pheno_test$pheno)
}

Trait_pheno_train_id = Trait_pheno_train$eid
Trait_pheno_test_id = Trait_pheno_test$eid

Trait_pheno_train = Trait_pheno_train[,-1]
Trait_pheno_test = Trait_pheno_test[,-1]

multi_trait_PRS_train = multi_trait_PRS[match(Trait_pheno_train_id,multi_trait_PRS$IID),-1]
multi_trait_PRS_train = as.data.table(scale(multi_trait_PRS_train))
multi_trait_PRS_test = multi_trait_PRS[match(Trait_pheno_test_id,multi_trait_PRS$IID),-1]
multi_trait_PRS_test = as.data.table(scale(multi_trait_PRS_test))

# Step2: PRS model training
library(glmnet)  # For elastic net
library(data.table)

# For Train
X_train <- model.matrix( ~ 0 + ., data = multi_trait_PRS_train)
y_train <- Trait_pheno_train$pheno

if (trait %in% continuous_trait) {
  family_choice <- "gaussian"
} else {
  family_choice <- "binomial"  # logistic
}

set.seed(123)  # reproducibility
cvfit <- cv.glmnet(
  x = X_train,
  y = y_train,
  alpha = 0.5,
  family = family_choice,
  intercept = FALSE  # no intercept in the penalized model
)
best_lambda <- cvfit$lambda.min

# For Test
X_test <- model.matrix(~ 0 + ., data = multi_trait_PRS_test)
MTAG_PRS_test <- predict(cvfit, newx = X_test, s = best_lambda, type = "link")
MTAG_PRS_table = data.table(IID = Trait_pheno_test_id, MTAG_PRS = as.numeric(MTAG_PRS_test))

write.table(MTAG_PRS_table,paste0("/gpfs/gibbs/pi/zhao/lx94/EEPRS/result/MTAG_meta_PRS/",trait,"/UKB_test_1over3_",trait,"_MTAG_PRS.txt"),quote=F,sep='\t',row.names=F,col.names=T)

vim train_MTAG_PRS.R

job_file="/gpfs/gibbs/pi/zhao/lx94/EEPRS/code/PRS/train_MTAG_PRS.txt"
> $job_file  # Empty the job file if it already exists

trait_list="HDL LDL TC logTG BMI WHR AgeSmk SmkI SmkC CigDay DrnkWk Glu2h HbA1c eGFR logCp AD Angina Asthma AF ADHD ASD BIP BrC CKD CAD HF IBD IS LuC MDD Osteoporosis OvC PAD PBC PSC PrC PAH RA SCZ SLE T2D"

for trait in ${trait_list}; do

echo "module load miniconda; conda activate r_env; Rscript /gpfs/gibbs/pi/zhao/lx94/EEPRS/code/PRS/train_MTAG_PRS.R ${trait}" >> $job_file

done

module load dSQ
cd /gpfs/gibbs/pi/zhao/lx94/EEPRS/code/PRS/

dsq --job-file train_MTAG_PRS.txt --partition=day --requeue --mem=20G --cpus-per-task=1 --ntasks=1 --nodes=1 --time=24:00:00 --mail-type=ALL
sbatch dsq-train_MTAG_PRS-$(date +%Y-%m-%d).sh
