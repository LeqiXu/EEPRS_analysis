## 0. PRScs
pop=EUR

trait_list="HDL LDL TC logTG BMI WHR AgeSmk SmkI SmkC CigDay DrnkWk Glu2h HbA1c eGFR logCp AD Angina Asthma AF ADHD ASD BIP BrC CKD CAD HF IBD IS LuC MDD Osteoporosis OvC PAD PBC PSC PrC PAH RA SCZ SLE T2D"

for trait in ${trait_list}; do

if [[ ! -e "/gpfs/gibbs/pi/zhao/lx94/EEPRS/result/Trait_PRS/UKB_test_1over3_${trait}_PRScsx_prs_${pop}.sscore" ]]; then
        sbatch <<EOT
#!/bin/bash
#SBATCH --partition=scavenge,day
#SBATCH --requeue
#SBATCH --mem=100G
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1 --nodes=1
#SBATCH --time=24:00:00
#SBATCH --mail-type=ALL
#SBATCH --job-name=PRS_${trait}_PRScsx_prs_${pop}
#SBATCH --output=out_PRS_${trait}_PRScsx_prs_${pop}.txt

module load PLINK/2

cd /gpfs/gibbs/pi/zhao/lx94/EEPRS/result/Trait_PRS/

plink2 --bfile /gpfs/gibbs/pi/zhao/lx94/JointPRS/data/ukbb_data/geno_data/${pop}_hm3 \
--double-id \
--keep /gpfs/gibbs/pi/zhao/lx94/JointPRS/data/ukbb_data/ancestry_info/EUR_test_1over3_doubleid.tsv \
--threads 1 \
--score ${trait}_${pop}_PRScsx.txt \
--out UKB_test_1over3_${trait}_PRScsx_prs_${pop}

EOT
fi

done

## 1. MIXPRS_word2vec100
pop=EUR

trait_list="HDL LDL TC logTG BMI WHR AgeSmk SmkI SmkC CigDay DrnkWk Glu2h HbA1c eGFR logCp AD Angina Asthma AF ADHD ASD BIP BrC CKD CAD HF IBD IS LuC MDD Osteoporosis OvC PAD PBC PSC PrC PAH RA SCZ SLE T2D"

for trait in ${trait_list}; do

if [[ ! -e "/gpfs/gibbs/pi/zhao/lx94/EEPRS/result/MIXPRS/UKB_test_1over3_${trait}_EEPRS_word2vec100_prs_${pop}.sscore" ]]; then
        sbatch <<EOT
#!/bin/bash
#SBATCH --partition=scavenge,day
#SBATCH --requeue
#SBATCH --mem=100G
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1 --nodes=1
#SBATCH --time=24:00:00
#SBATCH --mail-type=ALL
#SBATCH --job-name=PRS_${trait}_EEPRS_word2vec100_prs_${pop}
#SBATCH --output=out_PRS_${trait}_EEPRS_word2vec100_prs_${pop}.txt

module load PLINK/2

cd /gpfs/gibbs/pi/zhao/lx94/EEPRS/result/MIXPRS/

plink2 --bfile /gpfs/gibbs/pi/zhao/lx94/JointPRS/data/ukbb_data/geno_data/${pop}_hm3 \
--double-id \
--keep /gpfs/gibbs/pi/zhao/lx94/JointPRS/data/ukbb_data/ancestry_info/EUR_test_1over3_doubleid.tsv \
--threads 1 \
--score ${trait}_word2vec100_${pop}_MIXPRS.txt \
--out UKB_test_1over3_${trait}_EEPRS_word2vec100_prs_${pop}

EOT
fi

done


## 2. MIXPRS_word2vec100_PCA30
pop=EUR

trait_list="HDL LDL TC logTG BMI WHR AgeSmk SmkI SmkC CigDay DrnkWk Glu2h HbA1c eGFR logCp AD Angina Asthma AF ADHD ASD BIP BrC CKD CAD HF IBD IS LuC MDD Osteoporosis OvC PAD PBC PSC PrC PAH RA SCZ SLE T2D"

for trait in ${trait_list}; do

if [[ ! -e "/gpfs/gibbs/pi/zhao/lx94/EEPRS/result/MIXPRS/UKB_test_1over3_${trait}_EEPRS_word2vec100_PCA30_prs_${pop}.sscore" ]]; then
        sbatch <<EOT
#!/bin/bash
#SBATCH --partition=scavenge,day
#SBATCH --requeue
#SBATCH --mem=100G
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1 --nodes=1
#SBATCH --time=24:00:00
#SBATCH --mail-type=ALL
#SBATCH --job-name=PRS_${trait}_EEPRS_word2vec100_PCA30_prs_${pop}
#SBATCH --output=out_PRS_${trait}_EEPRS_word2vec100_PCA30_prs_${pop}.txt

module load PLINK/2

cd /gpfs/gibbs/pi/zhao/lx94/EEPRS/result/MIXPRS/

plink2 --bfile /gpfs/gibbs/pi/zhao/lx94/JointPRS/data/ukbb_data/geno_data/${pop}_hm3 \
--double-id \
--keep /gpfs/gibbs/pi/zhao/lx94/JointPRS/data/ukbb_data/ancestry_info/EUR_test_1over3_doubleid.tsv \
--threads 1 \
--score ${trait}_word2vec100_PCA30_${pop}_MIXPRS.txt \
--out UKB_test_1over3_${trait}_EEPRS_word2vec100_PCA30_prs_${pop}

EOT
fi

done


## 2. MIXPRS_word2vec100_ICA30
pop=EUR

trait_list="HDL LDL TC logTG BMI WHR AgeSmk SmkI SmkC CigDay DrnkWk Glu2h HbA1c eGFR logCp AD Angina Asthma AF ADHD ASD BIP BrC CKD CAD HF IBD IS LuC MDD Osteoporosis OvC PAD PBC PSC PrC PAH RA SCZ SLE T2D"

for trait in ${trait_list}; do

if [[ ! -e "/gpfs/gibbs/pi/zhao/lx94/EEPRS/result/MIXPRS/UKB_test_1over3_${trait}_EEPRS_word2vec100_ICA30_prs_${pop}.sscore" ]]; then
        sbatch <<EOT
#!/bin/bash
#SBATCH --partition=scavenge,day
#SBATCH --requeue
#SBATCH --mem=100G
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1 --nodes=1
#SBATCH --time=24:00:00
#SBATCH --mail-type=ALL
#SBATCH --job-name=PRS_${trait}_EEPRS_word2vec100_ICA30_prs_${pop}
#SBATCH --output=out_PRS_${trait}_EEPRS_word2vec100_ICA30_prs_${pop}.txt

module load PLINK/2

cd /gpfs/gibbs/pi/zhao/lx94/EEPRS/result/MIXPRS/

plink2 --bfile /gpfs/gibbs/pi/zhao/lx94/JointPRS/data/ukbb_data/geno_data/${pop}_hm3 \
--double-id \
--keep /gpfs/gibbs/pi/zhao/lx94/JointPRS/data/ukbb_data/ancestry_info/EUR_test_1over3_doubleid.tsv \
--threads 1 \
--score ${trait}_word2vec100_ICA30_${pop}_MIXPRS.txt \
--out UKB_test_1over3_${trait}_EEPRS_word2vec100_ICA30_prs_${pop}

EOT
fi

done


## 5. MIXPRS_gpticd3072_PCA67
pop=EUR

trait_list="HDL LDL TC logTG BMI WHR AgeSmk SmkI SmkC CigDay DrnkWk Glu2h HbA1c eGFR logCp AD Angina Asthma AF ADHD ASD BIP BrC CKD CAD HF IBD IS LuC MDD Osteoporosis OvC PAD PBC PSC PrC PAH RA SCZ SLE T2D"

for trait in ${trait_list}; do

if [[ ! -e "/gpfs/gibbs/pi/zhao/lx94/EEPRS/result/MIXPRS/UKB_test_1over3_${trait}_EEPRS_gpticd3072_PCA67_prs_${pop}.sscore" ]]; then
        sbatch <<EOT
#!/bin/bash
#SBATCH --partition=scavenge,day
#SBATCH --requeue
#SBATCH --mem=100G
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1 --nodes=1
#SBATCH --time=24:00:00
#SBATCH --mail-type=ALL
#SBATCH --job-name=PRS_${trait}_EEPRS_gpticd3072_PCA67_prs_${pop}
#SBATCH --output=out_PRS_${trait}_EEPRS_gpticd3072_PCA67_prs_${pop}.txt

module load PLINK/2

cd /gpfs/gibbs/pi/zhao/lx94/EEPRS/result/MIXPRS/

plink2 --bfile /gpfs/gibbs/pi/zhao/lx94/JointPRS/data/ukbb_data/geno_data/${pop}_hm3 \
--double-id \
--keep /gpfs/gibbs/pi/zhao/lx94/JointPRS/data/ukbb_data/ancestry_info/EUR_test_1over3_doubleid.tsv \
--threads 1 \
--score ${trait}_gpticd3072_PCA67_${pop}_MIXPRS.txt \
--out UKB_test_1over3_${trait}_EEPRS_gpticd3072_PCA67_prs_${pop}

EOT
fi

done


## use tmp
pop=EUR

trait_list="HDL LDL TC logTG BMI WHR AgeSmk SmkI SmkC CigDay DrnkWk Glu2h HbA1c eGFR logCp AD Angina Asthma AF ADHD ASD BIP BrC CKD CAD HF IBD IS LuC MDD Osteoporosis OvC PAD PBC PSC PrC PAH RA SCZ SLE T2D"

for trait in ${trait_list}; do

if [[ ! -e "/gpfs/gibbs/pi/zhao/lx94/EEPRS/result/MIXPRS/UKB_test_1over3_${trait}_EEPRS_gpticd3072_PCA67_prs_${pop}.sscore" ]]; then
        sbatch <<EOT
#!/bin/bash
#SBATCH --partition=scavenge,day
#SBATCH --requeue
#SBATCH --mem=100G
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1 --nodes=1
#SBATCH --time=24:00:00
#SBATCH --mail-type=ALL
#SBATCH --job-name=PRS_${trait}_EEPRS_gpticd3072_PCA67_prs_${pop}
#SBATCH --output=out_PRS_${trait}_EEPRS_gpticd3072_PCA67_prs_${pop}.txt

module load PLINK/2

mkdir -p /tmp/lx94/\${SLURM_JOBID}
cd /tmp/lx94/\${SLURM_JOBID}
cp /gpfs/gibbs/pi/zhao/lx94/JointPRS/data/ukbb_data/geno_data/EUR_hm3.* .; \

plink2 --bfile ${pop}_hm3 \
--double-id \
--keep /gpfs/gibbs/pi/zhao/lx94/JointPRS/data/ukbb_data/ancestry_info/EUR_test_1over3_doubleid.tsv \
--threads 1 \
--score /gpfs/gibbs/pi/zhao/lx94/EEPRS/result/MIXPRS/${trait}_gpticd3072_PCA67_${pop}_MIXPRS.txt \
--out /gpfs/gibbs/pi/zhao/lx94/EEPRS/result/MIXPRS/UKB_test_1over3_${trait}_EEPRS_gpticd3072_PCA67_prs_${pop}

EOT
fi

done

## 6. MIXPRS_gpticd3072_ICA67
pop=EUR

trait_list="HDL LDL TC logTG BMI WHR AgeSmk SmkI SmkC CigDay DrnkWk Glu2h HbA1c eGFR logCp AD Angina Asthma AF ADHD ASD BIP BrC CKD CAD HF IBD IS LuC MDD Osteoporosis OvC PAD PBC PSC PrC PAH RA SCZ SLE T2D"

for trait in ${trait_list}; do

if [[ ! -e "/gpfs/gibbs/pi/zhao/lx94/EEPRS/result/MIXPRS/UKB_test_1over3_${trait}_EEPRS_gpticd3072_ICA67_prs_${pop}.sscore" ]]; then
        sbatch <<EOT
#!/bin/bash
#SBATCH --partition=scavenge,day
#SBATCH --requeue
#SBATCH --mem=100G
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1 --nodes=1
#SBATCH --time=24:00:00
#SBATCH --mail-type=ALL
#SBATCH --job-name=PRS_${trait}_EEPRS_gpticd3072_ICA67_prs_${pop}
#SBATCH --output=out_PRS_${trait}_EEPRS_gpticd3072_ICA67_prs_${pop}.txt

module load PLINK/2

cd /gpfs/gibbs/pi/zhao/lx94/EEPRS/result/MIXPRS/

plink2 --bfile /gpfs/gibbs/pi/zhao/lx94/JointPRS/data/ukbb_data/geno_data/${pop}_hm3 \
--double-id \
--keep /gpfs/gibbs/pi/zhao/lx94/JointPRS/data/ukbb_data/ancestry_info/EUR_test_1over3_doubleid.tsv \
--threads 1 \
--score ${trait}_gpticd3072_ICA67_${pop}_MIXPRS.txt \
--out UKB_test_1over3_${trait}_EEPRS_gpticd3072_ICA67_prs_${pop}

EOT
fi

done


## use tmp
pop=EUR

trait_list="HDL LDL TC logTG BMI WHR AgeSmk SmkI SmkC CigDay DrnkWk Glu2h HbA1c eGFR logCp AD Angina Asthma AF ADHD ASD BIP BrC CKD CAD HF IBD IS LuC MDD Osteoporosis OvC PAD PBC PSC PrC PAH RA SCZ SLE T2D"

for trait in ${trait_list}; do

if [[ ! -e "/gpfs/gibbs/pi/zhao/lx94/EEPRS/result/MIXPRS/UKB_test_1over3_${trait}_EEPRS_gpticd3072_ICA67_prs_${pop}.sscore" ]]; then
        sbatch <<EOT
#!/bin/bash
#SBATCH --partition=scavenge,day
#SBATCH --requeue
#SBATCH --mem=100G
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1 --nodes=1
#SBATCH --time=24:00:00
#SBATCH --mail-type=ALL
#SBATCH --job-name=PRS_${trait}_EEPRS_gpticd3072_ICA67_prs_${pop}
#SBATCH --output=out_PRS_${trait}_EEPRS_gpticd3072_ICA67_prs_${pop}.txt

module load PLINK/2

mkdir -p /tmp/lx94/\${SLURM_JOBID}
cd /tmp/lx94/\${SLURM_JOBID}
cp /gpfs/gibbs/pi/zhao/lx94/JointPRS/data/ukbb_data/geno_data/EUR_hm3.* .; \

plink2 --bfile ${pop}_hm3 \
--double-id \
--keep /gpfs/gibbs/pi/zhao/lx94/JointPRS/data/ukbb_data/ancestry_info/EUR_test_1over3_doubleid.tsv \
--threads 1 \
--score /gpfs/gibbs/pi/zhao/lx94/EEPRS/result/MIXPRS/${trait}_gpticd3072_ICA67_${pop}_MIXPRS.txt \
--out /gpfs/gibbs/pi/zhao/lx94/EEPRS/result/MIXPRS/UKB_test_1over3_${trait}_EEPRS_gpticd3072_ICA67_prs_${pop}

EOT
fi

done