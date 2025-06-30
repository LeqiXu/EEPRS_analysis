## Word2Vec_ICA
pop=EUR

for i in {1..30}; do

if [[ ! -e "/gpfs/gibbs/pi/zhao/lx94/EEPRS/result/Embedding_PRS/UKB_test_1over3_word2vec100_ICA30_train_${pop}_UKB_V${i}_PRScsx_prs_${pop}.sscore" ]]; then
        sbatch <<EOT
#!/bin/bash
#SBATCH --partition=scavenge,day
#SBATCH --requeue
#SBATCH --mem=100G
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1 --nodes=1
#SBATCH --time=24:00:00
#SBATCH --mail-type=ALL
#SBATCH --job-name=PRS_word2vec100_ICA30_train_${pop}_UKB_V${i}_PRScsx_prs_${pop}
#SBATCH --output=out_PRS_word2vec100_ICA30_train_${pop}_UKB_V${i}_PRScsx_prs_${pop}.txt

module load PLINK/2

cd /gpfs/gibbs/pi/zhao/lx94/EEPRS/result/Embedding_PRS/

plink2 --bfile /gpfs/gibbs/pi/zhao/lx94/JointPRS/data/ukbb_data/geno_data/${pop}_hm3 \
--double-id \
--keep /gpfs/gibbs/pi/zhao/lx94/JointPRS/data/ukbb_data/ancestry_info/EUR_test_1over3_doubleid.tsv \
--threads 1 \
--score word2vec100_ICA30_train_${pop}_UKB_V${i}_PRScsx.txt \
--out UKB_test_1over3_word2vec100_ICA30_train_${pop}_UKB_V${i}_PRScsx_prs_${pop}

EOT
fi

done


## Word2Vec_ICA
pop=EUR

for i in {1..67}; do

if [[ ! -e "/gpfs/gibbs/pi/zhao/lx94/EEPRS/result/Embedding_PRS/UKB_test_1over3_gpticd3072_ICA67_train_${pop}_UKB_V${i}_PRScsx_prs_${pop}.sscore" ]]; then
        sbatch <<EOT
#!/bin/bash
#SBATCH --partition=scavenge,day
#SBATCH --requeue
#SBATCH --mem=100G
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1 --nodes=1
#SBATCH --time=24:00:00
#SBATCH --mail-type=ALL
#SBATCH --job-name=PRS_gpticd3072_ICA67_train_${pop}_UKB_V${i}_PRScsx_prs_${pop}
#SBATCH --output=out_PRS_gpticd3072_ICA67_train_${pop}_UKB_V${i}_PRScsx_prs_${pop}.txt

module load PLINK/2

cd /gpfs/gibbs/pi/zhao/lx94/EEPRS/result/Embedding_PRS/

plink2 --bfile /gpfs/gibbs/pi/zhao/lx94/JointPRS/data/ukbb_data/geno_data/${pop}_hm3 \
--double-id \
--keep /gpfs/gibbs/pi/zhao/lx94/JointPRS/data/ukbb_data/ancestry_info/EUR_test_1over3_doubleid.tsv \
--threads 1 \
--score gpticd3072_ICA67_train_${pop}_UKB_V${i}_PRScsx.txt \
--out UKB_test_1over3_gpticd3072_ICA67_train_${pop}_UKB_V${i}_PRScsx_prs_${pop}

EOT
fi

done