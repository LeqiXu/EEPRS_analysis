## obtain hm3 genetic data for EUR
#!/bin/bash
#SBATCH --partition=bigmem
#SBATCH --requeue
#SBATCH --mem=100G
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1 --nodes=1
#SBATCH --time=24:00:00
#SBATCH --mail-type=ALL
#SBATCH --job-name=ukbb_hm3_EUR
#SBATCH --output=out_ukbb_hm3_EUR.txt
module load PLINK/2

pop=EUR

plink2 --bfile /gpfs/gibbs/pi/zhao/lx94/JointPRS/data/ukbb_data/geno_data/${pop} \
--double-id \
--extract /gpfs/gibbs/pi/zhao/lx94/JointPRS/data/snplist_data/PRScsx_snplist.txt \
--keep /gpfs/gibbs/pi/zhao/lx94/JointPRS/data/ukbb_data/ancestry_info/${pop}_id.tsv \
--make-bed \
--out /gpfs/gibbs/pi/zhao/lx94/JointPRS/data/ukbb_data/geno_data/${pop}_hm3
