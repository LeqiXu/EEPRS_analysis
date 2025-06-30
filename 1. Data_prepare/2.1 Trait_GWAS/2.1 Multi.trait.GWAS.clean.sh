# 0. Create folder
trait_list="HDL LDL TC logTG BMI WHR AgeSmk SmkI SmkC CigDay DrnkWk Glu2h HbA1c eGFR logCp AD Angina Asthma AF ADHD ASD BIP BrC CKD CAD HF IBD IS LuC MDD Osteoporosis OvC PAD PBC PSC PrC PAH RA SCZ SLE T2D"

for trait in ${trait_list}; do
mkdir /gpfs/gibbs/pi/zhao/lx94/EEPRS/data/summary_data/PRScsx/MTAG_meta/${trait}
done

# 1.  MTAG GWAS
library(data.table)
library(stringr)

args <- commandArgs(trailingOnly = TRUE)

pop = "EUR"
trait = args[1]

if (trait %in% c("HDL","LDL","TC","logTG","IS")) {
  jiaqi_trait = tolower(trait)
} else if (trait == "WHR") {
  jiaqi_trait = "WHRadjBMI"
} else if (trait == "SmkI") {
  jiaqi_trait = "smk2"
} else if (trait == "SmkC") {
  jiaqi_trait = "SmkCes"
} else if (trait == "Glu2h") {
  jiaqi_trait = "2hGlu"
} else if (trait == "logCp") {
  jiaqi_trait = "CRP"
} else if (trait == "Angina") {
  jiaqi_trait = "I9_ANGINA"
} else if (trait == "AF") {
  jiaqi_trait = "af2"
} else if (trait == "BrC") {
  jiaqi_trait = "breast_cancer"
} else if (trait == "CAD") {
  jiaqi_trait = "cad_cardiogramc4a"
} else if (trait == "HF") {
  jiaqi_trait = "I9_HEARTFAIL"
} else if (trait == "LuC") {
  jiaqi_trait = "lung_cancer"
} else if (trait == "Osteoporosis") {
  jiaqi_trait = "M13_OSTEOPOROSIS"
} else if (trait == "OvC") {
  jiaqi_trait = "ovarian_cancer"
} else if (trait == "PrC") {
  jiaqi_trait = "prostate_cancer"
} else if (trait == "T2D") {
  jiaqi_trait = "t2d_mahajan2018b"
} else {
  jiaqi_trait = trait
}

input_dir <- paste0("/gpfs/gibbs/pi/zhao/jh2875/shared/prs/sumstats/MTAG_meta/",jiaqi_trait,"/")
output_dir <- paste0("/gpfs/gibbs/pi/zhao/lx94/EEPRS/data/summary_data/PRScsx/MTAG_meta/",trait,"/")

# Create output directory if it doesn't exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}
file_list <- list.files(input_dir, pattern = paste0(".*",jiaqi_trait,".*\\.txt$"), full.names = TRUE)

for (file in file_list) {

  # Read the file
  sumstats <- fread(file, header = TRUE)

  # Columns we need
  sumstat_PRScsx <- sumstats[, c("SNP", "A1", "A2", "mtag_beta", "mtag_pval")]
  colnames(sumstat_PRScsx) <- c("SNP", "A1", "A2", "BETA", "P")

  # Construct output filename
  base_name <- tools::file_path_sans_ext(basename(file))
  base_name <- str_replace(base_name,jiaqi_trait,trait)
  output_file <- paste0(output_dir, base_name, "_", pop, "_PRScsx.txt")

  # Write to disk
  write.table(sumstat_PRScsx, file = output_file, row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

  message("Processed ", basename(file), " -> ", output_file)
}

vim run_MATG_PRScsx.R


job_file="/gpfs/gibbs/pi/zhao/lx94/EEPRS/code/PRS/MATG_PRScsx.txt"
> $job_file  # Empty the job file if it already exists

trait_list="HDL LDL TC logTG BMI WHR AgeSmk SmkI SmkC CigDay DrnkWk Glu2h HbA1c eGFR logCp AD Angina Asthma AF ADHD ASD BIP BrC CKD CAD HF IBD IS LuC MDD Osteoporosis OvC PAD PBC PSC PrC PAH RA SCZ SLE T2D"

for trait in ${trait_list}; do

echo "module load miniconda; conda activate r_env; Rscript /gpfs/gibbs/pi/zhao/lx94/EEPRS/code/PRS/run_MATG_PRScsx.R ${trait}" >> $job_file

done


module load dSQ
cd /gpfs/gibbs/pi/zhao/lx94/EEPRS/code/PRS/

dsq --job-file MATG_PRScsx.txt --partition=day --requeue --mem=20G --cpus-per-task=1 --ntasks=1 --nodes=1 --time=24:00:00 --mail-type=ALL
sbatch dsq-MATG_PRScsx-$(date +%Y-%m-%d).sh