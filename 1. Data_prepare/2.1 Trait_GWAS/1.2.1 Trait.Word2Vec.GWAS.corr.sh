## Step0: # 1. Correlation calculation
job_file="/gpfs/gibbs/pi/zhao/lx94/EEPRS/code/glm/trait_word2vec100_corr.txt"
> $job_file  # Empty the job file if it already exists

pop=EUR
train_type=train

trait_list="HDL LDL TC logTG BMI WHR AgeSmk SmkI SmkC CigDay DrnkWk Glu2h HbA1c eGFR logCp AD Angina Asthma AF ADHD ASD BIP BrC CKD CAD HF IBD IS LuC MDD Osteoporosis OvC PAD PBC PSC PrC PAH RA SCZ SLE T2D homocysteine fibrinogen BUN IFCadjBMI ISIadjBMI FPadjBMI FG FI HOMA_BadjBMI HOMA_IRadjBMI"

for trait in ${trait_list}; do
for emd_i in {1..100}; do

# Construct the sumstats files for embeddings i and j (comma-separated)
file_i="/gpfs/gibbs/pi/zhao/lx94/EEPRS/data/summary_data/munge_data/${trait}_${pop}_ldsc.sumstats.gz"
file_j="/gpfs/gibbs/pi/zhao/lx94/EEPRS/data/embedding_data/munge_data/word2vec100_${train_type}_EUR_UKB_Embedding${emd_i}_ldsc.sumstats.gz"

out_prefix="/gpfs/gibbs/pi/zhao/lx94/EEPRS/result/h2_corr/corr/trait_embedding_corr/${trait}_word2vec100_${train_type}_EUR_UKB_Embedding${emd_i}_rg"

if [[ ! -e ${out_prefix}.log ]]; then
echo "module load miniconda; conda activate ldsc; python /gpfs/gibbs/pi/zhao/lx94/Software/ldsc/ldsc.py --rg "${file_i},${file_j}" --ref-ld-chr /gpfs/gibbs/pi/zhao/lx94/Data/ldsc/eur_w_ld_chr/ --w-ld-chr /gpfs/gibbs/pi/zhao/lx94/Data/ldsc/eur_w_ld_chr/ --out "${out_prefix}"" >> $job_file
fi

done
done

module load dSQ
cd /gpfs/gibbs/pi/zhao/lx94/EEPRS/code/glm/

dsq --job-file /gpfs/gibbs/pi/zhao/lx94/EEPRS/code/glm/trait_word2vec100_corr.txt --partition=scavenge,day --requeue --mem=30G --cpus-per-task=1 --ntasks=1 --nodes=1 --time=1-00:00:00 --mail-type=ALL
sbatch dsq-trait_word2vec100_corr-$(date +%Y-%m-%d).sh

## Step1: Get correlation function
get_trait_correlations <- function(trait, directory) {
  # Go to directory where the log files are
  setwd(directory)
  
  # Pattern to match files for this trait, e.g. CKD_word2vec100_train_EUR_UKB_EmbeddingNN_rgs.log
  pattern <- paste0("^", trait, "_word2vec100_train_EUR_UKB_Embedding[0-9]+_rg\\.log$")
  
  # List all files that match this trait pattern
  files <- list.files(
    pattern = pattern,
    full.names = TRUE
  )
  
  # A small helper to extract the embedding number from the filename
  get_embedding_number <- function(filename) {
    # For example from "CKD_word2vec100_train_EUR_UKB_Embedding12_rg.log" -> "12"
    sub(paste0("^", trait, "_word2vec100_train_EUR_UKB_Embedding([0-9]+)_rg\\.log$"),
        "\\1",
        basename(filename))
  }
  
  # Prepare a data frame to store all the results
  results <- data.frame(
    file      = character(),
    trait     = character(),
    emd_i     = integer(),
    rg        = numeric(),
    rg_se     = numeric(),
    z         = numeric(),
    p         = numeric(),
    stringsAsFactors = FALSE
  )
  
  for (f in files) {
    lines <- readLines(f, warn = FALSE)
    
    # Locate the header line: p1 p2 rg se z p h2_obs ...
    header_idx <- grep("^p1\\s+p2\\s+rg\\s+se\\s+z\\s+p\\s+h2_obs", lines)
    data_line_idx <- header_idx + 1  # The data line is usually the next line
    
    # Extract the numeric results if the line exists
    if (length(header_idx) == 1 && data_line_idx <= length(lines)) {
      columns <- strsplit(lines[data_line_idx], "\\s+")[[1]]
      
      rg_val    <- as.numeric(columns[3])
      rg_se_val <- as.numeric(columns[4])
      z_val     <- as.numeric(columns[5])
      p_val     <- as.numeric(columns[6])
    } else {
      rg_val    <- NA_real_
      rg_se_val <- NA_real_
      z_val     <- NA_real_
      p_val     <- NA_real_
    }
    
    # Extract embedding index from filename
    emd_str <- get_embedding_number(f)
    emd_val <- as.integer(emd_str)
    
    # Store in results data.frame
    results <- rbind(
      results,
      data.frame(
        file  = f,
        trait = trait,
        emd_i = emd_val,
        rg    = rg_val,
        rg_se = rg_se_val,
        z     = z_val,
        p     = p_val,
        stringsAsFactors = FALSE
      )
    )
  }

  return(results)
}

## Step2: Obtain all correlation results
directory_path <- "/gpfs/gibbs/pi/zhao/lx94/EEPRS/result/h2_corr/corr/trait_embedding_corr"

trait_list = c("HDL","LDL","TC","logTG","BMI","WHR","AgeSmk","SmkI","SmkC","CigDay","DrnkWk","Glu2h","HbA1c","eGFR","logCp",
                "AD","Angina","Asthma","AF","ADHD","ASD","BIP","BrC","CKD","CAD","HF","IBD","IS","LuC",
                "MDD","Osteoporosis","OvC","PAD","PBC","PSC","PrC","PAH","RA","SCZ","SLE","T2D",
                "homocysteine","fibrinogen","BUN","IFCadjBMI","ISIadjBMI","FPadjBMI","FG","FI","HOMA_BadjBMI","HOMA_IRadjBMI")

all_results <- data.frame()
for (trait in trait_list) {
  out <- get_trait_correlations(trait, directory_path)
  
  # Combine with a master table
  all_results <- rbind(all_results, out)
}

write.table(all_results, file=paste0("/gpfs/gibbs/pi/zhao/lx94/EEPRS/result/h2_corr/corr/trait_embedding_corr/UKB_word2vec100_train_trait_corr.txt"), row.names=F, col.names=T, quote=F, append=F, sep = "\t")


## Step3: Plot the correlation heatmap
library(dplyr)
library(tidyr)
library(pheatmap)

# obtain rg and p result
rg_data <- all_results %>%
  select(trait, emd_i, rg) %>%
  pivot_wider(names_from = emd_i, values_from = rg)

p_data <- all_results %>%
  select(trait, emd_i, p) %>%
  pivot_wider(names_from = emd_i, values_from = p)

# convert to data frames, then to matrices
rg_data <- as.data.frame(rg_data)
rownames(rg_data) <- rg_data$trait
rg_data$trait <- NULL
rg_mat <- as.matrix(rg_data)  # [51 x 100]

p_data <- as.data.frame(p_data)
rownames(p_data) <- p_data$trait
p_data$trait <- NULL
p_mat <- as.matrix(p_data)    # [51 x 100]

rg_rounded <- formatC(rg_mat, format = "f", digits = 2)
sig_star   <- ifelse(p_mat < 0.05, "*", "")

display_mat <- matrix(
  sig_star,
  nrow = nrow(rg_mat),
  ncol = ncol(rg_mat)
)
dimnames(display_mat) <- dimnames(rg_mat)


my_colors <- colorRampPalette(c("blue","white","red"))(50)

pheatmap(
  rg_mat,
  cluster_rows    = TRUE,
  cluster_cols    = TRUE,
  color           = my_colors,
  breaks          = seq(-1, 1, length.out = 51),
  display_numbers = display_mat,
  number_color    = "black",
  main            = "Genetic Correlation Estimates (51 Trait x Word2Vec)",
  treeheight_row  = 0,
  treeheight_col  = 0
)

## width: 1600 height: 800

## Step4: Obtain significant embeddings for each trait
library(dplyr)

all_traits <- data.frame(trait = unique(all_results$trait))

sig_table <- all_results %>%
  filter(p < 0.05) %>%
  group_by(trait) %>%
  summarize(
    sig_list = paste(sort(unique(emd_i)), collapse = " ")
  ) %>%
  right_join(all_traits, by = "trait") %>%
  mutate(sig_list = ifelse(is.na(sig_list), "", sig_list)) %>%
  ungroup()

colnames(sig_table) <- c("Trait", "Significant_Embeddings")
head(sig_table)

write.table(sig_table, file=paste0("/gpfs/gibbs/pi/zhao/lx94/EEPRS/result/h2_corr/corr/trait_embedding_corr/trait_significant_word2vec_embeddings.txt"), row.names=F, col.names=T, quote=F, append=F, sep = "\t")
