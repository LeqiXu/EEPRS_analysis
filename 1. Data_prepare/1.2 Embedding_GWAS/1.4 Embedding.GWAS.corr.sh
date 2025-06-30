## Word2Vec
# 1. Correlation calculation
job_file="/gpfs/gibbs/pi/zhao/lx94/EEPRS/code/glm/word2vec100_corr.txt"
> $job_file  # Empty the job file if it already exists

train_type=train

for i in {1..99}; do
for j in $(seq $((i+1)) 100); do

# Construct the sumstats files for embeddings i and j (comma-separated)
file_i="/gpfs/gibbs/pi/zhao/lx94/EEPRS/data/embedding_data/munge_data/word2vec100_${train_type}_EUR_UKB_Embedding${i}_ldsc.sumstats.gz"
file_j="/gpfs/gibbs/pi/zhao/lx94/EEPRS/data/embedding_data/munge_data/word2vec100_${train_type}_EUR_UKB_Embedding${j}_ldsc.sumstats.gz"
    
out_prefix="/gpfs/gibbs/pi/zhao/lx94/EEPRS/result/h2_corr/corr/embedding_embedding_corr/word2vec100_${train_type}_EUR_UKB_Embedding${i}_${j}_rg"

if [[ ! -e ${out_prefix}.log ]]; then
echo "module load miniconda; conda activate ldsc; python /gpfs/gibbs/pi/zhao/lx94/Software/ldsc/ldsc.py --rg "${file_i},${file_j}" --ref-ld-chr /gpfs/gibbs/pi/zhao/lx94/Data/ldsc/eur_w_ld_chr/ --w-ld-chr /gpfs/gibbs/pi/zhao/lx94/Data/ldsc/eur_w_ld_chr/ --out "${out_prefix}"" >> $job_file
fi
    
done
done

module load dSQ
cd /gpfs/gibbs/pi/zhao/lx94/EEPRS/code/glm/

dsq --job-file /gpfs/gibbs/pi/zhao/lx94/EEPRS/code/glm/word2vec100_corr.txt --partition=scavenge,day --requeue --mem=30G --cpus-per-task=1 --ntasks=1 --nodes=1 --time=1-00:00:00 --mail-type=ALL
sbatch dsq-word2vec100_corr-$(date +%Y-%m-%d).sh

# 2. Obtain correlation table
library(pheatmap)

setwd("/gpfs/gibbs/pi/zhao/lx94/EEPRS/result/h2_corr/corr/embedding_embedding_corr")

files <- list.files(
  pattern = "^word2vec100_train_EUR_UKB_Embedding[0-9]+_[0-9]+_rg\\.log$",
  full.names = TRUE
)

# Helper: extract i, j from filename
get_i_j_from_filename <- function(fname) {
  bname <- basename(fname)
  pattern <- "^word2vec100_train_EUR_UKB_Embedding([0-9]+)_([0-9]+)_rg\\.log$"
  captures <- sub(pattern, "\\1 \\2", bname)  # e.g. "13 27"
  parts <- strsplit(captures, " ")[[1]]
  i_val <- as.integer(parts[1])
  j_val <- as.integer(parts[2])
  return(c(i_val, j_val))
}

# Helper: read `rg` and `p` from the log
get_rg_p_from_logfile <- function(fname) {
  lines <- readLines(fname, warn = FALSE)
  header_idx <- grep("^p1\\s+p2\\s+rg\\s+se\\s+z\\s+p", lines)
  if (length(header_idx) == 1) {
    data_line_idx <- header_idx + 1
    if (data_line_idx <= length(lines)) {
      columns <- strsplit(lines[data_line_idx], "\\s+")[[1]]
      rg_val <- as.numeric(columns[3])  # 3rd column = rg
      p_val  <- as.numeric(columns[6])  # 6th column = p
      return(c(rg_val, p_val))
    }
  }
  return(c(NA_real_, NA_real_))
}

# Build a long data frame with i, j, rg, p
results <- data.frame(
  i  = integer(),
  j  = integer(),
  rg = numeric(),
  p  = numeric(),
  stringsAsFactors = FALSE
)

for (f in files) {
  ij       <- get_i_j_from_filename(f)
  i_val    <- ij[1]
  j_val    <- ij[2]
  rg_p_val <- get_rg_p_from_logfile(f)
  
  results <- rbind(
    results,
    data.frame(
      i  = i_val,
      j  = j_val,
      rg = rg_p_val[1],
      p  = rg_p_val[2],
      stringsAsFactors = FALSE
    )
  )
}

results$embedding_approach = "Word2Vec"
results

write.table(results, file=paste0("/gpfs/gibbs/pi/zhao/lx94/EEPRS/result/h2_corr/corr/embedding_embedding_corr/UKB_word2vec100_train_corr.txt"), row.names=F, col.names=T, quote=F, append=F, sep = "\t")

# Find maximum embedding index
max_dim <- max(results$i, results$j)

# Initialize matrices for rg and p
rg_matrix <- matrix(NA_real_, nrow = max_dim, ncol = max_dim)
p_matrix <- matrix(NA_real_, nrow = max_dim, ncol = max_dim)

# Fill the upper and lower triangles
for (k in seq_len(nrow(results))) {
  i_val <- results$i[k]
  j_val <- results$j[k]
  rg_val <- results$rg[k]
  p_val  <- results$p[k]
  
  rg_matrix[i_val, j_val] <- rg_val
  rg_matrix[j_val, i_val] <- rg_val
  
  p_matrix[i_val, j_val] <- p_val
  p_matrix[j_val, i_val] <- p_val
}

# Set diagonal to 1 for rg, 0 for p (or NA -- your choice)
diag(rg_matrix) <- 1
diag(p_matrix) <- 0

rownames(rg_matrix) <- paste0(1:max_dim)
colnames(rg_matrix) <- paste0(1:max_dim)
rownames(p_matrix)  <- paste0(1:max_dim)
colnames(p_matrix)  <- paste0(1:max_dim)

# Plot genetic correlation
sig_level <- 0.05  # adjust if you want a stricter threshold
sig_labels <- ifelse(p_matrix < sig_level, "*", "")

# Now we create a heatmap of rg_matrix, and supply `sig_labels` to display_numbers.
pheatmap(
  rg_matrix[1:100, 1:100],                
  cluster_rows    = TRUE,                 
  cluster_cols    = TRUE,                 
  display_numbers = sig_labels[1:100, 1:100],  
  number_color    = "black",
  main            = "Genetic Correlation Estimates (Word2Vec)",
  legend          = TRUE,
  # Set tree heights to zero to hide dendrogram lines
  treeheight_row  = 0,
  treeheight_col  = 0
)

# width: 1800 height:1800


## GPT
# 1. Correlation calculation
job_file="/gpfs/gibbs/pi/zhao/lx94/EEPRS/code/glm/gpticd3072_corr.txt"
> $job_file  # Empty the job file if it already exists

train_type=train

for i in {1..3071}; do
for j in $(seq $((i+1)) 3072); do

# Construct the sumstats files for embeddings i and j (comma-separated)
file_i="/gpfs/gibbs/pi/zhao/lx94/EEPRS/data/embedding_data/munge_data/gpticd3072_${train_type}_EUR_UKB_Embedding${i}_ldsc.sumstats.gz"
file_j="/gpfs/gibbs/pi/zhao/lx94/EEPRS/data/embedding_data/munge_data/gpticd3072_${train_type}_EUR_UKB_Embedding${j}_ldsc.sumstats.gz"
    
out_prefix="/gpfs/gibbs/pi/zhao/lx94/EEPRS/result/h2_corr/corr/embedding_embedding_corr/gpticd3072_${train_type}_EUR_UKB_Embedding${i}_${j}_rg"

if [[ ! -e ${out_prefix}.log ]]; then
echo "module load miniconda; conda activate ldsc; python /gpfs/gibbs/pi/zhao/lx94/Software/ldsc/ldsc.py --rg "${file_i},${file_j}" --ref-ld-chr /gpfs/gibbs/pi/zhao/lx94/Data/ldsc/eur_w_ld_chr/ --w-ld-chr /gpfs/gibbs/pi/zhao/lx94/Data/ldsc/eur_w_ld_chr/ --out "${out_prefix}"" >> $job_file
fi
    
done
done

module load dSQ
cd /gpfs/gibbs/pi/zhao/lx94/EEPRS/code/glm/

dsq --job-file /gpfs/gibbs/pi/zhao/lx94/EEPRS/code/glm/gpticd3072_corr.txt --partition=scavenge,day --requeue --mem=30G --cpus-per-task=1 --ntasks=1 --nodes=1 --time=1-00:00:00 --mail-type=ALL
sbatch dsq-gpticd3072_corr-$(date +%Y-%m-%d).sh

# 2. Obtain correlation table
library(pheatmap)

setwd("/gpfs/gibbs/pi/zhao/lx94/EEPRS/result/h2_corr/corr/embedding_embedding_corr")

files <- list.files(
  pattern = "^gpticd3072_train_EUR_UKB_Embedding[0-9]+_[0-9]+_rg\\.log$",
  full.names = TRUE
)

# Helper: extract i, j from filename
get_i_j_from_filename <- function(fname) {
  bname <- basename(fname)
  pattern <- "^gpticd3072_train_EUR_UKB_Embedding([0-9]+)_([0-9]+)_rg\\.log$"
  captures <- sub(pattern, "\\1 \\2", bname)  # e.g. "13 27"
  parts <- strsplit(captures, " ")[[1]]
  i_val <- as.integer(parts[1])
  j_val <- as.integer(parts[2])
  return(c(i_val, j_val))
}

# Helper: read `rg` and `p` from the log
get_rg_p_from_logfile <- function(fname) {
  lines <- readLines(fname, warn = FALSE)
  header_idx <- grep("^p1\\s+p2\\s+rg\\s+se\\s+z\\s+p", lines)
  if (length(header_idx) == 1) {
    data_line_idx <- header_idx + 1
    if (data_line_idx <= length(lines)) {
      columns <- strsplit(lines[data_line_idx], "\\s+")[[1]]
      rg_val <- as.numeric(columns[3])  # 3rd column = rg
      p_val  <- as.numeric(columns[6])  # 6th column = p
      return(c(rg_val, p_val))
    }
  }
  return(c(NA_real_, NA_real_))
}

# Build a long data frame with i, j, rg, p
results <- data.frame(
  i  = integer(),
  j  = integer(),
  rg = numeric(),
  p  = numeric(),
  stringsAsFactors = FALSE
)

for (f in files) {
  ij       <- get_i_j_from_filename(f)
  i_val    <- ij[1]
  j_val    <- ij[2]
  rg_p_val <- get_rg_p_from_logfile(f)
  
  results <- rbind(
    results,
    data.frame(
      i  = i_val,
      j  = j_val,
      rg = rg_p_val[1],
      p  = rg_p_val[2],
      stringsAsFactors = FALSE
    )
  )
}

results$embedding_approach = "GPT"
results

write.table(results, file=paste0("/gpfs/gibbs/pi/zhao/lx94/EEPRS/result/h2_corr/corr/embedding_embedding_corr/UKB_gpticd3072_train_corr.txt"), row.names=F, col.names=T, quote=F, append=F, sep = "\t")

# Find maximum embedding index
max_dim <- max(results$i, results$j)

# Initialize matrices for rg and p
rg_matrix <- matrix(NA_real_, nrow = max_dim, ncol = max_dim)
p_matrix <- matrix(NA_real_, nrow = max_dim, ncol = max_dim)

# Fill the upper and lower triangles
for (k in seq_len(nrow(results))) {
  i_val <- results$i[k]
  j_val <- results$j[k]
  rg_val <- results$rg[k]
  p_val  <- results$p[k]
  
  rg_matrix[i_val, j_val] <- rg_val
  rg_matrix[j_val, i_val] <- rg_val
  
  p_matrix[i_val, j_val] <- p_val
  p_matrix[j_val, i_val] <- p_val
}

# Set diagonal to 1 for rg, 0 for p (or NA -- your choice)
diag(rg_matrix) <- 1
diag(p_matrix) <- 0

rownames(rg_matrix) <- paste0(1:max_dim)
colnames(rg_matrix) <- paste0(1:max_dim)
rownames(p_matrix)  <- paste0(1:max_dim)
colnames(p_matrix)  <- paste0(1:max_dim)

# Plot genetic correlation
sig_level <- 0.05  # adjust if you want a stricter threshold
sig_labels <- ifelse(p_matrix < sig_level, "*", "")

# Now we create a heatmap of rg_matrix, and supply `sig_labels` to display_numbers.
pheatmap(
  rg_matrix[1:3072, 1:3072],                
  cluster_rows    = TRUE,                 
  cluster_cols    = TRUE,                 
  display_numbers = sig_labels[1:3072, 1:3072],  
  number_color    = "black",
  main            = "Genetic Correlation Estimates (GPT)",
  legend          = TRUE,
  # Set tree heights to zero to hide dendrogram lines
  treeheight_row  = 0,
  treeheight_col  = 0
)

# width: 1800 height:1800