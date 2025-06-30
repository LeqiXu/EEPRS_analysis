# 1. Correlation calculation
## PCA
job_file="/gpfs/gibbs/pi/zhao/lx94/EEPRS/code/glm/gpticd3072_PCA67_corr.txt"
> $job_file  # Empty the job file if it already exists

train_type=train

for i in {1..66}; do
for j in $(seq $((i+1)) 67); do

# Construct the sumstats files for embeddings i and j (comma-separated)
file_i="/gpfs/gibbs/pi/zhao/lx94/EEPRS/data/embedding_data/munge_data/gpticd3072_PCA67_${train_type}_EUR_UKB_PC${i}_ldsc.sumstats.gz"
file_j="/gpfs/gibbs/pi/zhao/lx94/EEPRS/data/embedding_data/munge_data/gpticd3072_PCA67_${train_type}_EUR_UKB_PC${j}_ldsc.sumstats.gz"
    
out_prefix="/gpfs/gibbs/pi/zhao/lx94/EEPRS/result/h2_corr/corr/embedding_embedding_corr/gpticd3072_PCA67_${train_type}_EUR_UKB_PC${i}_${j}_rg"

if [[ ! -e ${out_prefix}.log ]]; then
echo "module load miniconda; conda activate ldsc; python /gpfs/gibbs/pi/zhao/lx94/Software/ldsc/ldsc.py --rg "${file_i},${file_j}" --ref-ld-chr /gpfs/gibbs/pi/zhao/lx94/Data/ldsc/eur_w_ld_chr/ --w-ld-chr /gpfs/gibbs/pi/zhao/lx94/Data/ldsc/eur_w_ld_chr/ --out "${out_prefix}"" >> $job_file
fi 

done
done

module load dSQ
cd /gpfs/gibbs/pi/zhao/lx94/EEPRS/code/glm/

dsq --job-file /gpfs/gibbs/pi/zhao/lx94/EEPRS/code/glm/gpticd3072_PCA67_corr.txt --partition=scavenge,day --requeue --mem=30G --cpus-per-task=1 --ntasks=1 --nodes=1 --time=1-00:00:00 --mail-type=ALL
sbatch dsq-gpticd3072_PCA67_corr-$(date +%Y-%m-%d).sh

## ICA
job_file="/gpfs/gibbs/pi/zhao/lx94/EEPRS/code/glm/gpticd3072_ICA67_corr.txt"
> $job_file  # Empty the job file if it already exists

train_type=train

for i in {1..66}; do
for j in $(seq $((i+1)) 67); do

# Construct the sumstats files for embeddings i and j (comma-separated)
file_i="/gpfs/gibbs/pi/zhao/lx94/EEPRS/data/embedding_data/munge_data/gpticd3072_ICA67_${train_type}_EUR_UKB_V${i}_ldsc.sumstats.gz"
file_j="/gpfs/gibbs/pi/zhao/lx94/EEPRS/data/embedding_data/munge_data/gpticd3072_ICA67_${train_type}_EUR_UKB_V${j}_ldsc.sumstats.gz"
    
out_prefix="/gpfs/gibbs/pi/zhao/lx94/EEPRS/result/h2_corr/corr/embedding_embedding_corr/gpticd3072_ICA67_${train_type}_EUR_UKB_V${i}_${j}_rg"

if [[ ! -e ${out_prefix}.log ]]; then
echo "module load miniconda; conda activate ldsc; python /gpfs/gibbs/pi/zhao/lx94/Software/ldsc/ldsc.py --rg "${file_i},${file_j}" --ref-ld-chr /gpfs/gibbs/pi/zhao/lx94/Data/ldsc/eur_w_ld_chr/ --w-ld-chr /gpfs/gibbs/pi/zhao/lx94/Data/ldsc/eur_w_ld_chr/ --out "${out_prefix}"" >> $job_file
fi
    
done
done

module load dSQ
cd /gpfs/gibbs/pi/zhao/lx94/EEPRS/code/glm/

dsq --job-file /gpfs/gibbs/pi/zhao/lx94/EEPRS/code/glm/gpticd3072_ICA67_corr.txt --partition=scavenge,day --requeue --mem=30G --cpus-per-task=1 --ntasks=1 --nodes=1 --time=1-00:00:00 --mail-type=ALL
sbatch dsq-gpticd3072_ICA67_corr-$(date +%Y-%m-%d).sh

# 2. Obtain correlation table
setwd("/gpfs/gibbs/pi/zhao/lx94/EEPRS/result/h2_corr/corr/embedding_embedding_corr/")

## Load libraries
library(pheatmap)

## 1. Define a helper to read rg/p from a single logfile
get_rg_p_from_logfile <- function(fname) {
  lines <- readLines(fname, warn = FALSE)
  hdr   <- grep("^p1\\s+p2\\s+rg\\s+se\\s+z\\s+p", lines)
  if (length(hdr) == 1 && hdr + 1 <= length(lines)) {
    cols <- strsplit(lines[hdr + 1], "\\s+")[[1]]
    return(c(as.numeric(cols[3]), as.numeric(cols[6])))  # rg, p
  }
  return(c(NA_real_, NA_real_))
}

## 2. Define a single function to parse correlation logs (for PCA or ICA)
read_correlation_logs <- function(pattern, prefix) {
  files <- list.files(pattern = pattern, full.names = TRUE)
  out   <- data.frame(label_i = character(), label_j = character(),
                      rg = numeric(), p = numeric(),
                      stringsAsFactors = FALSE)
  
  for (f in files) {
    # Extract i, j from filename using capturing groups
    # pattern must have (\\d+) capturing the two numbers
    # e.g. "gpticd3072_PCA67_train_EUR_UKB_PC(\\d+)_(\\d+)_rg\\.log"
    m <- regexec(pattern, basename(f))
    x <- regmatches(basename(f), m)[[1]]
    if (length(x) == 3) {
      i_val <- as.integer(x[2])
      j_val <- as.integer(x[3])
      rgp   <- get_rg_p_from_logfile(f)
      out   <- rbind(out, data.frame(
        label_i = paste0(prefix, "_", i_val),
        label_j = paste0(prefix, "_", j_val),
        rg      = rgp[1],
        p       = rgp[2],
        stringsAsFactors = FALSE
      ))
    }
  }
  return(out)
}

## 3. Read PCA and ICA logs, then row-bind
## Adjust your patterns as needed
pca_results <- read_correlation_logs(
  "^gpticd3072_PCA67_train_EUR_UKB_PC(\\d+)_(\\d+)_rg\\.log$",
  "PCA"
)
ica_results <- read_correlation_logs(
  "^gpticd3072_ICA67_train_EUR_UKB_V(\\d+)_(\\d+)_rg\\.log$",
  "ICA"
)

pca_results$embedding_approach = "GPT_PCA"
pca_results

write.table(pca_results, file=paste0("/gpfs/gibbs/pi/zhao/lx94/EEPRS/result/h2_corr/corr/embedding_embedding_corr/UKB_gpticd3072_pca_train_corr.txt"), row.names=F, col.names=T, quote=F, append=F, sep = "\t")

ica_results$embedding_approach = "GPT_ICA"
ica_results

write.table(ica_results, file=paste0("/gpfs/gibbs/pi/zhao/lx94/EEPRS/result/h2_corr/corr/embedding_embedding_corr/UKB_gpticd3072_ica_train_corr.txt"), row.names=F, col.names=T, quote=F, append=F, sep = "\t")

## Build correlation matrix for all unique labels
build_correlation_matrices <- function(corr_df) {
  all_labels <- unique(c(corr_df$label_i, corr_df$label_j))
  all_labels <- sort(all_labels)
  N <- length(all_labels)

  rg_mat <- matrix(NA_real_, nrow = N, ncol = N,
                   dimnames = list(all_labels, all_labels))
  p_mat  <- rg_mat

  for (k in seq_len(nrow(corr_df))) {
    i_idx <- match(corr_df$label_i[k], all_labels)
    j_idx <- match(corr_df$label_j[k], all_labels)
    rg_val <- corr_df$rg[k]
    p_val  <- corr_df$p[k]
    rg_mat[i_idx, j_idx] <- rg_val
    rg_mat[j_idx, i_idx] <- rg_val
    p_mat[i_idx, j_idx] <- p_val
    p_mat[j_idx, i_idx] <- p_val
  }

  diag(rg_mat) <- 1
  diag(p_mat)  <- 0  # or NA if preferred

  list(rg = rg_mat, p = p_mat)
}

pca_matrices <- build_correlation_matrices(pca_results)
ica_matrices <- build_correlation_matrices(ica_results)

pca_sig_embedding <- paste0("PCA_",c(2,3,4,12,6,1,9,15,14,5,22,21,38,34,10,18,11,20,16,42,40,32,19,24,23,33,26,37,17,7,44,55,56,8,43,54,45,29,48,39,25,35,65,46,51,31,28,67,13,41,30,60))
ica_sig_embedding <- paste0("ICA_",c(54,1,32,15,22,9,56,67,43,25,18,38,12,2,50,20,57,28,37,52,30,49,29,23,21,11,3,51,40,55,42,19,48,5,63,24,34,45,44,35,31,4,7,46))

pca_rg_matrix <- pca_matrices$rg[pca_sig_embedding,pca_sig_embedding]
pca_p_matrix <- pca_matrices$p[pca_sig_embedding,pca_sig_embedding]

ica_rg_matrix <- ica_matrices$rg[ica_sig_embedding,ica_sig_embedding]
ica_p_matrix <- ica_matrices$p[ica_sig_embedding,ica_sig_embedding]

## Plot the heatmap (mark significant p-values with "*")
library(pheatmap)
library(grid)
library(ggpubr)

sig_level <- 0.05

## PCA heatmap as a grob
pca_sig_labels <- ifelse(pca_p_matrix < sig_level, "*", "")
pca_plot <- grid::grid.grabExpr(
  pheatmap(
    pca_rg_matrix,
    cluster_rows    = TRUE,
    cluster_cols    = TRUE,
    display_numbers = pca_sig_labels,
    number_color    = "black",
    main            = "Genetic Correlation Estimates (GPT_PCA)",
    legend          = TRUE,
    treeheight_row  = 0,
    treeheight_col  = 0
  )
)

## ICA heatmap as a grob
ica_sig_labels <- ifelse(ica_p_matrix < sig_level, "*", "")
ica_plot <- grid::grid.grabExpr(
  pheatmap(
    ica_rg_matrix,
    cluster_rows    = TRUE,
    cluster_cols    = TRUE,
    display_numbers = ica_sig_labels,
    number_color    = "black",
    main            = "Genetic Correlation Estimates (GPT_ICA)",
    legend          = TRUE,
    treeheight_row  = 0,
    treeheight_col  = 0
  )
)

## Combine with ggarrange
ggarrange(
  ggplotify::as.ggplot(pca_plot),
  ggplotify::as.ggplot(ica_plot),
  ncol = 2,
  labels = c("A", "B")
)

## Width: 1600 Height: 800