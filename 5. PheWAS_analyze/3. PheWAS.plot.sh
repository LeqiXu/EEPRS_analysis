## Step1 plot the PheWAS Manhatten
library(data.table)
library(vroom)
library(dplyr)
library(stringr)
library(PheWAS)
library(ggplot2)
library(ggpubr)

make_phewas_plot <- function(type, i) {
  title_type <- ifelse(type == "word2vec100_ICA30", "Word2Vec_ICA", "GPT_ICA")
  component <- paste0("V", i)
  this_label <- paste0(title_type, "_", component)

  results <- fread(paste0("/gpfs/gibbs/pi/zhao/lx94/EEPRS/result/Embedding_PRS/PheWAS/UKB_", type, "_", component, "_PheWAS.csv"))

  new.annotate.phenotype.description <- vroom::vroom(
    "/gpfs/gibbs/pi/zhao/lx94/EEPRS/data/phecode/original_phecodes_pheinfo.csv",
    .name = janitor::make_clean_names,
    delim = ",",
    col_types = c(phecode = "c", description = "c", groupnum = "double", group = "c", color = "c")
  )
  colnames(new.annotate.phenotype.description)[1] <- "phenotype"

  results <- results %>%
    mutate(p = p.adjust(p, method = "BH"))
  
  min_p <- 1e-20
  cap_value <- -log10(min_p)
  results$p[results$p < min_p] <- min_p
  colnames(results)[1] <- "phenotype"
  results$phenotype <- as.character(results$phenotype)

  results <- results %>%
    mutate(
      logp = -log10(p),
      logp_capped = ifelse(logp > cap_value, cap_value, logp),
      prs_label = this_label
    )

  my_theme <- theme(
    axis.title.x = element_text(size = 8, face = "bold"),
    axis.title.y = element_text(size = 8, face = "bold"),
    axis.text.x  = element_text(angle = 45, hjust = 1, size = 5, face = "bold"),
    axis.text.y  = element_text(size = 6, face = "bold"),
    plot.title   = element_text(size = 10, face = "bold"),
    strip.text   = element_text(size = 10, face = "bold"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_line(color = "grey", size = 0.5),
    panel.grid.minor.y = element_line(color = "lightgrey", size = 0.25),
    legend.position  = "top",
    legend.title     = element_text(size = 9, face = "bold"),
    legend.text      = element_text(size = 7),
    plot.margin = unit(c(0.1, 0.1, 0.5, 0.1), "cm")
  )

  cutoff_p <- 0.05

  p <- phewasManhattan(
    results,
    annotate.phenotype.description = new.annotate.phenotype.description,
    annotate.angle = 0,
    title = this_label,
    point.size = 0.5,
    annotate.size = 2.5,
    significant.line = cutoff_p,
    suggestive.line = cutoff_p,
    annotate.level = cutoff_p
  ) +
    theme_classic() + my_theme +
    scale_y_continuous(
      breaks = c(0, 5, 10, 15, cap_value),
      labels = c("0", "5", "10", "15", paste0(">", cap_value))
    )

  return(p)
}

top_components <- list(
  word2vec100_ICA30 = c(25, 19, 12, 2, 15),
  gpticd3072_ICA67  = c(54, 1, 32, 15, 22)
)

# Generate the plots
top = 3

idx_1 = top_components[["word2vec100_ICA30"]][top]
p1 <- make_phewas_plot("word2vec100_ICA30", idx_1)

idx_2 = top_components[["gpticd3072_ICA67"]][top]
p2 <- make_phewas_plot("gpticd3072_ICA67", idx_2)

# Arrange vertically using ggarrange
combined_plot <- ggarrange(p1, p2,
                           ncol = 2, nrow = 1,
                           labels = c("A", "B"),  # optional labels
                           heights = c(1, 1))     # equal height

# Display
print(combined_plot)

## Step2 summarize the significant disease for each method
library(data.table)
library(dplyr)
library(ggplot2)
library(vroom)

#---------------------------------------------
# 1) Define a function to process a single type
#---------------------------------------------
process_type <- function(type) {
  if (type == "word2vec100_ICA30") {
    sig_embedding <- c(25, 19, 12)
    title_type    <- "Word2Vec_ICA"
  } else if (type == "gpticd3072_ICA67") {
    sig_embedding <- c(54, 1, 32)
    title_type    <- "GPT_ICA"
  } else {
    stop("Unknown type: ", type)
  }
  
  # 1a) Load PheWAS annotation
  new.annotate.phenotype.description <- vroom::vroom(
    "/gpfs/gibbs/pi/zhao/lx94/EEPRS/data/phecode/original_phecodes_pheinfo.csv",
    .name = janitor::make_clean_names,
    delim = ",",
    col_types = c(phecode = "c", description = "c", groupnum = "double", group = "c", color = "c")
  )
  colnames(new.annotate.phenotype.description)[1] <- "phenotype"

  # 1b) Read all embedding results, combine
  all_results_list <- list()
  
  for (i in sig_embedding) {
    component  <- paste0("V", i)
    this_label <- paste0(title_type, "_", component)
    file_in    <- paste0("/gpfs/gibbs/pi/zhao/lx94/EEPRS/result/Embedding_PRS/PheWAS/UKB_",
                         type, "_", component, "_PheWAS.csv")
    
    if (!file.exists(file_in)) {
      message("File not found: ", file_in, " -- skipping")
      next
    }
    
    tmp <- fread(file_in)

    # Ensure a "phenotype" column
    colnames(tmp)[1] <- "phenotype"
    tmp$phenotype <- as.character(tmp$phenotype)
    
    # BH adjust
    tmp <- tmp %>%
      mutate(p = p.adjust(p, method = "BH"))

    # Clean up p-values
    min_p     <- 1e-20
    cap_value <- -log10(min_p)
    tmp$p[tmp$p < min_p] <- min_p

    tmp <- tmp %>%
    mutate(
      logp = -log10(p),
      logp_capped = ifelse(logp > cap_value, cap_value, logp),
      prs_label = this_label
    )
    
    all_results_list[[this_label]] <- tmp
  }
  
  # Combine into one data frame
  all_results <- dplyr::bind_rows(all_results_list)
  
  #---------------------------------------
  # 2) Join with phenotype annotation
  #---------------------------------------
  all_results_joined <- all_results %>%
    select(phenotype, prs_label, predictor, p, logp_capped) %>%
    left_join(new.annotate.phenotype.description, by = "phenotype") %>%
    na.omit()
  
  all_results_joined$title_type <- title_type
  
  return(all_results_joined)
}

#-----------------------------------------------------------------
# 2) Loop over both types, combine them, then facet by type/label
#-----------------------------------------------------------------
types <- c("gpticd3072_ICA67", "word2vec100_ICA30")

big_list <- list()
for (t in types) {
  cat("Processing type:", t, "...\n")
  out_df <- process_type(t)
  big_list[[t]] <- out_df
}

final_df <- dplyr::bind_rows(big_list)
final_df$title_type <- factor(final_df$title_type, levels = c("Word2Vec_ICA","GPT_ICA"))

#-----------------------------------------------------------------
# 3) Find significant pehnotypes
#-----------------------------------------------------------------
sig_df = final_df[which(final_df$p < 0.05),]
sig_df_Word2Vec_ICA = sig_df[which(sig_df$title_type == "Word2Vec_ICA"),c("prs_label","p","description","group")]
sig_df_GPT_ICA = sig_df[which(sig_df$title_type == "GPT_ICA"),c("prs_label","p","description","group")]