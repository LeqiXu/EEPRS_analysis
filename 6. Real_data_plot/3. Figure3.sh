library(data.table)
library(dplyr)
library(ggplot2)
library(vroom)
library(stringr)
library(RColorBrewer)

my_theme <- theme(
  # shrink axis titles
  axis.title.x = element_text(size = 8, face = "bold"),
  axis.title.y = element_text(size = 8, face = "bold"),
  # shrink tick labels
  axis.text.x  = element_text(size = 5, face = "bold"),
  axis.text.y  = element_text(size = 6, face = "bold"),
  # leave all other text at your preferred size
  plot.title   = element_text(size = 10, face = "bold"),
  strip.text   = element_text(size = 10, face = "bold"),
  panel.grid.major.x = element_blank(),
  panel.grid.minor.x = element_blank(),
  panel.grid.major.y = element_line(color = "grey", size = 0.5),
  panel.grid.minor.y = element_line(color = "lightgrey", size = 0.25),
  legend.position  = "top",
  legend.title     = element_text(size = 9, face = "bold"),     # title size
  legend.text      = element_text(size = 7),
  plot.margin = unit(c(0.1, 0.1, 0.5, 0.1), "cm")
)

setwd("/gpfs/gibbs/pi/zhao/lx94/EEPRS/result/plot")

#---------------------------------------------
# 1) Define a function to process a single type
#---------------------------------------------
process_type <- function(type) {
  if (type == "word2vec100_ICA30") {
    sig_embedding <- c(25,19,12,2,15,11,22,27,4,14,24,5,8,1,10,20,16,13,9,6,29,7)
    title_type    <- "Word2Vec_ICA"
  } else if (type == "gpticd3072_ICA67") {
    sig_embedding <- c(54,1,32,15,22,9,56,67,43,25,18,38,12,2,50,20,57,28,37,52,30,49,29,23,21,11,3,51,40,55,42,19,48,5,63,24,34,45,44,35,31,4,7,46)
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
  
  #---------------------------------------
  # 3) Summarize by (group, prs_label)
  #    max_logp for coloring
  #    min p for significance check
  #---------------------------------------
  group_summary_multi <- all_results_joined %>%
    group_by(group, prs_label) %>%
    summarize(
      max_logp = max(logp_capped, na.rm = TRUE),
      min_p = min(p, na.rm = TRUE),
      .groups  = "drop"
    ) %>%
    mutate(
      is_sig = min_p < 0.05,
      sig_label = ifelse(is_sig, "*", "")
    )
  
  #---------------------------------------
  # 4) Make sure x-axis is in the same order
  #---------------------------------------
  group_summary_multi <- group_summary_multi %>%
    mutate(prs_label = prs_label %>%
             str_replace("^Word2Vec_ICA_V", "ICA_") %>%
             str_replace("^GPT_ICA_V", "ICA_"))

  ordered_labels <- sapply(sig_embedding, function(ii) paste0("ICA_", ii))
  group_summary_multi$prs_label <- factor(
    group_summary_multi$prs_label,
    levels = ordered_labels
  )
  
  #---------------------------------------
  # 6) Add a column "title_type" or "type" for facet
  #---------------------------------------
  group_summary_multi$title_type <- title_type
  
  return(group_summary_multi)
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

#--------------------------------------------------------------
# 3) Single ggplot: facet by "title_type"
#--------------------------------------------------------------
combined_plot <- ggplot(final_df, aes(x = prs_label, y = group, fill = max_logp)) +
  geom_tile(color = "white", linewidth = 0.2) +
  geom_text(aes(label = sig_label), color = "black", size = 2, fontface = "bold") +
  scale_fill_gradient(
    low = "white",
    high = "#D73027",
    name = expression(-log[10](BH[p]))
  ) +
  labs(
    x = "Embedding PRS",
    y = "PheWAS Group"
  ) +
  facet_wrap(~ title_type, scales = "free", ncol = 1) +
  theme_classic(base_size = 12) +
  my_theme + 
  theme(
    axis.text.x  = element_text(size = 5, face = "bold",angle = 90, hjust = 1),
    strip.text = element_text(size = 7, face = "bold"),
    legend.position = "right"
  )

print(combined_plot)

# Save the figure to match the journal guidelines
ggsave(filename = "Figure3.pdf",
       plot = combined_plot,
       width = 6.5, 
       height = 8,
       units = "in",
       dpi = 300,
       device   = cairo_pdf)  # high-resolution for publication