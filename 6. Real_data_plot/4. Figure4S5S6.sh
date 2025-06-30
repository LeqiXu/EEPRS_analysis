## Step1 plot the PheWAS Manhatten
library(data.table)
library(vroom)
library(dplyr)
library(stringr)
library(PheWAS)
library(ggplot2)
library(ggpubr)

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
  panel.grid.major.x = element_line(color = "grey", size = 0.5),
  panel.grid.minor.x = element_line(color = "lightgrey", size = 0.25),
  panel.grid.major.y = element_blank(),
  panel.grid.minor.y = element_blank(),
  legend.position  = "top",
  legend.title     = element_text(size = 9, face = "bold"),     # title size
  legend.text      = element_text(size = 7),
  plot.margin = unit(c(0.1, 0.1, 0.5, 0.1), "cm")
)

setwd("/gpfs/gibbs/pi/zhao/lx94/EEPRS/result/plot")

make_phewas_plot <- function(type, i) {
  title_type <- ifelse(type == "word2vec100_ICA30", "Word2Vec_ICA: ICA", "GPT_ICA: ICA")
  component <- i
  this_label <- paste0(title_type, "_", component)

  results <- fread(paste0("/gpfs/gibbs/pi/zhao/lx94/EEPRS/result/Embedding_PRS/PheWAS/UKB_", type, "_V", component, "_PheWAS.csv"))

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

  cutoff_p <- 0.05
  
  p <- phewasManhattan(
  results,
  annotate.phenotype.description = new.annotate.phenotype.description,
  annotate.angle = 0,
  title = this_label,
  point.size = 0.5,
  annotate.size = 2,
  OR.direction = TRUE,
  switch.axis = TRUE,
  significant.line = cutoff_p,
  suggestive.line = cutoff_p,
  annotate.level = cutoff_p,
  x.axis.label = "PheWAS Group"
  ) +
  theme_classic() + my_theme +
  scale_x_continuous(
    breaks = c(0, 5, 10, 15, cap_value),
    labels = c("0", "5", "10", "15", paste0(">", cap_value))
  ) +  
  scale_shape_manual(
    name = "Direction",
    values = c("FALSE" = 25, "TRUE" = 24),
    labels = c("Negative", "Positive")
  ) +
  theme(
    strip.text = element_text(size = 7, face = "bold"),
    legend.position = "right"
  )  +
  labs(
    x = expression(-log[10](BH[p]))
  )

return(p)


}

top_components <- list(
  word2vec100_ICA30 = c(25, 19, 12, 2, 15),
  gpticd3072_ICA67  = c(54, 1, 32, 15, 22)
)

# combined_plot1
top = 1

idx_1 = top_components[["word2vec100_ICA30"]][top]
p1 <- make_phewas_plot("word2vec100_ICA30", idx_1)

idx_2 = top_components[["gpticd3072_ICA67"]][top]
p2 <- make_phewas_plot("gpticd3072_ICA67", idx_2)

combined_plot1 <- ggarrange(p1, p2,
                           ncol = 1, nrow = 2,
                           heights = c(1, 1))     # equal height
print(combined_plot1)

ggsave(filename = "Figure4.pdf",
       plot = combined_plot1,
       width = 6.5, 
       height = 8,
       units = "in",
       dpi = 300,
       device   = cairo_pdf)  # high-resolution for publication

# combined_plot2
top = 2

idx_1 = top_components[["word2vec100_ICA30"]][top]
p1 <- make_phewas_plot("word2vec100_ICA30", idx_1)

idx_2 = top_components[["gpticd3072_ICA67"]][top]
p2 <- make_phewas_plot("gpticd3072_ICA67", idx_2)

combined_plot2 <- ggarrange(p1, p2,
                           ncol = 1, nrow = 2,
                           heights = c(1, 1))     # equal height
print(combined_plot2)

ggsave(filename = "FigureS5.pdf",
       plot = combined_plot2,
       width = 6.5, 
       height = 8,
       units = "in",
       dpi = 300,
       device   = cairo_pdf)  # high-resolution for publication

# combined_plot3
top = 3

idx_1 = top_components[["word2vec100_ICA30"]][top]
p1 <- make_phewas_plot("word2vec100_ICA30", idx_1)

idx_2 = top_components[["gpticd3072_ICA67"]][top]
p2 <- make_phewas_plot("gpticd3072_ICA67", idx_2)

combined_plot3 <- ggarrange(p1, p2,
                           ncol = 1, nrow = 2,
                           heights = c(1, 1))     # equal height
print(combined_plot3)

ggsave(filename = "FigureS6.pdf",
       plot = combined_plot3,
       width = 6.5, 
       height = 8,
       units = "in",
       dpi = 300,
       device   = cairo_pdf)  # high-resolution for publication