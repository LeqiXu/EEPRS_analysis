# Plot1:
## Steup
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(forcats)
library(data.table)

my_theme <- theme(
  # shrink axis titles
  axis.title.x = element_text(size = 6, face = "bold"),
  axis.title.y = element_text(size = 6, face = "bold"),
  # shrink tick labels
  axis.text.x  = element_text(size = 5, face = "bold"),
  axis.text.y  = element_text(size = 4, face = "bold"),
  # leave all other text at your preferred size
  plot.title   = element_text(size = 8, face = "bold"),
  strip.text   = element_text(size = 6, face = "bold"),
  panel.grid.major.x = element_line(color = "grey", size = 0.5),
  panel.grid.minor.x = element_line(color = "lightgrey", size = 0.25),
  panel.grid.major.y = element_blank(),
  panel.grid.minor.y = element_blank(),
  legend.position  = "top",
  legend.title     = element_text(size = 9, face = "bold"),     # title size
  legend.text      = element_text(size = 7),
  plot.margin = unit(c(0.1, 0.1, 0.5, 0.1), "cm")
)

method_colors <- c(
  "Word2Vec"     = "#8DA0CB",
  "Word2Vec_PCA" = "#E78AC3",
  "Word2Vec_ICA" = "#FC8D62",
  "GPT_PCA"      = "#66C2A5",
  "GPT_ICA"      = "#A6D854"
)

setwd("/gpfs/gibbs/pi/zhao/lx94/EEPRS/result/plot")

## Data prepare
prs_table = fread(paste0("/gpfs/gibbs/pi/zhao/lx94/EEPRS/result/evaluation/EEPRS_benchmarking_R2_AUC.csv"))
prs_table <- prs_table %>%
  mutate(trait = fct_reorder(trait, PRScs, .desc = TRUE))

prs_table_continuous = prs_table[which(prs_table$type == "R2"),c("pop","trait","PRScs","EEPRS_Word2Vec","EEPRS_Word2Vec_PCA","EEPRS_Word2Vec_ICA","EEPRS_GPT_PCA","EEPRS_GPT_ICA")]
long_table_continuous <- melt(prs_table_continuous, id.vars = c("pop", "trait"), variable.name = "method", value.name = "R2")
long_table_continuous$method = factor(long_table_continuous$method, levels = c("PRScs","EEPRS_Word2Vec","EEPRS_Word2Vec_PCA","EEPRS_Word2Vec_ICA","EEPRS_GPT_PCA","EEPRS_GPT_ICA"))

prs_table_binary = prs_table[which(prs_table$type == "AUC"),c("pop","trait","PRScs","EEPRS_Word2Vec","EEPRS_Word2Vec_PCA","EEPRS_Word2Vec_ICA","EEPRS_GPT_PCA","EEPRS_GPT_ICA")]
long_table_binary <- melt(prs_table_binary, id.vars = c("pop", "trait"), variable.name = "method", value.name = "AUC")
long_table_binary$method = factor(long_table_binary$method, levels = c("PRScs","EEPRS_Word2Vec","EEPRS_Word2Vec_PCA","EEPRS_Word2Vec_ICA","EEPRS_GPT_PCA","EEPRS_GPT_ICA"))

## Calcualte the relative change of other methods over PRScs
PRScs_ref_continuous <- long_table_continuous[long_table_continuous$method == "PRScs", c("pop", "trait", "R2")]
colnames(PRScs_ref_continuous)[3] <- "PRScs_auto_R2"
long_table_continuous_with_PRScs  <- merge(long_table_continuous, PRScs_ref_continuous, by = c("pop", "trait"))
long_table_continuous_with_PRScs$relative_change <- with(long_table_continuous_with_PRScs, (R2 - PRScs_auto_R2) / PRScs_auto_R2)
long_table_continuous_with_PRScs <- long_table_continuous_with_PRScs[long_table_continuous_with_PRScs$method != "PRScs",]
long_table_continuous_with_PRScs <- long_table_continuous_with_PRScs[,c("pop", "trait","method","relative_change")]
long_table_continuous_with_PRScs$method = factor(long_table_continuous_with_PRScs$method, levels = c("EEPRS_Word2Vec","EEPRS_Word2Vec_PCA","EEPRS_Word2Vec_ICA","EEPRS_GPT_PCA","EEPRS_GPT_ICA"))

PRScs_ref_binary <- long_table_binary[long_table_binary$method == "PRScs", c("pop", "trait", "AUC")]
colnames(PRScs_ref_binary)[3] <- "PRScs_auto_AUC"
PRScs_ref_binary <- PRScs_ref_binary[which(PRScs_ref_binary$PRScs_auto_AUC > 0.5),]
long_table_binary_with_PRScs  <- merge(long_table_binary, PRScs_ref_binary, by = c("pop", "trait"))
long_table_binary_with_PRScs$relative_change <- with(long_table_binary_with_PRScs, (AUC - PRScs_auto_AUC) / (PRScs_auto_AUC - 0.5))
long_table_binary_with_PRScs <- long_table_binary_with_PRScs[long_table_binary_with_PRScs$method != "PRScs",]
long_table_binary_with_PRScs <- long_table_binary_with_PRScs[,c("pop", "trait","method","relative_change")]
long_table_binary_with_PRScs$method = factor(long_table_binary_with_PRScs$method, levels = c("EEPRS_Word2Vec","EEPRS_Word2Vec_PCA","EEPRS_Word2Vec_ICA","EEPRS_GPT_PCA","EEPRS_GPT_ICA"))

## Plot
## continuous
# 1. Subset data
df_sub <- long_table_continuous_with_PRScs %>%
  filter(method %in% c("EEPRS_Word2Vec",
                       "EEPRS_Word2Vec_PCA",
                       "EEPRS_Word2Vec_ICA",
                       "EEPRS_GPT_PCA",
                       "EEPRS_GPT_ICA"))

# 2. Compute average improvement *by* method
df_avg <- df_sub %>%
  group_by(method) %>%
  summarize(avg_improve = mean(relative_change, na.rm = TRUE))

# 3. Single plot with facets
continuous_plot <- ggplot(df_sub, 
       aes(x = reorder(trait, relative_change),
           y = 100 * relative_change,
           fill = interaction(method, relative_change > 0))) +
  geom_col() +
  coord_flip() +
  scale_y_continuous(labels = function(x) paste0(x, "%")) +
  labs(
    x = "Trait",
    y = "Relative Improvement (%)",
    title = "Continuous Traits"
  ) +
  theme_classic() +
  my_theme +
  guides(fill = FALSE) +
  facet_wrap(~ method, ncol = 1, scales = "free_y") +
  # 4. Add average line & text for each method:
  geom_hline(
    data = df_avg,
    aes(yintercept = 100 * avg_improve),
    color = "purple",
    linetype = "dashed",
    size = 0.5
  ) +
  geom_text(
    data = df_avg,
    aes(
      x = 1,
      y = 100 * avg_improve,
      label = paste0("Avg improve: ", sprintf("%.2f%%", 100 * avg_improve))
    ),
    inherit.aes = FALSE,
    hjust = -0.1,
    vjust = -0.3,
    fontface = "bold",
    color = "purple",
    size = 2
  ) +
  scale_fill_manual(
    values = c(
      # Negative changes (all grey)
      "EEPRS_Word2Vec.FALSE" = "grey",
      "EEPRS_Word2Vec_PCA.FALSE" = "grey",
      "EEPRS_Word2Vec_ICA.FALSE" = "grey",
      "EEPRS_GPT_PCA.FALSE" = "grey",
      "EEPRS_GPT_ICA.FALSE" = "grey",
      # Positive changes (method-specific colors)
      "EEPRS_Word2Vec.TRUE"     = "#8DA0CB",
      "EEPRS_Word2Vec_PCA.TRUE" = "#E78AC3",
      "EEPRS_Word2Vec_ICA.TRUE" = "#FC8D62",
      "EEPRS_GPT_PCA.TRUE"      = "#66C2A5",
      "EEPRS_GPT_ICA.TRUE"      = "#A6D854"
    )
  )

print(continuous_plot)

## binary
# 1. Subset data
df_sub <- long_table_binary_with_PRScs %>%
  filter(method %in% c("EEPRS_Word2Vec",
                       "EEPRS_Word2Vec_PCA",
                       "EEPRS_Word2Vec_ICA",
                       "EEPRS_GPT_PCA",
                       "EEPRS_GPT_ICA"))

# 2. Compute average improvement *by* method
df_avg <- df_sub %>%
  group_by(method) %>%
  summarize(avg_improve = mean(relative_change, na.rm = TRUE))

# 3. Single plot with facets
binary_plot <- ggplot(df_sub, 
       aes(x = reorder(trait, relative_change),
           y = 100 * relative_change,
           fill = interaction(method, relative_change > 0))) +
  geom_col() +
  coord_flip() +
  scale_y_continuous(labels = function(x) paste0(x, "%")) +
  labs(
    x = "Trait",
    y = "Relative Improvement (%)",
    title = "Binary Traits"
  ) +
  theme_classic() +
  my_theme +
  guides(fill = FALSE) +
  facet_wrap(~ method, ncol = 1, scales = "free_y") +
  geom_hline(
    data = df_avg,
    aes(yintercept = 100 * avg_improve),
    color = "purple",
    linetype = "dashed",
    size = 0.5
  ) +
  geom_text(
    data = df_avg,
    aes(
      x = 1,
      y = 100 * avg_improve,
      label = paste0("Avg improve: ", sprintf("%.2f%%", 100 * avg_improve))
    ),
    inherit.aes = FALSE,
    hjust = -0.1,
    vjust = -0.3,
    fontface = "bold",
    color = "purple",
    size = 2
  ) +
  scale_fill_manual(
    values = c(
      # Negative changes (all grey)
      "EEPRS_Word2Vec.FALSE" = "grey",
      "EEPRS_Word2Vec_PCA.FALSE" = "grey",
      "EEPRS_Word2Vec_ICA.FALSE" = "grey",
      "EEPRS_GPT_PCA.FALSE" = "grey",
      "EEPRS_GPT_ICA.FALSE" = "grey",
      # Positive changes (method-specific colors)
      "EEPRS_Word2Vec.TRUE"     = "#8DA0CB",
      "EEPRS_Word2Vec_PCA.TRUE" = "#E78AC3",
      "EEPRS_Word2Vec_ICA.TRUE" = "#FC8D62",
      "EEPRS_GPT_PCA.TRUE"      = "#66C2A5",
      "EEPRS_GPT_ICA.TRUE"      = "#A6D854"
    )
  )

print(binary_plot)

combined_plot <- ggarrange(
  continuous_plot, binary_plot,
  ncol = 2,
  labels = c("A", "B"),
  font.label = list(size = 9)
)

print(combined_plot)

# Save the figure to match the journal guidelines
ggsave(filename = "Figure5.pdf",
       plot = combined_plot,
       width = 6.5, 
       height = 8.5,
       units = "in",
       dpi = 300,
       device   = cairo_pdf)  # high-resolution for publication