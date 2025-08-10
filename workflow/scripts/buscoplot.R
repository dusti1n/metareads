#!/usr/bin/env Rscript

# Load libraries and set conflict preferences
library(conflicted)

# Load required packages
conflict_prefer("filter", "dplyr")
conflict_prefer("lag", "dplyr")

# Load required libraries for data manipulation and visualization
library(ggplot2)
library(reshape2)
library(tidyverse)

# Read command-line arguments
args <- commandArgs(trailingOnly = TRUE)
# Input BUSCO CSV file
input_csv <- args[1]
# Output PDF path
output_plot <- args[2]
# Sample name (currently unused)
wildcard_sample <- args[3] 

# Read BUSCO summary table
busco_df <- read.csv(input_csv, header = TRUE)
busco_df$Strain <- as.factor(busco_df$Strain)

# Reshape data to long format for ggplot
busco_df.melted <- melt(busco_df, id.vars = "Strain")
busco_df.melted$variable <- relevel(busco_df.melted$variable, "Missing")

# Generate stacked bar plot
busco_plot <- ggplot(busco_df.melted, aes(x = Strain, fill = fct_rev(variable), y = value)) +
  geom_bar(position = "stack", width = 0.7, stat = "identity") +
  labs(x = "Strain", y = "BUSCO", fill = "Type") +
  scale_y_continuous(labels = scales::comma) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 12))

# Save plot as PDF
pdf(output_plot, width = 8, height = 5, paper = 'special')
print(busco_plot)
dev.off()
