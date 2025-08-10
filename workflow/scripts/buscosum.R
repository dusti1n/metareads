#!/usr/bin/env Rscript

# Load required libraries
library(ggplot2)
library(reshape2)
library(tidyverse)
library(gridExtra)
library(conflicted)

# Resolve function conflicts
conflict_prefer("filter", "dplyr")
conflict_prefer("lag", "dplyr")

# Read command-line arguments
args <- commandArgs(trailingOnly = TRUE)
# all input CSVs
input_csvs <- args[1:(length(args)-1)]
# last argument is output PDF
output_pdf <- args[length(args)]     

# Store ggplot objects
plot_list <- list()

# Loop through each CSV and create a BUSCO bar plot
for (csv_file in input_csvs) {
  busco_df <- read.csv(csv_file, header = TRUE)
  busco_df$Strain <- as.factor(busco_df$Strain)

  # Convert to long format
  busco_df.melted <- melt(busco_df, id.vars = "Strain")
  busco_df.melted$variable <- relevel(busco_df.melted$variable, "Missing")

  # Generate ggplot for each input
  p <- ggplot(busco_df.melted, aes(x = Strain, fill = fct_rev(variable), y = value)) +
    geom_bar(position = "stack", width = 0.7, stat = "identity") +
    labs(x = "Strain", y = "BUSCO", fill = "Type") +
    scale_y_continuous(labels = scales::comma) +
    theme_bw() +
    ggtitle(str_remove(basename(csv_file), "_busco_sum.*\\.csv$")) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      axis.text.y = element_text(size = 10),
      axis.title = element_text(size = 11),
      plot.title = element_text(size = 10)
    )

  # Add plot to list
  plot_list[[length(plot_list) + 1]] <- p
}

# Arrange all plots into one PDF (2 plots per row)
pdf(output_pdf, width = 11, height = 8)
cols <- if (length(plot_list) == 1) 1 else 2
do.call("grid.arrange", c(plot_list, ncol = cols))
dev.off()
