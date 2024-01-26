#!/usr/bin/env Rscript
library("tidyverse")

args <- commandArgs(trailingOnly = TRUE)

# Check if arguments are provided
if(length(args)==0) {
  stop("No arguments provided.")
}

# Parse arguments
input_file <- args[1]
output_file <- args[2]

df = read.csv(input_file, header = TRUE, sep = ',', stringsAsFactors = FALSE)

df = df %>%
    group_by(annotation_id, gene, ko, ontology) %>%
    summarise(count = n())

write.csv(df, output_file, row.names = FALSE)