# Logger config
library(logging)
basicConfig()
loginfo('Start')

# Load data
cm_args <- commandArgs(trailingOnly = TRUE)
count <- cm_args[1]
cols <- cm_args[2]
if (!file.exists(count)) {
  stop("File ", count, " does not exists!")
}
if (!file.exists(cols)) {
  stop("File ", cols, " does not exists!")
}

# Read count matrix table and correct its form
counts<- read.table(file=count, header = TRUE, sep = ",", dec = ".")
gene_names <- counts$X
counts <- counts[,-c(1)]
rownames(counts) <- gene_names

# Read annotation and make table with useful form
coldata <- read.table(file=cols, sep = ",", dec = ".")
colnames(coldata) <- c("names", "condition")
coldata <- coldata[-c(1), ]

library(BiocManager)
library(dplyr)
library(tximeta)
library(SummarizedExperiment)
library(NOISeq)

# Prepare data
factors <- as.data.frame(coldata$condition)

noiseq_data <- readData(data = counts, factors = factors)

head(noiseq_data)
# Calculate diff expression
#noiseq_res <- noiseq(noiseq_data, factor = condition)