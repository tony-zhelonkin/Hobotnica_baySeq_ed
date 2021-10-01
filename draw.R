library("BiocManager")
library("dplyr")
library("tximeta")
library("SummarizedExperiment")

gse <- readRDS(file = "data/matrix.rds")
head(assays(gse)$counts)
#assays$counts - matrix
#    abundance ?
#   length ?
#rowData -
#rowRanges
#colData