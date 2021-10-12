source("calculate_distmatrix_utils.R")
args = commandArgs(trailingOnly=TRUE)
countMatrixFile <- args[1]
calculate_distmatrixes(countMatrixFile)

