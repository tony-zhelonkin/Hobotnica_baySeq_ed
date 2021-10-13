source("calculate_distmatrix_utils.R")
args = commandArgs(trailingOnly=TRUE)
countMatrixFile <- args[1]
if (!file.exists(countMatrixFile)) {
  stop("File ", countMatrixFile, " does not exists!")
}
calculate_distmatrixes(countMatrixFile)

