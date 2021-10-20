#read annotation
cm_args <- commandArgs(trailingOnly = TRUE)
cols <- cm_args[1]
if (!file.exists(cols)) {
  stop("File ", cols, " does not exists!")
}
coldata <- read.table(file=cols, sep = ",", dec = ".")
colnames(coldata) <- c("names", "condition")
coldata <- coldata[-c(1), ]

#load dist matrixes

distmat <- read.table(file="data/DESeq_sig.txt.distmatrix", header = TRUE, sep = ",", dec = ".")

source("hobotnica/R/Hobotnica.R")

cat(Hobotnica(distmat, coldata$condition))