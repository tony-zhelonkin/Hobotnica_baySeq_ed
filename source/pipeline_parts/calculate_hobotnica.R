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

if (length(cm_args) <= 2) {
    stop("Enter output directory name")
}
out <- cm_args[3]
if (!file.exists(out)) {
    dir.create(out)
}
deseq_distmat <- readRDS(file = file.path(out, "DESeq_sig.txt.distmatrix"))
ebseq_distmat <- readRDS(file = file.path(out, "EBSeq_sig.txt.distmatrix"))
edger_distmat <- readRDS(file = file.path(out, "edgeR_sig.txt.distmatrix"))
voom_distmat <- readRDS(file = file.path(out, "voom_sig.txt.distmatrix"))
noiseq_distmat <- readRDS(file = file.path(out, "NOISeq_sig.txt.distmatrix"))
loginfo('Start scoring Hobotnica')
source("source/hobotnica_using.R")
use_hobotnica(deseq_distmat, ebseq_distmat, edger_distmat, voom_distmat, noiseq_distmat, condition, out)
loginfo('Hobotnica scores are ready')