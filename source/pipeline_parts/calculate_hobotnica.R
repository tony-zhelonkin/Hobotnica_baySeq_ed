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

coldata <- read.table(file=cols, sep = ",", dec = ".")
colnames(coldata) <- c("names", "condition")
coldata <- coldata[-c(1), ]

deseq_distmat <- read.table(file =file.path(out, "DESeq_sig.txt.distmatrix"), header=TRUE, sep=",")
deseq_distmat <- as.matrix(data.frame(deseq_distmat[, -1], row.names = deseq_distmat[, 1]))
colnames(deseq_distmat) <- rownames(deseq_distmat)

ebseq_distmat <- read.table(file =file.path(out, "EBSeq_sig.txt.distmatrix"), header=TRUE, sep=",")
ebseq_distmat <- as.matrix(data.frame(ebseq_distmat[, -1], row.names = ebseq_distmat[, 1]))
colnames(ebseq_distmat) <- rownames(ebseq_distmat)

edger_distmat <- read.table(file =file.path(out, "edgeR_sig.txt.distmatrix"), header=TRUE, sep=",")
edger_distmat <- as.matrix(data.frame(edger_distmat[, -1], row.names = edger_distmat[, 1]))
colnames(edger_distmat) <- rownames(edger_distmat)

voom_distmat <- read.table(file =file.path(out, "voom_sig.txt.distmatrix"), header=TRUE, sep=",")
voom_distmat <- as.matrix(data.frame(voom_distmat[, -1], row.names = voom_distmat[, 1]))
colnames(voom_distmat) <- rownames(voom_distmat)

noiseq_distmat <- read.table(file =file.path(out, "NOISeq_sig.txt.distmatrix"), header=TRUE, sep=",")
noiseq_distmat <- as.matrix(data.frame(noiseq_distmat[, -1], row.names = noiseq_distmat[, 1]))
colnames(noiseq_distmat) <- rownames(noiseq_distmat)

loginfo('Start scoring Hobotnica')
source("source/hobotnica_using.R")
use_hobotnica(deseq_distmat, ebseq_distmat, edger_distmat, voom_distmat, noiseq_distmat, coldata$condition, out)
loginfo('Hobotnica scores are ready')