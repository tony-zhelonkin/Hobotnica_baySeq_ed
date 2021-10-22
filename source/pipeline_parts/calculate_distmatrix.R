# Logger config
library(logging)
basicConfig()
loginfo('Start')
sink(type="output")
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
loginfo('Calculate dist matrixes for tools')
# Calculate dist matrices for tools
source("source/calculate_distmatrix_utils.R")
deseq_distmat <- calculate_distmatrix (count, file.path(out, "DESeq_sig.txt"))
ebseq_distmat <- calculate_distmatrix (count, file.path(out, "EBSeq_sig.txt"))
edger_distmat <- calculate_distmatrix (count, file.path(out, "edgeR_sig.txt"))
voom_distmat <- calculate_distmatrix (count, file.path(out, "voom_sig.txt"))
noiseq_distmat <- calculate_distmatrix (count, file.path(out, "NOISeq_sig.txt"))

# Save results of distance matrices
loginfo('Save results')
write.csv(as.data.frame(as.matrix(deseq_distmat)), file="DESeq_sig.txt.distmatrix")
write.csv(as.data.frame(as.matrix(ebseq_distmat)), file="EBSeq_sig.txt.distmatrix")
write.csv(as.data.frame(as.matrix(edger_distmat)), file="edgeR_sig.txt.distmatrix")
write.csv(as.data.frame(as.matrix(voom_distmat)), file="voom_sig.txt.distmatrix")
write.csv(as.data.frame(as.matrix(noiseq_distmat)), file="NOISeq_sig.txt.distmatrix")


