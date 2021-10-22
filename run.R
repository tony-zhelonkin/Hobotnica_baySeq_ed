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
save_flag <- 0
if (length(cm_args)  >= 4) {
    option <- cm_args[4]
    if (option == "-s") {
        save_flag <- 1
    }
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

loginfo('Run 5 differential analysis tools and score with Hobotnica')
loginfo(paste0('Count matrix is in ', count))
loginfo(paste0('Coldata is in ', cols))
loginfo(paste0('Output will be written in ', out))

# Run tools
# Import functions from modules
source("source/DESeq.R")
source("source/EBSeq.R")
source("source/edgeR.R")
source("source/voom.R")
source("source/NOISeq.R")
deseq2_res <- deseq2_f(counts, coldata)
ebseq_res <- ebseq_f(counts, coldata)
edger_res <- edger_f(counts, coldata)
voom_res <- voom_f(counts, coldata)
noiseq_res <- noiseq_f(counts, coldata)

if (save_flag == 1) {
    # Save results of differential expression
    loginfo('Save results of differential expression')
    saveRDS(deseq2_res, file = file.path(out, "DESeq_res.rds"))
    saveRDS(ebseq_res, file = file.path(out, "EBSeq_res.rds"))
    saveRDS(edger_res, file = file.path(out, "edgeR_res.rds"))
    saveRDS(voom_res, file = file.path(out, "voom_res.rds"))
    saveRDS(noiseq_res, file = file.path(out, "NOISeq_res.rds"))
}


# Visualize results of differential expression
loginfo('Visualize results of differential expression')
deseq2_v(deseq2_res, out)
ebseq_v(ebseq_res, out)
edger_v(edger_res, out)
voom_v(voom_res, out)
noiseq_v(noiseq_res, out)

loginfo('Make signatures of differential expression analysis')
source("source/signatures_utils.R")

top_signature(deseq2_res, ebseq_res, edger_res, voom_res, noiseq_res, 30, out)

loginfo('Visualize signature crossing')
# Read signatures from files

# Draw a Venn diagram
draw_venn_diag(out)

loginfo('Calculate dist matrixes for tools')
# Calculate dist matrices for tools
source("source/calculate_distmatrix_utils.R")
deseq_distmat <- calculate_distmatrix (count, file.path(out, "DESeq_sig.txt"))
ebseq_distmat <- calculate_distmatrix (count, file.path(out, "EBSeq_sig.txt"))
edger_distmat <- calculate_distmatrix (count, file.path(out, "edgeR_sig.txt"))
voom_distmat <- calculate_distmatrix (count, file.path(out, "voom_sig.txt"))
noiseq_distmat <- calculate_distmatrix (count, file.path(out, "NOISeq_sig.txt"))

if (save_flag == 1) {
    # Save results of distance matrices
    loginfo('Save results')
    write.csv(as.data.frame(as.matrix(deseq_distmat)), file=file.path(out,"DESeq_sig.txt.distmatrix"))
    write.csv(as.data.frame(as.matrix(ebseq_distmat)), file=file.path(out,"EBSeq_sig.txt.distmatrix"))
    write.csv(as.data.frame(as.matrix(edger_distmat)), file=file.path(out,"edgeR_sig.txt.distmatrix"))
    write.csv(as.data.frame(as.matrix(voom_distmat)), file=file.path(out,"voom_sig.txt.distmatrix"))
    write.csv(as.data.frame(as.matrix(noiseq_distmat)), file=file.path(out,"NOISeq_sig.txt.distmatrix"))
}


loginfo('Ready for Hobotnica using')
loginfo('Start scoring Hobotnica')
source("source/hobotnica_using.R")
use_hobotnica(deseq_distmat, ebseq_distmat, edger_distmat, voom_distmat, noiseq_distmat, coldata$condition, out)
loginfo('Hobotnica scores are ready')
