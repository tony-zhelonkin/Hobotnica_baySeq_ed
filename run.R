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

loginfo('Run 5 differential analysis tools and score with Hobotnica')
loginfo(paste0('Count matrix is in ', count))
loginfo(paste0('Coldata is in ', cols))

# Run tools
# Import functions from modules
source("DESeq.R")
source("EBSeq.R")
source("edgeR.R")
source("voom.R")
deseq2_res <- deseq2_f(counts, coldata)
ebseq_res <- ebseq_f(counts, coldata)
edger_res <- edger_f(counts, coldata)
voom_res <- voom_f(counts, coldata)


# Save results
loginfo('Save results')
saveRDS(deseq2_res, file = "data/DESeq_res.rds")
saveRDS(ebseq_res, file = "data/EBSeq_res.rds")
saveRDS(edger_res, file = "data/edgeR_res.rds")
saveRDS(voom_res, file = "data/voom_res.rds")

# Visualize results of differential expression
loginfo('Visualize results of differential expression')
deseq2_v(deseq2_res)
ebseq_v(ebseq_res)
edger_v(edger_res)
voom_v(voom_res)

loginfo('Make signatures of differential expression analysis')
source("signature_making_utils.R")
# Load differential expression analysis results
results_deseq2 <- readRDS(file = "data/DESeq_res.rds")
results_ebseq <- readRDS(file = "data/EBSeq_res.rds")
results_edger <- readRDS(file = "data/edgeR_res.rds")
results_voom <- readRDS(file = "data/voom_res.rds")

top_signature(results_deseq2, results_ebseq, results_edger, results_voom)

loginfo('Visualize signature crossing')
# Read signatures from files
sig_edgeR <- as.vector(unlist(read.delim(file = "data/edgeR_sig.txt", header = FALSE)))
sig_DeSeq2 <- as.vector(unlist(read.delim(file = "data/DESeq_sig.txt", header = FALSE)))
sig_VoomLimma <- as.vector(unlist(read.delim(file = "data/voom_sig.txt", header = FALSE)))
sig_EBSeq <- as.vector(unlist(read.delim(file = "data/EBSeq_sig.txt", header = FALSE)))

# Draw a Venn diagram
draw_venn_diag(sig_edgeR, sig_DeSeq2, sig_VoomLimma, sig_EBSeq)

loginfo('Calculate dist matrixes for tools')
# Calculate dist matrixes for tools
source("calculate_distmatrix_utils.R")
calculate_distmatrixes(count)

loginfo('Ready for Hobotnica using')