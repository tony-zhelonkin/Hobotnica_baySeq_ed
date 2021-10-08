library(logging)
source("signature_making_utils.R")

# Logger config
basicConfig()

loginfo('Make 5 signatures')
# Load differential expression analysis results
results_deseq2 <- readRDS(file = "data/DESeq_res.rds")
results_ebseq <- readRDS(file = "data/EBSeq_res.rds")
results_edger <- readRDS(file = "data/edgeR_res.rds")
results_voom <- readRDS(file = "data/voom_res.rds")

top_signature(results_deseq2, results_ebseq, results_edger, results_voom)

loginfo('Visualize results')
# read signatures from files
sig_edgeR <- as.vector(unlist(read.delim(file = "data/edgeR_sig.txt", header = FALSE)))
sig_DeSeq2 <- as.vector(unlist(read.delim(file = "data/DESeq_sig.txt", header = FALSE)))
sig_VoomLimma <- as.vector(unlist(read.delim(file = "data/voom_sig.txt", header = FALSE)))
sig_EBSeq <- as.vector(unlist(read.delim(file = "data/EBSeq_sig.txt", header = FALSE)))

# draw a Venn diagram
draw_venn_diag(sig_edgeR, sig_DeSeq2, sig_VoomLimma, sig_EBSeq)