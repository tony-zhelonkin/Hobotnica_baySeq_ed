library(logging)
source("signatures_utils.R")

# Logger config
basicConfig()

loginfo('Make 5 signatures')
# Load differential expression analysis results
results_deseq2 <- readRDS(file = "data/DESeq_res.rds")
results_ebseq <- readRDS(file = "data/EBSeq_res.rds")
results_edger <- readRDS(file = "data/edgeR_res.rds")
results_voom <- readRDS(file = "data/voom_res.rds")

top_signature(results_deseq2, results_ebseq, results_edger, results_voom, 1000)
#filtered_signature(results_deseq2, results_ebseq, results_edger, results_voom)

loginfo('Visualize results')

# draw a Venn diagram
draw_venn_diag()