# Logger config
library(logging)
basicConfig()
loginfo('Start')

# Load data
cm_args <- commandArgs(trailingOnly = TRUE)
count <- cm_args[1]
cols <- cm_args[2]
#option <- cm_args[4]
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

source("source/signatures_utils.R")

# Logger config
basicConfig()
loginfo('Make 5 signatures')
# Load differential expression analysis results
results_deseq2 <- readRDS(file = file.path(out, "DESeq_res.rds"))
results_ebseq <- readRDS(file = file.path(out, "EBSeq_res.rds"))
results_edger <- readRDS(file = file.path(out, "edgeR_res.rds"))
results_voom <- readRDS(file = file.path(out, "voom_res.rds"))
results_noiseq <- readRDS(file = file.path(out, "NOISeq_res.rds"))

top_signature(results_deseq2, results_ebseq, results_edger, results_voom, results_noiseq, 30, out)
#filtered_signature(results_deseq2, results_ebseq, results_edger, results_voom, results_noiseq, out)

loginfo('Visualize results')

# draw a Venn diagram
draw_venn_diag(out)