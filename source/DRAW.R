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

# Read annotation and make table with useful form
coldata <- read.table(file=cols, sep = ",", dec = ".")
colnames(coldata) <- c("names", "condition")
coldata <- coldata[-c(1), ]
condition <- coldata$condition
n <- 30

# Import functions from modules
source("source/DESeq.R")
source("source/EBSeq.R")
source("source/edgeR.R")
source("source/voom.R")
source("source/NOISeq.R")

deseq2_res <- readRDS(file = file.path(out, "DESeq_res.rds"))
ebseq_res <- readRDS(file = file.path(out, "EBSeq_res.rds"))
edger_res <- readRDS(file = file.path(out, "edgeR_res.rds"))
voom_res <- readRDS(file = file.path(out, "voom_res.rds"))
noiseq_res <- readRDS(file = file.path(out, "NOISeq_res.rds"))

pdf(file.path(out, "heatmaps.pdf"))
# Visualize using heat maps
loginfo('Calculate heat maps for tools')
source("source/calculate_distmatrix_utils.R")
loginfo('Calculate DESeq heat map')
deseq_distmat <- heatmap_v (count, "DESeq", condition, out)
loginfo('DESeq heat map is done')
loginfo('Calculate EBSeq heat map')
ebseq_distmat <- heatmap_v (count, "EBSeq", condition, out)
loginfo('EBSeq heat map is done')
loginfo('Calculate edgeR heat map')
edger_distmat <- heatmap_v (count, "edgeR", condition, out)
loginfo('edgeR heat map is done')
loginfo('Calculate voom heat map')
voom_distmat <- heatmap_v (count, "voom", condition, out)
loginfo('voom heat map is done')
loginfo('Calculate NOISeq heat map')
noiseq_distmat <- heatmap_v (count, "NOISeq", condition, out)
loginfo('NOISeq heat map is done')
loginfo('Heat maps for tools are made')
dev.off()

loginfo('The end!')
#pdf(file.path(out, "de_plots.pdf"))
#deseq2_v(deseq2_res, out)
#ebseq_v(ebseq_res, out)
#edger_v(edger_res, out)
#voom_v(voom_res, out)
#noiseq_v(noiseq_res, out)
#dev.off()
#loginfo('Visualization results of differential expression is done')
