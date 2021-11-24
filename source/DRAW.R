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
source("source/calculate_distmatrix_utils.R")

pdf(file.path(out, "heatmaps.pdf"))
# Visualize using heat maps
loginfo('Calculate heat maps for tools')
loginfo('Calculate DESeq heat map')
deseq_hm <- heatmap_v (count, "DESeq", condition, ', |log2FC| < 0.25, sorted by p-value', out)
loginfo('DESeq heat map is done')
loginfo('Calculate EBSeq heat map')
ebseq_hm <- heatmap_v (count, "EBSeq", condition, ', sorted by PostFC', out)
loginfo('EBSeq heat map is done')
loginfo('Calculate edgeR heat map')
edger_hm <- heatmap_v (count, "edgeR", condition, ', |log2FC| < 0.25, sorted by p-value', out)
loginfo('edgeR heat map is done')
loginfo('Calculate voom heat map')
voom_hm <- heatmap_v (count, "voom", condition, ', |log2FC| < 0.25, sorted by p-value', out)
loginfo('voom heat map is done')
loginfo('Calculate NOISeq heat map')
noiseq_hm <- heatmap_v (count, "NOISeq", condition, ', q > 0.8, sorted by probability', out)
loginfo('NOISeq heat map is done')
loginfo('Heat maps for tools are made')
dev.off()