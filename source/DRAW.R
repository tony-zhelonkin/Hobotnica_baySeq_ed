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

# Visualize results of differential expression
loginfo('Visualize results of differential expression')
pdf(file.path(out, "de_plots.pdf"))
library(gridExtra)
library(grid)
library(ggplot2)
library(lattice)
deseq2_v(deseq2_res, out)
#a <- ebseq_v(ebseq_res, out)
b <- edger_v(edger_res, out)
c <- voom_v(voom_res, out)
#d <- noiseq_v(noiseq_res, out)
grid.arrange(b, c, ncol=2,     top = textGrob('EnhancedVolcano',
      just = c('center'),
      gp = gpar(fontsize = 32)))
dev.off()
loginfo('Visualization results of differential expression is done')
