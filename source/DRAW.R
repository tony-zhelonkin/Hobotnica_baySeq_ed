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

h_results <- read.table(file=file.path(out, "hobotnica_scores.txt"), sep = " ", dec = ".")
colnames(h_results) <- c("names", "scores")
h_results <- h_results[order(h_results$scores, decreasing = TRUE), ]

best_genes <- readLines(file.path(out,paste0(h_results$names[1], "_sig.txt")))
top2_cross <- best_genes %in% readLines(file.path(out,paste0(h_results$names[2], "_sig.txt")))
top3_cross <- best_genes %in% readLines(file.path(out,paste0(h_results$names[3], "_sig.txt")))
top4_cross <- best_genes %in% readLines(file.path(out,paste0(h_results$names[4], "_sig.txt")))
top5_cross <- best_genes %in% readLines(file.path(out,paste0(h_results$names[5], "_sig.txt")))
cross_results <- data.frame(best_genes, ifelse(top2_cross, '+', '-'), ifelse(top3_cross, '+', '-'),
    ifelse(top4_cross, '+', '-'), ifelse(top5_cross, '+', '-'))
colnames(cross_results) <- h_results$names

