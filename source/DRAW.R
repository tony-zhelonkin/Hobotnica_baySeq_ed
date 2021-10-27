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

n <- 30
deseq_sig <- readLines(file.path(out,paste0("DESeq_sig.txt")))
ebseq_sig <- readLines(file.path(out,paste0("EBSeq_sig.txt")))
edger_sig <- readLines(file.path(out,paste0("edgeR_sig.txt")))
noiseq_sig <- readLines(file.path(out,paste0("NOISeq_sig.txt")))
voom_sig <- readLines(file.path(out,paste0("voom_sig.txt")))
sigs <- unique(c(deseq_sig, ebseq_sig, edger_sig, noiseq_sig, voom_sig))
deseq_cross <- ifelse(sigs %in% deseq_sig, '+', '-')
ebseq_cross <- ifelse(sigs %in% ebseq_sig, '+', '-')
edger_cross <- ifelse(sigs %in% edger_sig, '+', '-')
noiseq_cross <- ifelse(sigs %in% noiseq_sig, '+', '-')
voom_cross <- ifelse(sigs %in% voom_sig, '+', '-')
count_cross <- ifelse(deseq_cross == '+', 1, 0) + ifelse(ebseq_cross == '+', 1, 0) +
                ifelse(edger_cross == '+', 1, 0) + ifelse(noiseq_cross == '+', 1, 0) +
                ifelse(voom_cross == '+', 1, 0)
cross_results <- data.frame(deseq_cross, ebseq_cross, edger_cross,
    noiseq_cross, voom_cross, count_cross)
colnames(cross_results) <- c("DESeq", "EBSeq", "edgeR", "NOISeq", "voom", "count")
rownames(cross_results) <- sigs
cross_results <- cross_results[order(cross_results$count, decreasing = TRUE), ]
cross_results
write.table(cross_results, file=file.path(out,"crossing.csv"), row.names = FALSE)