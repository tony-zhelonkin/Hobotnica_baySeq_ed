
# Calculate hobotnica score for one tool
calculate_hobotnica <- function(distmat_filename, condition) {
    distmat <- read.table(file = distmat_filename, header=TRUE, sep=",")
    distmat <- as.matrix(data.frame(distmat[, -1], row.names = distmat[, 1]))
    colnames(distmat) <- rownames(distmat)
    source("hobotnica/R/Hobotnica.R")
    return(Hobotnica(distmat, condition))
}

# Calculate hobotnica scores for all tools and write it in the files
use_hobotnica <- function(condition, out) {
    results <- c()
    results <- c(results, paste0("DESeq2 ", calculate_hobotnica(file.path(out, "DESeq_sig.txt.distmatrix"), condition)))
    results <- c(results, paste0("EBSeq ", calculate_hobotnica(file.path(out, "EBSeq_sig.txt.distmatrix"), condition)))
    results <- c(results, paste0("edgeR ", calculate_hobotnica(file.path(out, "edgeR_sig.txt.distmatrix"), condition)))
    results <- c(results, paste0("voom ", calculate_hobotnica(file.path(out, "voom_sig.txt.distmatrix"), condition)))
    results <- c(results, paste0("NOISeq ", calculate_hobotnica(file.path(out, "NOISeq_sig.txt.distmatrix"), condition)))
    writeLines(results, file.path(out, "hobotnica_scores.txt"))
    return(results)
}
