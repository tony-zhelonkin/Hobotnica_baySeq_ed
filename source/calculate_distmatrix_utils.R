# Calculate certain dist matrix
calculate_distmatrix <- function(countMatrixFile, signatureFile) {
    # Download library
    library("amap")

    # Download data from files
    signature <- readLines(signatureFile)
    cm <- read.table(countMatrixFile, header=T, sep=",")
    cm <- data.frame(cm[, -1], row.names=cm[, 1])

    # Calculate dist matrix
    cm_subset <- cm[signature, ]
    cat("SUBSET SHAPE: ")
    cat(dim(cm_subset))
    cat("\n")
    distMatrix <- Dist(t(cm_subset), method="kendall", nbproc=6)

    # Save dist matrix
    write.csv(as.data.frame(as.matrix(distMatrix)), file=paste0(signatureFile, ".distmatrix"))

    return(0)
}

# Calculate dist matrixes for all tools
calculate_distmatrixes <- function(countMatrixFile, out) {
    # Calculate for all tools
    calculate_distmatrix (countMatrixFile, file.path(out, "DESeq_sig.txt"))
    calculate_distmatrix (countMatrixFile, file.path(out, "EBSeq_sig.txt"))
    calculate_distmatrix (countMatrixFile, file.path(out, "edgeR_sig.txt"))
    calculate_distmatrix (countMatrixFile, file.path(out, "voom_sig.txt"))
    calculate_distmatrix (countMatrixFile, file.path(out, "NOISeq_sig.txt"))
    return(0)
}
