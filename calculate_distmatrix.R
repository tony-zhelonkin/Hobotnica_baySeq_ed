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
}

# Calculate dist matrixes for all tools
calculate_distmatrixes <- function(countMatrixFile) {
    # Calculate for all tools
    calculate_distmatrix (countMatrixFile, "data/DESeq_sig.txt")
    calculate_distmatrix (countMatrixFile, "data/EBSeq_sig.txt")
    calculate_distmatrix (countMatrixFile, "data/edgeR_sig.txt")
    calculate_distmatrix (countMatrixFile, "data/voom_sig.txt")
}


args = commandArgs(trailingOnly=TRUE)
countMatrixFile <- args[1]
calculate_distmatrixes(countMatrixFile)