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

    return(distMatrix)
}