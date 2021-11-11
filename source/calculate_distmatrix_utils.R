# Calculate certain dist matrix
calculate_distmatrix <- function(countMatrixFile, signatureFile) {
    # Download library
    suppressMessages(library("amap"))

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

heatmap_v <- function(countMatrixFile, tool_name, condition, out) {
    signatureFile <- file.path(out, paste0(tool_name, "_sig.txt"))
    signature <- readLines(signatureFile)
    cm <- read.table(countMatrixFile, header=T, sep=",")
    cm <- read.table(count, header=T, sep=",")
    cm <- data.frame(cm[, -1], row.names=cm[, 1])
    cm_subset <- cm[signature, ]

    my_sample_col <- data.frame(condition)
    row.names(my_sample_col) <- colnames(cm_subset)
    cm_subset <- cm_subset[, row.names(my_sample_col)]
    cal_z_score <- function(x){
      (x - min(x)) / (max(x) - min(x))

    }

    data_subset_norm <- t(apply(cm_subset, 1, cal_z_score))

    library('pheatmap')
    pheatmap(data_subset_norm, annotation_col = my_sample_col, annotation_names_col = FALSE, annotation_legend = TRUE,
            annotation_names_row = FALSE, drop_levels = FALSE, show_rownames = T, show_colnames = F,
            cluster_cols = F)
}