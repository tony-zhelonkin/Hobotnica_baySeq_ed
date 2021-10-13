# Calculate diff expression analysis
noiseq_f <- function(counts, coldata) {
    # Import libraries
    library(BiocManager)
    library(dplyr)
    library(tximeta)
    library(SummarizedExperiment)
    library(NOISeq)

    # Prepare data
    factors <- as.data.frame(coldata$condition)
    colnames(factors) <- c("condition")
    noiseq_data <- readData(data = counts, factors = factors)
    head(noiseq_data)

    # Calculate diff expression
    noiseq_res <- noiseqbio(noiseq_data, factor = "condition")

    # Summary for results
    #summary(res)

    # Return results of diff expression analysis
    return(noiseq_res)
}

# Make a signature of top-n genes
noiseq_top <- function(results, n) {
    library("NOISeq")
    library("biomaRt")

    results <- degenes(results)
    # Filter results by logFC > 1 or logFC < -1
    filtered_results <- results[!is.na(results$log2FC) > 0 && abs(results$log2FC) > 1, ]

    # Extract top-n differentially expressed genes ordered by p-value
    top <- head(filtered_results[order(filtered_results$prob), ], n)
    tmp <- gsub("\\..*","",row.names(top))

    # Write top-n genes with original (ENSEMBL) encoding
    return(tmp)

}
