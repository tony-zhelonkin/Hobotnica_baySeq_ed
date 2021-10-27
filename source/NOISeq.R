# Calculate diff expression analysis
noiseq_f <- function(counts, coldata) {
    # Import libraries
    suppressMessages(library(BiocManager))
    suppressMessages(library(dplyr))
    suppressMessages(library(tximeta))
    suppressMessages(library(SummarizedExperiment))
    suppressMessages(library(NOISeq))

    # Prepare data
    factors <- as.data.frame(coldata$condition)
    colnames(factors) <- c("condition")
    noiseq_data <- readData(data = counts, factors = factors)
    head(noiseq_data)

    # Calculate diff expression
    noiseq_res <- noiseq(noiseq_data, factor = "condition")

    # Summary for results
    #summary(res)

    # Return results of diff expression analysis
    return(noiseq_res)
}

# Visualize function
noiseq_v <- function(noiseq_res, out) {
    suppressMessages(library(BiocManager))
    suppressMessages(library(dplyr))
    suppressMessages(library(tximeta))
    suppressMessages(library(SummarizedExperiment))
    suppressMessages(library(NOISeq))
    suppressMessages(library(EnhancedVolcano))
    
    # Built in graphic function
    # DE.plot(noiseq_res, q = 0.8, graphic = "MD")
  
    noiseq_res <- noiseq_res@results[[1]]
    x = noiseq_res$M
    y = noiseq_res$D
    xycond = noiseq_res$prob
    group <- ifelse(xycond > 0.8, "DE", "NDE")
    df <- data.frame(x = x, y = y, group = factor(group))
    colors <- c("red", "gray")

    pdf(file.path(out, "noiseq_plot.pdf"))
    plot(df$x, df$y, col = colors[df$group], pch = 16,
        ylim = c(0,10000),
        main="NOISeq results",
        xlab="M", ylab="D")
    DEgplots <- data.frame(x = ifelse(xycond > 0.8, x, 0), y = ifelse(xycond > 0.8, y, 0),
        name = ifelse(xycond > 0.8, rownames(noiseq_res), ""))
    text(DEgplots$x, DEgplots$y,
         labels = DEgplots$name,
         cex = 1.0, pos = 4, col = "black")
    dev.off()
}


# Make a signature of top-n genes
noiseq_top <- function(results, n) {
    # Import libraries
    suppressMessages(library("NOISeq"))
    suppressMessages(library("biomaRt"))

    # Filter results by q threshold
    filtered_results <- degenes(results, q = 0.8, M = NULL)

    # Extract top-n differentially expressed genes ordered by prob
    top <- head(filtered_results[order(filtered_results$prob, decreasing = TRUE), ], n)
    tmp <- gsub("\\..*","",row.names(top))

    # Write top-n genes with original (ENSEMBL) encoding
    return(tmp)

}

# Filter genes by logFC and probability
noiseq_filtered <- function(results) {
  suppressMessages(library("NOISeq"))

  # Filter results by q threshold
  filtered_results <- degenes(results, q = 0.8, M = NULL)
  filtered_genes <- gsub("\\..*","",row.names(filtered_results))
  return(filtered_genes)
}