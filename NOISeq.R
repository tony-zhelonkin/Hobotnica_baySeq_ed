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
    noiseq_res <- noiseq(noiseq_data, factor = "condition")

    # Summary for results
    #summary(res)

    # Return results of diff expression analysis
    return(noiseq_res)
}

# Visualize function
noiseq_v <- function(noiseq_res) {
    library(BiocManager)
    library(dplyr)
    library(tximeta)
    library(SummarizedExperiment)
    library(NOISeq)
    library(EnhancedVolcano)

    noiseq_res <- noiseq_res@results[[1]]
    x = noiseq_res$M
    y = noiseq_res$prob

    xcond = 6  # Change group colors
    ycond = 0.90
    group <- ifelse((x < -xcond | x > xcond) & y > ycond , "DE", ifelse(y > ycond , "ADE2", ifelse((x < -xcond | x > xcond), "ADE1", "NDE")))
    df <- data.frame(x = x, y = y, group = factor(group))
    colors <- c("red", "blue", "green", "gray")

    pdf("data/noiseq_plot.pdf")   # Make plot
    plot(df$x, df$y, col = colors[df$group], pch = 16,
        xlim = c(-10,10),
        ylim = c(0,1),
        main="NOISeq results",
        xlab="NOISeq M", ylab="NOISeq probability")
    DEgplots <- data.frame(x = ifelse((x < -xcond | x > xcond) & y > ycond, x, 0), y = ifelse((x < -xcond | x > xcond) & y > ycond , y, 0),
        name = ifelse((x < -xcond | x > xcond) & y > ycond , rownames(noiseq_res), ""))
    text(DEgplots$x, DEgplots$y,
         labels = DEgplots$name,
         cex = 1.0, pos = 4, col = "black")
    dev.off()
}


# Make a signature of top-n genes
noiseq_top <- function(results, n) {
    library("NOISeq")
    library("biomaRt")

    results <- results@results[[1]]
    filtered_results <- results
    # Filter results by logFC > 1 or logFC < -1

    # Extract top-n differentially expressed genes ordered by p-value
    top <- head(filtered_results[order(filtered_results$prob, decreasing = TRUE), ], n)
    tmp <- gsub("\\..*","",row.names(top))

    # Write top-n genes with original (ENSEMBL) encoding
    return(tmp)

}

# Filter genes by logFC and probability
noiseq_filtered <- function(results) {
  library("NOISeq")

  # filtering results by probability > 0.9
  filtered_results <- degenes(results)
  filtered_genes <- gsub("\\..*","",row.names(filtered_results))
  return(filtered_genes)
}