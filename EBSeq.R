# Calculate diff expression analysis
ebseq_f <- function(counts, coldata) {
    # Import libraries
    library(BiocManager)
    library(dplyr)
    library(tximeta)
    library(SummarizedExperiment)
    library(EBSeq)

    # Prepare data
    condition <- coldata$condition

    # Calculate diff expression
    Sizes = MedianNorm(as.matrix(counts))
    EBOut = EBTest(Data = as.matrix(counts),
                    Conditions = as.factor(condition),
                    sizeFactors = Sizes, maxround = 5)
    Out = GetDEResults(EBOut)
    res <- Out$PPMat

    # Summary for results
    summary(res)

    # Return results of diff expression analysis
    return(EBOut)
}

# Visualize function
ebseq_v <- function(ebseq_res) {
    library(BiocManager)
    library(dplyr)
    library(tximeta)
    library(SummarizedExperiment)
    library(EBSeq)
    library(EnhancedVolcano)
    res <- GetDEResults(ebseq_res)$PPMat
    GeneFC <- PostFC(ebseq_res, SmallNum = 0.01)    # Create dataframe with groups
    x = log2(GeneFC$RealFC)
    y = ebseq_res$PPDE

    xcond = 2   # Change group colors
    ycond = 0.7
    group <- ifelse((x < -xcond | x > xcond) & y > ycond , "DE", ifelse(y > ycond , "ADE2", ifelse((x < -xcond | x > xcond), "ADE1", "NDE")))
    df <- data.frame(x = x, y = y, group = factor(group))
    colors <- c("red", "blue", "green", "gray")

    pdf("data/ebseq_plot.pdf")   # Make plot
    plot(df$x, df$y, col = colors[df$group], pch = 16,
        xlim =c(-5,5),
        ylim=c(0,1),
        main="EBSeq results",
        xlab="EBSeq Log2 Fold Change", ylab="EBSeq PPDE")
    DEgplots <- data.frame(x = ifelse((x < -xcond | x > xcond) & y > ycond, x, 0), y = ifelse((x < -xcond | x > xcond) & y > ycond , y, 0),
        name = ifelse((x < -xcond | x > xcond) & y > ycond , rownames(res), ""))
    text(DEgplots$x, DEgplots$y,
         labels = DEgplots$name,
         cex = 0.6, pos = 4, col = "black")
    dev.off()
}