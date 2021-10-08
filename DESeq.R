# Calculate diff expression analysis
deseq2_f <- function(counts, coldata) {
    # Import libraries
    library(BiocManager)
    library(dplyr)
    library(tximeta)
    library(SummarizedExperiment)
    library(DESeq2)

    # Prepare data
    condition <- coldata$condition
    contrast <- c("condition", unique(coldata[c("condition")])$condition[1], unique(coldata[c("condition")])$condition[2])
    dds <- DESeqDataSetFromMatrix(countData=as.matrix(counts), colData=coldata, design = ~condition)

    # Calculate diff expression
    dds <- DESeq(dds)
    res <- results(dds)
    res <- lfcShrink(dds, contrast=contrast,
                     res=res, type='normal')
    res <- results(dds, contrast=contrast)

    # To decrease number of genes use: res <- subset(res, padj < 0.5)

    # Summary for results
    summary(res)

    # Return results of diff expression analysis
    return(res)
}

# Visualize function
deseq2_v <- function(deseq2_res) {
    library(BiocManager)
    library(dplyr)
    library(tximeta)
    library(SummarizedExperiment)
    library(DESeq2)
    library(EnhancedVolcano)
    pdf("data/deseq2_plot.pdf")
    EnhancedVolcano(deseq2_res,
            lab = rownames(deseq2_res),
            x = 'log2FoldChange',
            y = 'pvalue',
            pCutoff = 10e-4,
            FCcutoff = 3,
            title = "DESeq2 results",
            subtitle = "Differential expression")
}

# Make a signature of top-20 genes
deseq2_top <- function(results) {
    library("DESeq2")
    library("biomaRt")

    # Filter results by logFC > 1 or logFC < -1
    filtered_results <- results[!is.na(results$log2FoldChange) > 0 && abs(results$log2FoldChange) > 0.25, ]

    # Extract top-20 differentially expressed genes ordered by p-value
    top <- head(filtered_results[order(filtered_results$pvalue), ], 20)
    tmp <- gsub("\\..*","",row.names(top))

    # Write top-20 genes with original (ENSEMBL) encoding
    return(tmp)

}