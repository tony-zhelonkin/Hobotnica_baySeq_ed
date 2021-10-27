# Calculate diff expression analysis
deseq2_f <- function(counts, coldata) {
    # Import libraries
    suppressMessages(library(BiocManager))
    suppressMessages(library(dplyr))
    suppressMessages(library(tximeta))
    suppressMessages(library(SummarizedExperiment))
    suppressMessages(library(DESeq2))

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
deseq2_v <- function(deseq2_res, out) {
    suppressMessages(library(BiocManager))
    suppressMessages(library(dplyr))
    suppressMessages(library(tximeta))
    suppressMessages(library(SummarizedExperiment))
    suppressMessages(library(DESeq2))
    suppressMessages(library(EnhancedVolcano))

    vis_res <- EnhancedVolcano(deseq2_res,
            lab = rownames(deseq2_res),
            x = 'log2FoldChange',
            y = 'pvalue',
            pCutoff = 10e-4,
            FCcutoff = 3,
            title = "DESeq2 results",
            subtitle = "Differential expression")
    return(vis_res)
}

# Make a signature of top-n genes
deseq2_top <- function(results, n) {
    # Import libraries
    suppressMessages(library("DESeq2"))
    suppressMessages(library("biomaRt"))

    # Filter results by logFC > 1 or logFC < -1
    filtered_results <- results[!is.na(results$log2FoldChange) > 0 && abs(results$log2FoldChange) > 0.25, ]

    # Extract top-n differentially expressed genes ordered by p-value
    top <- head(filtered_results[order(filtered_results$pvalue), ], n)
    tmp <- gsub("\\..*","",row.names(top))

    # Write top-n genes with original (ENSEMBL) encoding
    return(tmp)

}

# Filter genes by logFC and p-value
deseq2_filtered <- function(results) {
  suppressMessages(library("DESeq2"))

  # filtering results by log2FC >= 2 and p-value < 0.05
  filtered_results <- results[abs(results$log2FoldChange) >= 2 && results$pvalue < 0.05, ]
  filtered_genes <- gsub("\\..*","",row.names(filtered_results))
  return(filtered_genes)
}