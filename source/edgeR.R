# Calculate diff expression analysis
edger_f <- function(counts, coldata) {
    # Import libraries
    suppressMessages(library(BiocManager))
    suppressMessages(library(dplyr))
    suppressMessages(library(tximeta))
    suppressMessages(library(SummarizedExperiment))
    suppressMessages(library(edgeR))

    # Prepare data
    condition <- coldata$condition
    dgList <- DGEList(counts)

    # Filter data
    countsPerMillion <- cpm(dgList)
    countCheck <- countsPerMillion > 1
    keep <- which(rowSums(countCheck) >= 2)
    dgList <- dgList[keep,]

    # Normalize data
    dgList <- calcNormFactors(dgList, method="TMM")

    # Calculate diff expression
    design <- model.matrix(~condition, data = coldata)
    dgList <- estimateGLMCommonDisp(dgList, design=design)
    fit <- glmFit(dgList, design)
    lrt <- glmLRT(fit)
    res <- topTags(lrt, n = Inf)
    res <- res$table

    # To decrease number of genes use: res <- subset(res, PValue < 0.5)

    # Summary for results
    summary(res)

    # Return results of diff expression analysis
    return(res)
}

# Visualize function
edger_v <- function(edger_res, out) {
    suppressMessages(library(BiocManager))
    suppressMessages(library(dplyr))
    suppressMessages(library(tximeta))
    suppressMessages(library(SummarizedExperiment))
    suppressMessages(library(edgeR))
    suppressMessages(library(EnhancedVolcano))
    vis_res <-  EnhancedVolcano(edger_res,
            lab = rownames(edger_res),
            x = 'logFC',
            y = 'PValue',
            pCutoff = 10e-15,
            FCcutoff = 4,
            title = "edgeR results",
            subtitle = "|log2FC| < 0.25, sorted by p-value ")
    return(vis_res)
}

# Make a signature of top-n genes
edger_top <- function(results, n) {
    # Import libraries
    suppressMessages(library("edgeR"))
    suppressMessages(library("biomaRt"))

    # Filter results by logFC > 0.25 or logFC < -0.25
    filtered_results <- results[!is.na(results$logFC) > 0 && abs(results$logFC) > 0.25, ]

    # Extract top-n differentially expressed genes ordered by p-value
    top <- head(filtered_results[order(filtered_results$PValue), ], n)
    tmp <- gsub("\\..*","",row.names(top))

    # Write top-n genes with original (ENSEMBL) encoding
    return(tmp)

}

# Filter genes by logFC and p-value
edger_filtered <- function(results) {
    suppressMessages(library("edgeR"))

    # filtering results by log2FC >= 2 and p-value < 0.05
    filtered_results <- results[abs(results$logFC) >= 2 && results$PValue < 0.05, ]
    filtered_genes <- gsub("\\..*","",row.names(filtered_results))
    return(filtered_genes)
}