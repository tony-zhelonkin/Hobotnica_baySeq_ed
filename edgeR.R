# Calculate diff expression analysis
edger_f <- function(counts, coldata) {
    # Import libraries
    library(BiocManager)
    library(dplyr)
    library(tximeta)
    library(SummarizedExperiment)
    library(edgeR)

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
edger_v <- function(edger_res) {
    library(BiocManager)
    library(dplyr)
    library(tximeta)
    library(SummarizedExperiment)
    library(edgeR)
    library(EnhancedVolcano)
    pdf("data/edger_plot.pdf")
    EnhancedVolcano(edger_res,
            lab = rownames(edger_res),
            x = 'logFC',
            y = 'PValue',
            pCutoff = 10e-15,
            FCcutoff = 4,
            title = "edgeR results",
            subtitle = "Differential expression")
}

# Make a signature of top-20 genes
edger_top <- function(results) {
    library("edgeR")
    library("biomaRt")

    # Filter results by logFC > 1 or logFC < -1
    filtered_results <- results[abs(results$logFC) > 1, ]

    # Extract top-20 differentially expressed genes ordered by p-value
    top <- head(filtered_results[order(filtered_results$PValue), ], 20)
    tmp <- gsub("\\..*","",row.names(top))

    # Write top-20 genes with original (ENSEMBL) encoding
    return(tmp)

}