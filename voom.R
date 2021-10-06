# Calculate diff expression analysis
voom_f <- function(counts, coldata) {
    # Import libraries
    library(BiocManager)
    library(dplyr)
    library(tximeta)
    library(SummarizedExperiment)
    library(edgeR)
    library(limma)

    # Prepare data
    condition <- coldata$condition
    d <- DGEList(counts)

    #calculate diff expression
    design <- model.matrix(~condition + 0, data = coldata)
    y <- voom(d, design)
    fit <- lmFit(y, design)
    head(coef(fit))


    contrast <- c(unique(coldata[c("condition")])$condition[1], unique(coldata[c("condition")])$condition[2])
    contrast[1] <- paste("condition", contrast[1], sep="")
    contrast[2] <- paste("condition", contrast[2], sep="")
    contrast <- paste(contrast[1], " - ", contrast[2])
    contrasts <- makeContrasts(contrasts = contrast, levels = colnames(coef(fit)))
    contrasts

    fit2 <- contrasts.fit(fit, contrasts)
    fit2 <- eBayes(fit2, trend=TRUE)
    res <- topTable(fit2, sort.by = "P", n = Inf)
    head(res, 20)

    # To decrease number of genes use: res <- subset(res, adj.P.Val < 0.5)


    # Return results of diff expression analysis
    return(res)
}

# Visualize function
voom_v <- function(voom_res) {
    library(BiocManager)
    library(dplyr)
    library(tximeta)
    library(SummarizedExperiment)
    library(edgeR)
    library(limma)
    library(EnhancedVolcano)
    pdf("data/voom_plot.pdf")
    EnhancedVolcano(voom_res,
            lab = rownames(voom_res),
            x = 'logFC',
            y = 'P.Value',
            pCutoff = 10e-4,
            FCcutoff = 1.5,
            title = "limma-voom results",
            subtitle = "Differential expression")
}