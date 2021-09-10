#use tximeta
library("BiocManager")
library("dplyr")
library("tximeta")
library("SummarizedExperiment")

gse <- readRDS(file = "data/matrix.rds")
library("DESeq2")
dds <- DESeqDataSet(gse, design = ~condition)

#calculate diff expression
dds <- DESeq(dds)
res <- results(dds)
res <- lfcShrink(dds, contrast=c("condition", "Untreated", "Treated"),
                 res=res, type='normal')
res <- results(dds, contrast=c("condition", "Untreated", "Treated"))
resSig <- subset(res, padj < 0.5)
summary(resSig)

saveRDS(resSig, file = "data/DESeq-res.rds")

library(EnhancedVolcano)
EnhancedVolcano(resSig,
		lab = rownames(resSig),
		x = 'log2FoldChange',
		y = 'pvalue',
		pCutoff = 10e-4,
		FCcutoff = 3,
		title = "DESeq2 results",
		subtitle = "Differential expression")
