#use tximeta
library("BiocManager")
library("dplyr")
library("tximeta")
library("SummarizedExperiment")

gse <- readRDS(file = "matrix.rds")
library("DESeq2")
dds <- DESeqDataSet(gse, design = ~conditions_dexamethasone)

#calculate diff expression
dds <- DESeq(dds)
res <- results(dds)
res <- lfcShrink(dds, contrast=c("conditions_dexamethasone", "Untreated", "Treated"),
                 res=res, type='normal')
res <- results(dds, contrast=c("conditions_dexamethasone", "Untreated", "Treated"))
resSig <- subset(res, padj < 0.5)
summary(resSig)

library(EnhancedVolcano)
EnhancedVolcano(resSig,
		lab = rownames(resSig),
		x = 'log2FoldChange',
		y = 'pvalue',
		pCutoff = 10e-4,
		FCcutoff = 3,
		title = "DESeq2 results",
		subtitle = "Differential expression")
