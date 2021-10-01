#use tximeta
library("BiocManager")
library("dplyr")
library("tximeta")
library("SummarizedExperiment")
library("DESeq2")

counts<- read.table(file="data/TCGA_prostate_countmatrix.txt", header = TRUE, sep = ",", dec = ".")
gene_names <- counts$X
counts <- counts[,-c(1)]
rownames(counts) <- gene_names
coldata <- read.table(file="data/annotation_TCGA_prostate.txt", sep = ",", dec = ".")
colnames(coldata) <- c("names", "condition")
rownames(coldata) <- coldata$names
coldata <- coldata[-c(1), ]
condition <- coldata$condition

contrast <- c("condition", unique(coldata[c("condition")])$condition[1], unique(coldata[c("condition")])$condition[2])
contrast

dds <- DESeqDataSetFromMatrix(countData=as.matrix(counts), colData=coldata, design = ~condition)

#calculate diff expression
dds <- DESeq(dds)
res <- results(dds)
res <- lfcShrink(dds, contrast=contrast,
                 res=res, type='normal')
res <- results(dds, contrast=contrast)
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
