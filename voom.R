
#use tximeta
library("BiocManager")
library("dplyr")
library("tximeta")
library("SummarizedExperiment")

counts<- read.table(file="data/TCGA_prostate_countmatrix.txt", header = TRUE, sep = ",", dec = ".")
gene_names <- counts$X
counts <- counts[,-c(1)]
rownames(counts) <- gene_names
coldata <- read.table(file="data/annotation_TCGA_prostate.txt", sep = ",", dec = ".")
colnames(coldata) <- c("names", "condition")
coldata <- coldata[-c(1), ]
condition <- coldata$condition


library("edgeR")
library("limma")

d <- DGEList(counts)

#calculate diff expression
design <- model.matrix(~condition + 0, data = coldata)
y <- voom(d, design)
cat ("VOOMED\n")
fit <- lmFit(y, design)
head(coef(fit))


contrast <- c(unique(coldata[c("condition")])$condition[1], unique(coldata[c("condition")])$condition[2])
contrast[1] <- paste("condition", contrast[1], sep="")
contrast[2] <- paste("condition", contrast[2], sep="")
contrast <- paste(contrast[1], " - ", contrast[2])
contrasts <- makeContrasts(contrast, levels = colnames(coef(fit)))
contrasts

fit2 <- contrasts.fit(fit, contrasts)
fit2 <- eBayes(fit2, trend=TRUE)
res <- topTable(fit2, sort.by = "P", n = Inf)
head(res, 20)

resSig <- subset(res, adj.P.Val < 0.5)
saveRDS(resSig, file = "data/voom-res.rds")

library(EnhancedVolcano)
EnhancedVolcano(resSig,
		lab = rownames(resSig),
		x = 'logFC',
		y = 'P.Value',
		pCutoff = 10e-4,
		FCcutoff = 1.5,
		title = "limma-voom results",
		subtitle = "Differential expression")

