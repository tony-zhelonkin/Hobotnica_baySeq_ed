
#use tximeta
library("BiocManager")
library("dplyr")
library("tximeta")
library("SummarizedExperiment")

gse <- readRDS(file = "data/matrix.rds")
coldata <- readRDS(file = "data/coldata.rds")
library("edgeR")
library("limma")

condition <- coldata$condition

d <- makeDGEList(gse)

#calculate diff expression
design <- model.matrix(~condition + 0, data = coldata)
y <- voom(d, design)
cat ("VOOMED\n")
fit <- lmFit(y, design)
head(coef(fit))

contrasts <- makeContrasts(conditionUntreated - conditionTreated, levels = colnames(coef(fit)))
contrasts
fit2 <- contrasts.fit(fit, contrasts)
fit2 <- eBayes(fit2, trend=TRUE)
res <- topTable(fit2, sort.by = "P", n = Inf)
head(res, 20)

resSig <- subset(res, adj.P.Val < 0.5)

library(EnhancedVolcano)
EnhancedVolcano(resSig,
		lab = rownames(resSig),
		x = 'logFC',
		y = 'P.Value',
		pCutoff = 10e-4,
		FCcutoff = 1.5,
		title = "limma-voom results",
		subtitle = "Differential expression")

