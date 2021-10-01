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

dgList <- DGEList(counts)

#filtering
countsPerMillion <- cpm(dgList)
countCheck <- countsPerMillion > 1
keep <- which(rowSums(countCheck) >= 2)
dgList <- dgList[keep,]

#normalization
dgList <- calcNormFactors(dgList, method="TMM")



#calculate diff expression
design <- model.matrix(~condition, data = coldata)

dgList <- estimateGLMCommonDisp(dgList, design=design)

fit <- glmFit(dgList, design)
lrt <- glmLRT(fit)
cat ("EDGERed\n")
head(coef(lrt))

res <- topTags(lrt, n = 1000000)
head(res, 15)
resSig <- res$table
summary(resSig)
#resSig <- subset(resSig, PValue < 0.5)
saveRDS(resSig, file = "data/edgeR-res.rds")


library(EnhancedVolcano)
EnhancedVolcano(resSig,
		lab = rownames(resSig),
		x = 'logFC',
		y = 'PValue',
		pCutoff = 10e-15,
		FCcutoff = 4,
		title = "edgeR results",
		subtitle = "Differential expression")

