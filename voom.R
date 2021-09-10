#Find data
sample_name <- "SRR1039"
location <- paste("quants/", sample_name, sep = "", collapse = "")

#Create data array
names <- c()
files <- c()
conditions <- c()
cell_line <- c()
for (i in "508":"523") {
    names <- c(names, paste(sample_name, i, sep = "", collapse = ""))
    files <- c(files, paste(location, i, "_quant/quant.sf", sep = "", collapse = ""))
}
conditions_dexamethasone <- c(conditions, "Untreated", "Treated", "Untreated", "Treated", "Untreated", "Treated", "Untreated", "Treated", "Untreated", "Treated", "Untreated", "Treated", "Untreated", "Treated", "Untreated", "Treated")
conditions_albuterol <- c(conditions, "Untreated", "Untreated", "Treated", "Treated", "Untreated", "Untreated", "Treated", "Treated", "Untreated", "Untreated", "Treated", "Treated", "Untreated", "Untreated", "Treated", "Treated")
cell_line <- c(cell_line, "N61311", "N61311", "N61311", "N61311", "N052611", "N052611", "N052611", "N052611", "N080611", "N080611", "N080611", "N080611", "N061011", "N061011", "N061011", "N061011")
coldata <- data.frame(names = names, files = files, conditions_dexamethasone = conditions_dexamethasone, conditions_albuterol = conditions_albuterol, cell_line = cell_line, stringsAsFactors=FALSE)
coldata


#use tximeta
library("BiocManager")
library("dplyr")
library("tximeta")
library("SummarizedExperiment")

gse <- readRDS(file = "matrix.rds")

library("edgeR")
library("limma")

d <- makeDGEList(gse)

#calculate diff expression
design <- model.matrix(~conditions_dexamethasone + 0, data = coldata)
y <- voom(d, design)
cat ("VOOMED\n")
fit <- lmFit(y, design)
head(coef(fit))

contrasts <- makeContrasts(conditions_dexamethasoneUntreated - conditions_dexamethasoneTreated, levels = colnames(coef(fit)))
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

