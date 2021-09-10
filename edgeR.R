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

dgList <- makeDGEList(gse)

#filtering
countsPerMillion <- cpm(dgList)
countCheck <- countsPerMillion > 1
keep <- which(rowSums(countCheck) >= 2)
dgList <- dgList[keep,]

#normalization
dgList <- calcNormFactors(dgList, method="TMM")



#calculate diff expression
design <- model.matrix(~conditions_dexamethasone, data = coldata)

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

library(EnhancedVolcano)
EnhancedVolcano(resSig,
		lab = rownames(resSig),
		x = 'logFC',
		y = 'PValue',
		pCutoff = 10e-15,
		FCcutoff = 4,
		title = "edgeR results",
		subtitle = "Differential expression")

