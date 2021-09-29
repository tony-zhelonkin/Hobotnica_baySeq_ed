library("EBSeq")
library("biomaRt")

#extracting results of analysis
results <- readRDS(file = "data/EBSeq-res.rds")
ppmat = results$PPMat
GeneFC <- PostFC(results, SmallNum = 0.01)
postfc <- as.data.frame(GeneFC$PostFC)
colnames(postfc) <- c("PostFC")

#extracting top-20 differentially expressed genes ordered by PostFC
top <- head(postfc[order(postfc$PostFC, decreasing = TRUE), , drop = FALSE], 20)
tmp <- gsub("\\..*","",row.names(top))

#performing symbol annotation of top-20 genes 
mart <- useMart("ENSEMBL_MART_ENSEMBL")
mart <- useDataset("hsapiens_gene_ensembl", mart)
annot <- getBM(mart, attributes = c("ensembl_gene_id", "external_gene_name"), 
               filter = "ensembl_gene_id", values = tmp, uniqueRows = TRUE)
top_genes <- annot[annot$external_gene_name != "", ]$external_gene_name
top_genes

#writing results to file
writeLines(top_genes, "data/EBSeq-top.txt")
