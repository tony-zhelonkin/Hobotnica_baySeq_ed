library("DESeq2")
library("biomaRt")

#extracting results of analysis
results <- readRDS(file = "data/DESeq-res.rds")

#filtering results by logFC > 1 or logFC < -1
filtered_results <- results[abs(results$log2FoldChange) > 1, ]

#extracting top-20 differentially expressed genes ordered by p-value
top <- head(filtered_results[order(filtered_results$pvalue), ], 20)
tmp <- gsub("\\..*","",row.names(top))

#performing symbol annotation of top-20 genes 
mart <- useMart("ENSEMBL_MART_ENSEMBL")
mart <- useDataset("hsapiens_gene_ensembl", mart)
annot <- getBM(mart, attributes = c("ensembl_gene_id", "external_gene_name"), 
               filter = "ensembl_gene_id", values = tmp, uniqueRows = TRUE)
top_genes <- annot[annot$external_gene_name != "", ]$external_gene_name
top_genes

#writing results to file
writeLines(top_genes, "data/DESeq-top.txt")