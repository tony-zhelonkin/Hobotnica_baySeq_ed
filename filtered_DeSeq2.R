deseq2_filtered <- function(filename) {
  library("DESeq2")
  # reading results of diff. expression analysis from RDS file
  results <- readRDS(file = filename)
  
  # filtering results by log2FC >= 2 and p-value < 0.05
  filtered_results <- results[abs(results$log2FoldChange) >= 2 && results$pvalue < 0.05, ]
  
  return(filtered_results)
}

