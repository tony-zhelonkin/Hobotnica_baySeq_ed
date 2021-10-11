voom_filtered <- function(filename) {
  library("edgeR")
  # reading results of diff. expression analysis from RDS file
  results <- readRDS(file = filename)
  
  # filtering results by log2FC >= 2 and p-value < 0.05
  filtered_results <- results[abs(results$logFC) >= 2 && results$P.Value < 0.05, ]
  filtered_genes <- gsub("\\..*","",row.names(filtered_results))
  return(filtered_genes)
}