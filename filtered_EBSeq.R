ebseq_filtered <- function(filename) {
  library("EBSeq")
  
  # reading results of diff. expression analysis from RDS file
  results <- readRDS(file = filename)
  
  # filtering results by PPDE >= 0.95
  ppde <- results$PPDE
  ppde <- ppde[ppde >= 0.95]
  ppde <- as.data.frame(ppde)
  colnames(ppde) <- c("PPDE")
  filtered_genes <- gsub("\\..*","",row.names(ppde))
  return(filtered_genes)
}