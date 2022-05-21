# Calculate diff expression analysis
bayseq_f <- function(counts, coldata) {

    # Import libraries
    suppressMessages(library(BiocManager))
    suppressMessages(library(dplyr))
    suppressMessages(library(tximeta))
    suppressMessages(library(SummarizedExperiment))
    suppressMessages(library(baySeq))
    suppressMessages(library(snow))
    
    
    # Prepare data
    gene_id <- rownames(counts)
    regroup <- c(sort(unique(names(counts))[1]), sort(unique(names(counts))[-1]))
    
   
    counts <- counts[, regroup] # reorder groups
    
    replicates <- coldata$condition
    
    # Group 1 technical replicates
    n_ref <- length(replicates[replicates == unique(replicates)[1]]) # number of control replicates
    ref_list <- rep(x = unique(replicates)[1], times = n_ref) # vector of control replicates
    
    # Group 2 technical replicates
    n_exp <- length(replicates[replicates == unique(replicates)[2]]) # number of experimental replicates
    exp_list <- rep(x = unique(replicates)[2], times = n_exp) # number of experimental replicates
    
    # Create grouping variable
    groups <- list(NDE = c(ref_list, ref_list), DE = c(ref_list, exp_list)) # NDE = prior distribution, DE - posterior 
    
    cname <- rownames(counts) # gene names
    
    reads_matrix <- as.matrix(counts) # final matrix of reads for DE anaysis
    
    # Calculate diff expression
    CD <- new("countData", data = reads_matrix, replicates = replicates, groups = groups)
    
    libsizes(CD) <- getLibsizes(CD) # estimate library sizes
    
    CD@annotation <- as.data.frame(cname) # annotate the object with gene names
    
    cl <- makeCluster(8, "SOCK") # create a cluster of parallel computation with snow library
    
    CD <- getPriors.NB(CD, samplesize = 10000, estimation = "QL", cl = cl) # create apriori distribution
    CD <- getLikelihoods(CD, pET = 'BIC', cl = cl) # aposteriori distribution
    
    bay_res <- topCounts(CD, group="DE", number = length(cname)) # all baySeq statistics
    
    # Estimating log2fold change
    log_data <- log2(CD@data)
    control = apply(log_data[,1:3], 1, mean) #calculate the mean of each gene per control group
    test = apply(log_data[,4:6], 1, mean) #calcuate the mean of each gene per test group
    
    log2FoldChange <- control - test # log2(control / test) == log2FoldChange
    
    # Create results data frame
    fold_frame <- cbind(cname, log2FoldChange)
    fold_frame <- as.data.frame(fold_frame)
    
    # Preparing final results
    bay_res <- merge.data.frame(x = bay_res, y = fold_frame, by.y = "cname")
    bay_res$log2FoldChange <- as.numeric(bay_res$log2FoldChange)
    bay_res <- bay_res[order(bay_res$FDR.DE), ] # sort by FDR
    rownames(bay_res) <- bay_res$cname
    
    # substitute Inf with NA to make a volcano plot
    bay_res <- do.call(data_frame, lapply(bay_res, function(x) replace(x, is.infinite(x),NA)))
    
    
    # Summary for results
    summary(bay_res)
    
    # Return results of diff expression analysis
    return(bay_res)
    
  }
  
  # Visualize function
  bayseq_v <- function(bayseq_res, out) {
    suppressMessages(library(BiocManager))
    suppressMessages(library(dplyr))
    suppressMessages(library(tximeta))
    suppressMessages(library(SummarizedExperiment))
    suppressMessages(library(baySeq))
    suppressMessages(library(EnhancedVolcano))
    
    vis_res <- EnhancedVolcano(bayseq_res,
                               lab = bayseq_res$cname,
                               x = 'log2FoldChange',
                               y = 'FDR.DE',
                               pCutoff = 10e-4,
                               FCcutoff = 3,
                               title = "baySeq results",
                               subtitle = "|log2FC| < 0.25, sorted by p-value")
    return(vis_res)
  }
  
  
  # Make a signature of top-n genes
  bayseq_top <- function(results, n) {
    # Import libraries
    suppressMessages(library("baySeq"))
    
    # Filter results by logFC > 0.25 or logFC < -0.25
    filtered_results <- results[!is.na(results$log2FoldChange) > 0 && abs(results$log2FoldChange) > 0.25, ]
    
    # Extract top-n differentially expressed genes ordered by p-value
    top <- head(filtered_results, n)
    rownames(top) <- top$cname
    tmp <- gsub("\\..*","",row.names(top))
    
    # Write top-n genes with original (ENSEMBL) encoding
    return(tmp)
    
  }
  
  # Filter genes by logFC and p-value
  bayseq_filtered <- function(results) {
    suppressMessages(library("baySeq"))
    
    # filtering results by log2FC >= 2 and p-value < 0.05
    filtered_results <- results[abs(results$log2FoldChange) >= 2 && results$FDR.DE < 0.05, ]
    rownames(filtered_results) <- filtered_results$cname
    filtered_genes <- gsub("\\..*","",row.names(filtered_results))
    return(filtered_genes)
  }