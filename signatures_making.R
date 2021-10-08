library(logging)

# Logger config
basicConfig()

# Make a signature of top-20 genes
top_signature <- function(results_deseq2, results_ebseq, results_edger, results_voom) {
    # Import top-20 creating functions
    import::here(deseq2_top, .from = DESeq.R)
    import::here(ebseq_top, .from = EBSeq.R)
    import::here(edger_top, .from = edgeR.R)
    import::here(voom_top, .from = voom.R)

    #Download results of different tools
    writeLines(deseq2_top(results_deseq2), "data/DESeq_sig.txt")
    writeLines(ebseq_top(results_ebseq), "data/EBSeq_sig.txt")
    writeLines(edger_top(results_edger), "data/edgeR_sig.txt")
    writeLines(voom_top(results_voom), "data/voom_sig.txt")
}

# this function saves to file a Venn diagram based on top-20 differentially expressed
# genes extracted from edgeR, DeSeq2, voom+limma & EBSeq instruments
draw_venn_diag <- function(top_edgeR, top_DeSeq2, top_VoomLimma, top_EBSeq) {
    library("VennDiagram")
    library("RColorBrewer")
    myCol <- brewer.pal(4, "Pastel2")
    grid.newpage()
    venn_obj <- venn.diagram(
    x = list(top_edgeR, top_DeSeq2, top_VoomLimma, top_EBSeq),
    category.names = c("edgeR", "DeSeq2", "Voom+Limma", "EBSeq"),
    height = 480,
    width = 480,
    lwd = 1,
    fill = myCol,
    cex = 1.1,
    fontface = "bold",
    fontfamily = "sans",
    cat.cex = 0.8,
    cat.fontface = "bold",
    cat.default.pos = "outer",
    cat.fontfamily = "sans",
    #cat.pos = c(0, -30, -130, 150),
    filename = NULL
    )
    grid.draw(venn_obj)
}

loginfo('Make 5 signatures')
# Load differential expression analysis results
results_deseq2 <- readRDS(file = "data/DESeq_res.rds")
results_ebseq <- readRDS(file = "data/EBSeq_res.rds")
results_edger <- readRDS(file = "data/edgeR_res.rds")
results_voom <- readRDS(file = "data/voom_res.rds")

top_signature(results_deseq2, results_ebseq, results_edger, results_voom)

loginfo('Visualize results')
# read signatures from files
sig_edgeR <- as.vector(unlist(read.delim(file = "data/edgeR_sig.txt", header = FALSE)))
sig_DeSeq2 <- as.vector(unlist(read.delim(file = "data/DESeq_sig.txt", header = FALSE)))
sig_VoomLimma <- as.vector(unlist(read.delim(file = "data/voom_sig.txt", header = FALSE)))
sig_EBSeq <- as.vector(unlist(read.delim(file = "data/EBSeq_sig.txt", header = FALSE)))

# draw a Venn diagram
draw_venn_diag(sig_edgeR, sig_DeSeq2, sig_VoomLimma, sig_EBSeq)