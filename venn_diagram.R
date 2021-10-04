library("VennDiagram")
library("RColorBrewer")

# this function saves to file a Venn diagram based on top-20 differentially expressed 
# genes extracted from edgeR, DeSeq2, voom+limma & EBSeq instruments
draw_venn_diag <- function(top_edgeR, top_DeSeq2, top_VoomLimma, top_EBSeq) {
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
    #filename = NULL
    filename = "data/venn_diag.png"
  )
  grid.draw(venn_obj)
}

# read top-20 differentially expressed genes from files
top_edgeR <- as.vector(unlist(read.delim(file = "data/edgeR-top.txt", header = FALSE)))
top_DeSeq2 <- as.vector(unlist(read.delim(file = "data/DESeq-top.txt", header = FALSE)))
top_VoomLimma <- as.vector(unlist(read.delim(file = "data/voom-top.txt", header = FALSE)))
top_EBSeq <- as.vector(unlist(read.delim(file = "data/EBSeq-top.txt", header = FALSE)))

# draw a Venn diagram
draw_venn_diag(top_edgeR, top_DeSeq2, top_VoomLimma, top_EBSeq)