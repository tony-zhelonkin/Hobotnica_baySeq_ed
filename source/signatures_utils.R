# Make a signature of top-n genes
top_signature <- function(results_deseq2, results_ebseq, results_edger, results_voom, results_noiseq, n, out) {
    # Import top-n creating functions
    import::here(deseq2_top, .from = 'source/DESeq.R')
    import::here(ebseq_top, .from = 'source/EBSeq.R')
    import::here(edger_top, .from = 'source/edgeR.R')
    import::here(voom_top, .from = 'source/voom.R')
    import::here(noiseq_top, .from = 'source/NOISeq.R')

    #Download results of different tools
    writeLines(deseq2_top(results_deseq2, n), file.path(out, "DESeq_sig.txt"))
    writeLines(ebseq_top(results_ebseq, n), file.path(out, "EBSeq_sig.txt"))
    writeLines(edger_top(results_edger, n), file.path(out, "edgeR_sig.txt"))
    writeLines(voom_top(results_voom, n), file.path(out, "voom_sig.txt"))
    writeLines(noiseq_top(results_noiseq, n), file.path(out, "NOISeq_sig.txt"))
}

# Make a signature of filtered genes
filtered_signature <- function(results_deseq2, results_ebseq, results_edger, results_voom, results_noiseq, out) {
    # Import filtered creating functions
    import::here(deseq2_filtered, .from = 'source/DESeq.R')
    import::here(ebseq_filtered, .from = 'source/EBSeq.R')
    import::here(edger_filtered, .from = 'source/edgeR.R')
    import::here(voom_filtered, .from = 'source/voom.R')
    import::here(noiseq_filtered, .from = 'source/NOISeq.R')

    #Download results of different tools
    writeLines(deseq2_filtered(results_deseq2), file.path(out, "DESeq_sig.txt"))
    writeLines(ebseq_filtered(results_ebseq), file.path(out, "EBSeq_sig.txt"))
    writeLines(edger_filtered(results_edger), file.path(out, "edgeR_sig.txt"))
    writeLines(voom_filtered(results_voom), file.path(out, "voom_sig.txt"))
    writeLines(noiseq_filtered(results_noiseq), file.path(out, "NOISeq_sig.txt"))
}


# this function saves to file a Venn diagram based on signature of differentially expressed
# genes extracted from edgeR, DeSeq2, voom+limma, EBSeq & NOISeq instruments
draw_venn_diag <- function(out) {
    library(futile.logger)
    # Download results of signature making. Check on zero length subsets
    sig_vis <- list()
    sig_names <- c()
    number_of_elem <- 0
    if (file.info(file.path(out, "edgeR_sig.txt"))$size == 0) {
        cat('There is no gene in edgeR for that signature\n')
    } else {
        sig_vis [[length(sig_vis) + 1]] <- as.vector(unlist(read.delim(file = file.path(out, "edgeR_sig.txt"), header = FALSE)))
        sig_names <- c(sig_names, "edgeR")
        number_of_elem <- number_of_elem + 1
    }
    if (file.info(file.path(out, "DESeq_sig.txt"))$size == 0) {
        cat('There is no gene in DESeq for that signature\n')
    } else {
        sig_vis [[length(sig_vis) + 1]] <- as.vector(unlist(read.delim(file = file.path(out, "DESeq_sig.txt"), header = FALSE)))
        sig_names <- c(sig_names, "DeSeq2")
        number_of_elem <- number_of_elem + 1
    }
    if (file.info(file.path(out, "voom_sig.txt"))$size == 0) {
        cat('There is no gene in voom for that signature\n')
    } else {
        sig_vis [[length(sig_vis) + 1]] <- as.vector(unlist(read.delim(file = file.path(out, "voom_sig.txt"), header = FALSE)))
        sig_names <- c(sig_names, "Voom+Limma")
        number_of_elem <- number_of_elem + 1
    }
    if (file.info(file.path(out, "EBSeq_sig.txt"))$size == 0) {
        cat('There is no gene in EBSeq for that signature\n')
    } else {
        sig_vis [[length(sig_vis) + 1]] <- as.vector(unlist(read.delim(file = file.path(out, "EBSeq_sig.txt"), header = FALSE)))
        sig_names <- c(sig_names, "EBSeq")
        number_of_elem <- number_of_elem + 1
    }
    if (file.info(file.path(out, "NOISeq_sig.txt"))$size == 0) {
        cat('There is no gene in NOISeq for that signature\n')
    } else {
        sig_vis [[length(sig_vis) + 1]] <- as.vector(unlist(read.delim(file = file.path(out, "NOISeq_sig.txt"), header = FALSE)))
        sig_names <- c(sig_names, "NOISeq")
        number_of_elem <- number_of_elem + 1
    }


    import::here(venn.diagram, .from = VennDiagram)
    import::here(brewer.pal, .from = RColorBrewer)
    import::here(grid.newpage, grid.draw, .from = grid)
    if (number_of_elem >= 3) {
        myCol <- brewer.pal(number_of_elem, "Pastel2")   
    } else if (number_of_elem == 2) {
        myCol <- c("#B3E2CD", "#FDCDAC")
    } else if (number_of_elem == 1) {
        myCol <- c("#CBD5E8")
    }

    pdf(file.path(out, "venn_diagram.pdf"))
    grid.newpage()
    futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")
    if (number_of_elem == 5) {
        pos_ = c(0, -30, -130, 130, 30)
    } else if (number_of_elem == 4) {
        pos_ = c(-15, 15, 0, 0)
    } else if (number_of_elem == 3) {
        pos_ = c(-40, 40, 180)
    } else if (number_of_elem == 2) {
        pos_ = c(-50, 50)
    } else {
        pos_ = c(0)
    }
    venn_obj <- venn.diagram(
    x = sig_vis,
    category.names = sig_names,
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
    cat.pos = pos_,
    filename = NULL
    )
    grid.draw(venn_obj)
}
