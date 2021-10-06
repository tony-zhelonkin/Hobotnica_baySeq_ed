library(logging)
# Import functions from modules
import::from(DESeq.R, deseq2_f, deseq2_v)
import::from(EBSeq.R, ebseq_f, ebseq_v)
import::from(edgeR.R, edger_f, edger_v)
import::from(voom.R, voom_f, voom_v)

# Logger config
basicConfig()

# Load data
cm_args <- commandArgs(trailingOnly = TRUE)
count <- cm_args[1]
cols <- cm_args[2]
if (!file.exists(count)) {
  stop("File ", count, " does not exists!")
}
if (!file.exists(cols)) {
  stop("File ", cols, " does not exists!")
}

# Read count matrix table and correct its form
counts<- read.table(file=count, header = TRUE, sep = ",", dec = ".")
gene_names <- counts$X
counts <- counts[,-c(1)]
rownames(counts) <- gene_names

# Read annotation and make table with useful form
coldata <- read.table(file=cols, sep = ",", dec = ".")
colnames(coldata) <- c("names", "condition")
coldata <- coldata[-c(1), ]

loginfo('Run 5 differential analysis tools and score with Hobotnica')
loginfo(paste0('Count matrix is in ', count))
loginfo(paste0('Coldata is in ', cols))

# Run tools
deseq2_res <- deseq2_f(counts, coldata)
ebseq_res <- ebseq_f(counts, coldata)
edger_res <- edger_f(counts, coldata)
voom_res <- voom_f(counts, coldata)


# Save results
loginfo('Save results')
saveRDS(deseq2_res, file = "data/DESeq_res.rds")
saveRDS(ebseq_res, file = "data/EBSeq_res.rds")
saveRDS(edger_res, file = "data/edgeR_res.rds")
saveRDS(voom_res, file = "data/voom_res.rds")

# Visualize results
loginfo('Visualize results')
deseq2_v(deseq2_res)
ebseq_v(ebseq_res)
edger_v(edger_res)
voom_v(voom_res)
