library(logging)
# import functions from modules
import::from(DESeq.R, deseq2)
import::from(EBSeq.R, ebseq)
import::from(edgeR.R, edger)
import::from(voom.R, voom)

# Logger config
basicConfig()

# Load data
cm_args <- commandArgs(trailingOnly = TRUE)
count <- cm_args[1]
cols <- cm_args[2]
out <- cm_args[3]
if (!file.exists(count)) {
  stop("File ", count, " does not exists!")
}
if (!file.exists(cols)) {
  stop("File ", count, " does not exists!")
}
gse <- readRDS(file = count)
coldata <- readRDS(file = cols)
condition <- coldata$condition
if (out == "") {
  stop("Enter output file as third argument!")
}

loginfo('Run 5 differential analysis tools and score with Hobotnica')
loginfo(paste0('Count matrix is in ', count))
loginfo(paste0('Coldata is in ', cols))
loginfo(paste0('Output will be written in ', out))

# Run tools
deseq2_res <- deseq2(gse, condition)
ebseq_res <- ebseq(gse, condition)
edger_res <- edger(gse, condition)
voom_res <- voom(gse, condition)

# Output
loginfo('Save results')
fout <- file(out)
writeLines(c(
  paste0("DESeq score: ", deseq2_res),
  paste0("EBSeq score: ", ebseq_res),
  paste0("edgeR score: ", edger_res),
  paste0("voom score: ", voom_res)
), fout)

