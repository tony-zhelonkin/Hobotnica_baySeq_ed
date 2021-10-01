library("BiocManager")
library("dplyr")
library("tximeta")
library("SummarizedExperiment")
library("EBSeq")

counts<- read.table(file="data/TCGA_prostate_countmatrix.txt", header = TRUE, sep = ",", dec = ".")
gene_names <- counts$X
counts <- counts[,-c(1)]
rownames(counts) <- gene_names
coldata <- read.table(file="data/annotation_TCGA_prostate.txt", sep = ",", dec = ".")
colnames(coldata) <- c("names", "condition")
coldata <- coldata[-c(1), ]
condition <- coldata$condition

Sizes = MedianNorm(as.matrix(counts))

EBOut = EBTest(Data = as.matrix(counts),
    Conditions = as.factor(condition),
    sizeFactors = Sizes, maxround = 5)
saveRDS(EBOut, file = "data/EBSeq-res.rds")
Out = GetDEResults(EBOut)

res <- Out$PPMat

head(res)
cat ("EBSeqed\n")

#Visualize
library("EnhancedVolcano")
GeneFC <- PostFC(EBOut, SmallNum = 0.01)
x = log2(GeneFC$RealFC)
y = EBOut$PPDE
# Create dataframe with groups
xcond = 2
ycond = 0.7
group <- ifelse((x < -xcond | x > xcond) & y > ycond , "DE", ifelse(y > ycond , "ADE2", ifelse((x < -xcond | x > xcond), "ADE1", "NDE")))
df <- data.frame(x = x, y = y, group = factor(group))
# Change group colors
colors <- c("red", "blue", "green", "gray")
#Plotting
plot(df$x, df$y, col = colors[df$group], pch = 16,
    xlim =c(-5,5),
    ylim=c(0,1),
    main="EBSeq results",
    xlab="EBSeq Log2 Fold Change", ylab="EBSeq PPDE")
DEgplots <- data.frame(x = ifelse((x < -xcond | x > xcond) & y > ycond, x, 0), y = ifelse((x < -xcond | x > xcond) & y > ycond , y, 0),
    name = ifelse((x < -xcond | x > xcond) & y > ycond , rownames(res), ""))
text(DEgplots$x, DEgplots$y,
     labels = DEgplots$name,
     cex = 0.6, pos = 4, col = "black")
