library("BiocManager")
library("dplyr")
library("tximeta")
library("SummarizedExperiment")

library("EBSeq")
gse <- readRDS(file = "data/matrix.rds")
coldata <- readRDS(file = "data/coldata.rds")
condition <- coldata$condition
GeneMat <- assays(gse)$counts

cat ("ROWS\n")
ncol(GeneMat)
head(GeneMat, 10)
as.factor(condition)

#data(GeneMat)
#str(GeneMat)
Sizes = MedianNorm(GeneMat)
EBOut = EBTest(Data = GeneMat,
    Conditions = as.factor(condition),
    sizeFactors = Sizes, maxround = 5)
saveRDS(EBOut, file = "data/EBSeq-res.rds")
Out = GetDEResults(EBOut)

res <- Out$PPMat

head(res)
cat ("EBSeqed\n")

#Out
library("EnhancedVolcano")
GeneFC <- PostFC(EBOut, SmallNum = 0.01)
#heatmap.2(NormalizedMatrix[GenesOfInterest,], scale=”row”, trace=”none”, Colv=F)

x = log2(GeneFC$RealFC)
y = EBOut$PPDE
# Create dataframe with groups
xcond = 2
ycond = 0.7
group <- ifelse((x < -xcond | x > xcond) & y > ycond , "DE", ifelse(y > ycond , "ADE2", ifelse((x < -xcond | x > xcond), "ADE1", "NDE")))
df <- data.frame(x = x, y = y, group = factor(group))

# Change group colors
colors <- c("red", "blue", "green", "gray")
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
