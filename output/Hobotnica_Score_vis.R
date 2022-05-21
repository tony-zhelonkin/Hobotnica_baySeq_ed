
library(RColorBrewer)

h_scores <- read.csv(file = "/path/hobotnica_scores.txt", header = F, sep = " ", dec = ".")

colnames(h_scores) <- c("tool", "H_score")


counts <- h_scores$H_score
colors <- brewer.pal(6, "Set2") 

barplot(counts, main="H-scores across tools",
        names.arg = h_scores$tool,
        col= colors )

png("h_scores.png", 500, 450)
barplot(counts, main="H-scores across tools",
        names.arg = h_scores$tool,
        col= colors )
dev.off()


