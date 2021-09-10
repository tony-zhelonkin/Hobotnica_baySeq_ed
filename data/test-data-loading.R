#Find data
sample_name <- "SRR1039"

#Create data array
names <- c()
for (i in "508":"523") {
    names <- c(names, paste(sample_name, i, sep = "", collapse = ""))
}

condition<- c("Untreated", "Treated", "Untreated", "Treated", "Untreated", "Treated", "Untreated", "Treated", "Untreated", "Treated", "Untreated", "Treated", "Untreated", "Treated", "Untreated", "Treated")
coldata <- data.frame(names = names, condition = condition, stringsAsFactors=FALSE)
saveRDS(coldata, file = "coldata.rds")

