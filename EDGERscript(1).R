countdata <- read.table("Normalsample.txt", header=TRUE, row.names=1)
countdata <- as.matrix(countdata)
head(countdata)
# Assign condition (first four are controls, second four contain the expansion)
(condition <- factor(c("Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal",
                       "Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal")))
library(edgeR)
d <- DGEList(counts=countdata,group=condition)
d
##(coldata <- data.frame(row.names=colnames(countdata), condition))
##dds <- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~1)
##dds
d <- calcNormFactors(d,lib.size=NULL, method = "TMM")
d

y <- as.matrix(countdata)
calcNormFactors(y)

dge <- DGEList(counts=y)
cpm(dge)
dge <- calcNormFactors(dge)
a <- cpm(dge)
write.table(a, file = "TMM_normal.txt", sep="\t", quote=F, col.names=NA)
##tumor
countdata1 <- read.table("Tumorsample.txt", header=TRUE, row.names=1)
countdata1 <- as.matrix(countdata1)
head(countdata)

y <- as.matrix(countdata1)
calcNormFactors(y)

dge <- DGEList(counts=y)
cpm(dge)
dge <- calcNormFactors(dge)
b <- cpm(dge)
write.table(b, file = "TMM_tumor.txt", sep="\t", quote=F, col.names=NA)
