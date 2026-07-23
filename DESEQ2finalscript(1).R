countdata <- read.table("Normalsample.txt", header=TRUE, row.names=1)
countdata <- as.matrix(countdata)
head(countdata)
# Assign condition (first four are controls, second four contain the expansion)
(condition <- factor(c("Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal",
                           "Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal")))
library(DESeq2)
# Create a coldata frame and instantiate the DESeqDataSet. See ?DESeqDataSetFromMatrix
(coldata <- data.frame(row.names=colnames(countdata), condition))
dds <- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~1)
dds
?DESeqDataSetFromMatrix
dds <- DESeq(dds)
res <- results(dds)
res
summary(res)
res2 <- lfcShrink(dds, contrast=contrast, res=res)
head(res)
## Filtering to find significant genes
padj.cutoff <- 0.05 # False Discovery Rate cutoff
significant_results <- res[which(res$padj < padj.cutoff),]
## save results using customized file_name
file_name = 'significant_padj_0.05_normal_deseq2.txt'
write.table(significant_results, file_name, quote=FALSE)
##for tumor file
countdata1 <- read.table("Tumorsample.txt", header=TRUE, row.names=1)
countdata1 <- as.matrix(countdata1)
head(countdata1)
# Assign condition (first four are controls, second four contain the expansion)
#(condition <- factor(c("Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal",
 ##                      "Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal")))
(condition <- factor(c("Tumor",	"Tumor","Tumor","Tumor","Tumor","Tumor","Tumor","Tumor","Tumor","Tumor","Tumor",
                       "Tumor","Tumor","Tumor","Tumor","Tumor","Tumor","Tumor","Tumor","Tumor","Tumor",	"Tumor",	"Tumor","Tumor","Tumor",	
                       "Tumor",	"Tumor",	"Tumor",	"Tumor",	"Tumor",	"Tumor",	"Tumor",
                       "Tumor",	"Tumor",	"Tumor",	"Tumor",	"Tumor",	"Tumor",	"Tumor",	"Tumor",
                       "Tumor",	"Tumor",	"Tumor",	"Tumor",	"Tumor",	"Tumor",	"Tumor",	"Tumor",	"Tumor",	"Tumor",	"Tumor",
                       "Tumor",	"Tumor",	"Tumor",	"Tumor",	"Tumor",	"Tumor",	"Tumor",	"Tumor",	
                       "Tumor",	"Tumor",	"Tumor",	"Tumor",	"Tumor",	"Tumor",	"Tumor",	"Tumor",	
                       "Tumor",	"Tumor",	"Tumor",	"Tumor",	"Tumor",	"Tumor",	"Tumor",	"Tumor",	"Tumor",	"Tumor",	"Tumor",	"Tumor",	
                       "Tumor",	"Tumor",	"Tumor",	"Tumor",	"Tumor",	"Tumor",	"Tumor",	"Tumor",	"Tumor",	"Tumor",	"Tumor",	"Tumor",	
                       "Tumor",	"Tumor",	"Tumor",	"Tumor",	"Tumor",	"Tumor",	"Tumor",	"Tumor",	"Tumor",	"Tumor",	"Tumor",	
                      "Tumor",	"Tumor",	"Tumor",	"Tumor",	"Tumor","Tumor","Tumor","Tumor","Tumor","Tumor","Tumor","Tumor","Tumor","Tumor",
                      "Tumor","Tumor","Tumor","Tumor","Tumor","Tumor","Tumor","Tumor","Tumor","Tumor","Tumor","Tumor","Tumor","Tumor","Tumor","Tumor",
                      "Tumor","Tumor","Tumor","Tumor","Tumor","Tumor","Tumor","Tumor","Tumor","Tumor","Tumor","Tumor","Tumor","Tumor","Tumor","Tumor",
                      "Tumor","Tumor","Tumor","Tumor","Tumor","Tumor","Tumor","Tumor","Tumor","Tumor","Tumor","Tumor","Tumor",
                      "Tumor","Tumor","Tumor","Tumor","Tumor","Tumor","Tumor","Tumor","Tumor","Tumor","Tumor","Tumor","Tumor","Tumor","Tumor",
                      "Tumor","Tumor","Tumor","Tumor","Tumor","Tumor","Tumor","Tumor","Tumor","Tumor","Tumor","Tumor","Tumor","Tumor","Tumor","Tumor",
                      "Tumor","Tumor","Tumor","Tumor","Tumor","Tumor","Tumor","Tumor","Tumor","Tumor","Tumor","Tumor","Tumor","Tumor","Tumor","Tumor",
                      "Tumor","Tumor","Tumor","Tumor","Tumor","Tumor","Tumor","Tumor","Tumor","Tumor","Tumor","Tumor","Tumor","Tumor","Tumor","Tumor",
                      "Tumor","Tumor","Tumor","Tumor","Tumor","Tumor","Tumor","Tumor","Tumor","Tumor","Tumor","Tumor","Tumor","Tumor","Tumor","Tumor",
                      "Tumor","Tumor","Tumor","Tumor","Tumor","Tumor","Tumor","Tumor","Tumor","Tumor","Tumor","Tumor","Tumor","Tumor","Tumor","Tumor",
                      "Tumor","Tumor","Tumor","Tumor","Tumor","Tumor","Tumor","Tumor","Tumor","Tumor","Tumor","Tumor","Tumor","Tumor","Tumor","Tumor",
                      "Tumor","Tumor","Tumor","Tumor","Tumor","Tumor","Tumor","Tumor","Tumor","Tumor","Tumor","Tumor","Tumor","Tumor","Tumor","Tumor",
                      "Tumor","Tumor","Tumor","Tumor","Tumor","Tumor","Tumor","Tumor","Tumor","Tumor","Tumor","Tumor","Tumor","Tumor","Tumor","Tumor",
                      "Tumor","Tumor","Tumor","Tumor","Tumor","Tumor","Tumor","Tumor","Tumor","Tumor","Tumor","Tumor","Tumor","Tumor","Tumor","Tumor",
                      "Tumor","Tumor","Tumor","Tumor","Tumor","Tumor","Tumor","Tumor","Tumor","Tumor","Tumor","Tumor","Tumor","Tumor","Tumor","Tumor",
                      "Tumor","Tumor","Tumor","Tumor","Tumor","Tumor","Tumor","Tumor","Tumor","Tumor","Tumor","Tumor","Tumor","Tumor","Tumor","Tumor",
                      "Tumor","Tumor","Tumor","Tumor","Tumor","Tumor","Tumor","Tumor","Tumor","Tumor","Tumor","Tumor","Tumor","Tumor","Tumor","Tumor",
                      "Tumor","Tumor","Tumor","Tumor","Tumor","Tumor","Tumor","Tumor","Tumor","Tumor","Tumor","Tumor","Tumor","Tumor","Tumor","Tumor",
                      "Tumor","Tumor","Tumor","Tumor","Tumor","Tumor","Tumor","Tumor","Tumor","Tumor","Tumor","Tumor","Tumor","Tumor","Tumor","Tumor",
                      "Tumor","Tumor","Tumor","Tumor","Tumor","Tumor","Tumor","Tumor","Tumor","Tumor","Tumor","Tumor","Tumor","Tumor","Tumor","Tumor","Tumor","Tumor")))
library(DESeq2)
# Create a coldata frame and instantiate the DESeqDataSet. See ?DESeqDataSetFromMatrix
(coldata <- data.frame(row.names=colnames(countdata1), condition))
dds <- DESeqDataSetFromMatrix(countData=countdata1, colData=coldata, design=~1)
dds
?DESeqDataSetFromMatrix
dds <- DESeq(dds)
res <- results(dds)
res
summary(res)
res2 <- lfcShrink(dds, contrast=contrast, res=res)
head(res)
## Filtering to find significant genes
padj.cutoff <- 0.05 # False Discovery Rate cutoff
significant_results <- res[which(res$padj < padj.cutoff),]
## save results using customized file_name
file_name = 'significant_padj_0.05_tumor_deseq2.txt'
write.table(significant_results, file_name, quote=FALSE)
dds <- estimateSizeFactors(dds)
results(dds)
##new script
##metadata <- read.table("normal_metadata.txt", header = TRUE, row.names = 1)
##metadata <- as.matrix(metadata)
##dds <- DESeqDataSetFromMatrix(countData=countdata, 
                              ##colData=metadata, 
                              ##design= ~ condition, tidy = TRUE)
##write.table(dds, file="DESeq2_normal.txt", sep = "\t", row.names = TRUE)


###SIZE-FACTOR NORMALISATION
countdata <- read.table("Normalsample.txt", header=TRUE, row.names=1)
countdata <- as.matrix(countdata)
head(countdata)
# Assign condition (first four are controls, second four contain the expansion)
(condition <- factor(c("Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal",
                       "Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal")))
library(DESeq2)
# Create a coldata frame and instantiate the DESeqDataSet. See ?DESeqDataSetFromMatrix
(coldata <- data.frame(row.names=colnames(countdata), condition))
dds <- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~1)
dds
e <- sizeFactors(dds)
write.table(e, file = "sizefactor_normal.txt")
dds <- estimateSizeFactors(dds)
normalized_counts <- counts(dds, normalized=TRUE)
write.table(normalized_counts, file="sizefactorTUMOR.txt", sep="\t", quote=F, col.names=NA)
