workdir=getwd()
setwd(workdir)

gene_count=read.delim(file = "gene_expression.txt",header = T,sep = "\t",comment.char = "#",stringsAsFactors = F)
head(gene_count)
rownames(gene_count)=gene_count$Geneid
as.matrix(names(gene_count))
#gene_count=gene_count[,c(1,7:ncol(gene_count))]
gene_count=gene_count[,c(1,7:ncol(gene_count))]
#gene_count=gene_count[,c(1,7,10)]
as.matrix(names(gene_count))
head(gene_count)
##colnames(gene_count)=c("Id","Control","Test")
##as.matrix(names(gene_count))
TMP=gene_count
TMP=TMP[,-1]
head(TMP)
require(devtools)
library(devtools)
lapply(names(TMP),function(x){write.table(cbind(rownames(TMP),TMP[[x]]), file = paste("output", x, sep = ""),row.names = F, quote = F, col.names = F,sep = "\t")})
getwd()

##For DESeq analysis, we need a matrix as input so,transform countdata data frame into matrix of counts with transcripts ids as its row names.##
y=as.matrix(gene_count[,c(2:ncol(gene_count))])
head(y)
rownames(y)=gene_count$Geneid
head(y)
#storage.mode(y) = "integer"
##for DESEq2 analysis
y1=as.data.frame(y[,c(1,2)])
head(y1)
y2=as.data.frame(y[,c(1,5)])
head(y2)
y3=as.data.frame(y[,c(2,6)])
head(y3)
y4=as.data.frame(y[,c(5,6)])
head(y4)
y5=as.data.frame(y[,c(3,4)])
head(y5)
y6=as.data.frame(y[,c(3,7)])
head(y6)
y7=as.data.frame(y[,c(4,8)])
head(y7)
y8=as.data.frame(y[,c(7,8)])
head(y8)
##set the column name to condition##
##condition=factor(c("Control", "Patient" , "Patient" , "Patient" , "Patient"))
condition=factor(c("SAMPLE1" ,"SAMPLE2","SAMPLE3","SAMPLE4","SAMPLE5","SAMPLE6","SAMPLE7","SAMPLE8"))
condition1=factor(c("SAMPLE1" ,"SAMPLE2"))
condition2=factor(c("SAMPLE1","SAMPLE5"))
condition3=factor(c("SAMPLE2","SAMPLE6"))
condition4=factor(c("SAMPLE5","SAMPLE6"))
condition5=factor(c("SAMPLE3","SAMPLE4"))
condition6=factor(c("SAMPLE3","SAMPLE7"))
condition7=factor(c("SAMPLE4","SAMPLE8"))
condition8=factor(c("SAMPLE7","SAMPLE8"))

#batch=factor(c("Test_Control1" , "Test_Control2" , "Control2" , "Control2"))

##############create DESeqDataSet###########
library(DESeq2)
##require(DESeq2)
##dde=DESeqDataSetFromMatrix(countData = y, colData=DataFrame(condition),design = ~condition)
dde=DESeqDataSetFromMatrix(countData = y,colData = DataFrame(condition),design = ~condition)
dde1=DESeqDataSetFromMatrix(countData = y1,colData = DataFrame(condition1),design = ~condition1)
dde2=DESeqDataSetFromMatrix(countData = y2,colData = DataFrame(condition2),design = ~condition2)
dde3=DESeqDataSetFromMatrix(countData = y3,colData = DataFrame(condition3),design = ~condition3)
dde4=DESeqDataSetFromMatrix(countData = y4,colData = DataFrame(condition4),design = ~condition4)
dde5=DESeqDataSetFromMatrix(countData = y5,colData = DataFrame(condition5),design = ~condition5)
dde6=DESeqDataSetFromMatrix(countData = y6,colData = DataFrame(condition6),design = ~condition6)
dde7=DESeqDataSetFromMatrix(countData = y7,colData = DataFrame(condition7),design = ~condition7)
dde8=DESeqDataSetFromMatrix(countData = y8,colData = DataFrame(condition8),design = ~condition8)

#dde=DESeqDataSetFromMatrix(countData = y, colData=DataFrame(condition,batch),design = ~batch + condition)

##Retain genes which have at least 1 read in any of the samples##
dde=estimateSizeFactors(dde)
sizeFactors(dde)
dde1=estimateSizeFactors(dde1)
sizeFactors(dde1)
dde2=estimateSizeFactors(dde2)
sizeFactors(dde2)
dde3=estimateSizeFactors(dde3)
sizeFactors(dde3)
dde4=estimateSizeFactors(dde4)
sizeFactors(dde4)
dde5=estimateSizeFactors(dde5)
sizeFactors(dde5)
dde6=estimateSizeFactors(dde6)
sizeFactors(dde6)
dde7=estimateSizeFactors(dde7)
sizeFactors(dde7)
dde8=estimateSizeFactors(dde8)
sizeFactors(dde8)
#       AB-T1        AB-T2   AB-Control CL-7-Control      CL-7-T1      CL-7-T2      CL-7-T3      CL-7-T4 
#0.4079209    1.6541158    1.4671952    1.7484639    1.1823330    0.6441425    1.4242791    1.4153648 
#CL-9-Control    CL-9-Test 
#1.2409271    0.2822171 


#dds=dde[rowSums(counts(dde)>5)>=3]
dds=dde[rowSums(counts(dde))>3]
dds1=dde1[rowSums(counts(dde1))>3]
dds2=dde2[rowSums(counts(dde2))>3]
dds3=dde3[rowSums(counts(dde3))>3]
dds4=dde4[rowSums(counts(dde4))>3]
dds5=dde5[rowSums(counts(dde5))>3]
dds6=dde6[rowSums(counts(dde6))>3]
dds7=dde7[rowSums(counts(dde7))>3]
dds8=dde8[rowSums(counts(dde8))>3]

##dim: 4673 2 

##Now specify the condition that should be considered as reference for differential expression analysis.##
#dds$condition<- relevel(dds$condition, ref="VectorControl_Elute")
dds1$condition1= relevel(dds1$condition1, ref="SAMPLE1")
dds2$condition2= relevel(dds2$condition2, ref="SAMPLE1")
dds3$condition3= relevel(dds3$condition3, ref="SAMPLE2")
dds4$condition4= relevel(dds4$condition4, ref="SAMPLE5")
dds5$condition5= relevel(dds5$condition5, ref="SAMPLE3")
dds6$condition6= relevel(dds6$condition6, ref="SAMPLE3")
dds7$condition7= relevel(dds7$condition7, ref="SAMPLE4")
dds8$condition8= relevel(dds8$condition8, ref="SAMPLE7")
##Perform Differential expression (DE) analysis##
dds1=DESeq(dds1)
resultsNames(dds1)
dds2=DESeq(dds2)
resultsNames(dds2)
dds3=DESeq(dds3)
resultsNames(dds3)
dds4=DESeq(dds4)
resultsNames(dds4)
dds5=DESeq(dds5)
resultsNames(dds5)
dds6=DESeq(dds6)
resultsNames(dds6)
dds7=DESeq(dds7)
resultsNames(dds7)
dds8=DESeq(dds8)
resultsNames(dds8)

normcounts=counts(dde, normalized=TRUE)
normcounts1=counts(dde1, normalized=TRUE)
normcounts2=counts(dde2, normalized=TRUE)
normcounts3=counts(dde3, normalized=TRUE)
normcounts4=counts(dde4, normalized=TRUE)
normcounts5=counts(dde5, normalized=TRUE)
normcounts6=counts(dde6, normalized=TRUE)
normcounts7=counts(dde7, normalized=TRUE)
normcounts8=counts(dde8, normalized=TRUE)
#vsd <- vst(dds)
#plotPCA(vsd, "condition")
#assay(vsd) <- limma::removeBatchEffect(assay(vsd), vsd$batch)
#png("PCA_plot.png", type="cairo", res = 300, units="in", width=7, height=4)
#plotPCA(vsd, "condition")
#dev.off()

#normcounts=assay(vsd)
#plotPCA(vsd, "condition")
#assay(vsd) <- limma::removeBatchEffect(assay(vsd), vsd$batch)
#png("PCA_plot.png", type="cairo", res = 300, units="in", width=8, height=8)
#plotPCA(vsd, "condition")
#dev.off()

#normcounts=assay(dde)
#write.table(as.data.frame(normcounts), file="Normalized_counts.csv",sep = "\t",row.names = T,col.names = T, quote = F)

#calulate corrilation matrix
#corr_matrix=cor(forpca, method = "spearman")
corr_matrix=cor(normcounts, method = "spearman")
openxlsx::write.xlsx(corr_matrix, file = "03_Correlation_matrix.xlsx")

require(PerformanceAnalytics)
png("Correlation_plot.png", type="cairo", res = 600, units="in", width=ncol(corr_matrix), height=nrow(corr_matrix))
chart.Correlation(corr_matrix, histogram=TRUE, pch=19,method = "spearman")
dev.off()
detach("package:PerformanceAnalytics")

#Import target file with grouping information
require(readr)
target=read_delim("targetfile.txt","\t", escape_double = FALSE, trim_ws = TRUE)

require(ggfortify)
require(ggplot2)
require(ggrepel)
#Draw PCA plot
tpca=t(normcounts)
tpca=data.frame(tpca)
tpca$group=target$label
png("PCA_plot.png", type="cairo", res = 300, units="in", width=14, height=10)
autoplot(prcomp(tpca[,c(1:(ncol(tpca)-1))]), data=tpca, colour='group')+
  geom_text_repel(aes(label=rownames(tpca),color=factor(group)),size=6)
#autoplot(prcomp(tpca[,c(1:(ncol(tpca)-1))]), data=tpca, colour='group', label.size = 5) +
#+  geom_text_repel(aes(label = rownames(prcomp(tpca))))
dev.off()


# ##### Annotation to add with files ########
###Annotation
gff_anno=read.delim("GRCh38_gene_anno_all.New2021_gene.txt",stringsAsFactors = F, header = T, sep = "\t")
head(gff_anno)
#RAW COUNTS AND NORMALIZED COUNTS
norm_counts=as.data.frame(normcounts)
norm_counts$Id <- rownames(norm_counts)
head(norm_counts)
norm_counts=norm_counts[,c(ncol(norm_counts),1:(ncol(norm_counts)-1))]
norm_counts.anno=merge(norm_counts,gff_anno,by.x="Id",by.y="gene.id",all.x=T)
openxlsx::write.xlsx(norm_counts.anno, file = "01_Normalized_Expression_Counts.xlsx")
gene_count.anno=merge(gene_count,gff_anno,by.x="Geneid",by.y="Id",all.x=T)
openxlsx::write.xlsx(gene_count.anno, file = "02_Expression_Counts.xlsx")

norm_counts1=as.data.frame(normcounts1)
norm_counts1$Id <- rownames(norm_counts1)
head(norm_counts1)
norm_counts1=norm_counts1[,c(ncol(norm_counts1),1:(ncol(norm_counts1)-1))]
norm_counts1.anno=merge(norm_counts1,gff_anno,by.x="Id",by.y="Id",all.x=T)
#openxlsx::write.xlsx(norm_counts1.anno, file = "01_Normalized_Expression_Counts.xlsx")
gene_count.anno=merge(gene_count,gff_anno,by.x="Geneid",by.y="Id",all.x=T)
#openxlsx::write.xlsx(gene_count.anno, file = "02_Expression_Counts.xlsx")

norm_counts2=as.data.frame(normcounts2)
norm_counts2$Id <- rownames(norm_counts2)
head(norm_counts2)
norm_counts2=norm_counts2[,c(ncol(norm_counts2),1:(ncol(norm_counts2)-1))]
norm_counts2.anno=merge(norm_counts2,gff_anno,by.x="Id",by.y="Id",all.x=T)
#openxlsx::write.xlsx(norm_counts2.anno, file = "01_Normalized_Expression_Counts.xlsx")
gene_count.anno=merge(gene_count,gff_anno,by.x="Geneid",by.y="Id",all.x=T)
#openxlsx::write.xlsx(gene_count.anno, file = "02_Expression_Counts.xlsx")

norm_counts3=as.data.frame(normcounts3)
norm_counts3$Id <- rownames(norm_counts3)
head(norm_counts3)
norm_counts3=norm_counts3[,c(ncol(norm_counts3),1:(ncol(norm_counts3)-1))]
norm_counts3.anno=merge(norm_counts3,gff_anno,by.x="Id",by.y="gene_id",all.x=T)
#openxlsx::write.xlsx(norm_counts3.anno, file = "01_Normalized_Expression_Counts.xlsx")
gene_count.anno=merge(gene_count,gff_anno,by.x="Geneid",by.y="Id",all.x=T)
#openxlsx::write.xlsx(gene_count.anno, file = "02_Expression_Counts.xlsx")

norm_counts4=as.data.frame(normcounts4)
norm_counts4$Id <- rownames(norm_counts4)
head(norm_counts4)
norm_counts4=norm_counts4[,c(ncol(norm_counts4),1:(ncol(norm_counts4)-1))]
norm_counts4.anno=merge(norm_counts4,gff_anno,by.x="Id",by.y="gene_id",all.x=T)
#openxlsx::write.xlsx(norm_counts4.anno, file = "01_Normalized_Expression_Counts.xlsx")
gene_count.anno=merge(gene_count,gff_anno,by.x="Geneid",by.y="Id",all.x=T)
#openxlsx::write.xlsx(gene_count.anno, file = "02_Expression_Counts.xlsx")

norm_counts5=as.data.frame(normcounts5)
norm_counts5$Id <- rownames(norm_counts5)
head(norm_counts5)
norm_counts5=norm_counts5[,c(ncol(norm_counts5),1:(ncol(norm_counts5)-1))]
norm_counts5.anno=merge(norm_counts5,gff_anno,by.x="Id",by.y="Id",all.x=T)
#openxlsx::write.xlsx(norm_counts5.anno, file = "01_Normalized_Expression_Counts.xlsx")
gene_count.anno=merge(gene_count,gff_anno,by.x="Geneid",by.y="Id",all.x=T)
#openxlsx::write.xlsx(gene_count.anno, file = "02_Expression_Counts.xlsx")

norm_counts6=as.data.frame(normcounts6)
norm_counts6$Id <- rownames(norm_counts6)
head(norm_counts6)
norm_counts6=norm_counts6[,c(ncol(norm_counts6),1:(ncol(norm_counts6)-1))]
norm_counts6.anno=merge(norm_counts6,gff_anno,by.x="Id",by.y="Id",all.x=T)
#openxlsx::write.xlsx(norm_counts6.anno, file = "01_Normalized_Expression_Counts.xlsx")
gene_count.anno=merge(gene_count,gff_anno,by.x="Geneid",by.y="Id",all.x=T)
#openxlsx::write.xlsx(gene_count.anno, file = "02_Expression_Counts.xlsx")

norm_counts7=as.data.frame(normcounts7)
norm_counts7$Id <- rownames(norm_counts7)
head(norm_counts7)
norm_counts7=norm_counts7[,c(ncol(norm_counts7),1:(ncol(norm_counts7)-1))]
norm_counts7.anno=merge(norm_counts7,gff_anno,by.x="Id",by.y="Id",all.x=T)
#openxlsx::write.xlsx(norm_counts7.anno, file = "01_Normalized_Expression_Counts.xlsx")
gene_count.anno=merge(gene_count,gff_anno,by.x="Geneid",by.y="Id",all.x=T)
#openxlsx::write.xlsx(gene_count.anno, file = "02_Expression_Counts.xlsx")

norm_counts8=as.data.frame(normcounts8)
norm_counts8$Id <- rownames(norm_counts8)
head(norm_counts8)
norm_counts8=norm_counts8[,c(ncol(norm_counts8),1:(ncol(norm_counts8)-1))]
norm_counts8.anno=merge(norm_counts8,gff_anno,by.x="Id",by.y="Id",all.x=T)
#openxlsx::write.xlsx(norm_counts8.anno, file = "01_Normalized_Expression_Counts.xlsx")
gene_count.anno=merge(gene_count,gff_anno,by.x="Geneid",by.y="Id",all.x=T)
#openxlsx::write.xlsx(gene_count.anno, file = "02_Expression_Counts.xlsx")


###############################################
######## Controlvs Test. #######
###############################################
#PUT Control IN SECOND PLACE like AB-Control
MSWT3_MSGN2=results(dds1, contrast=c("condition1","MSGN2","MSWT3"))
MSWT3_MSGN2_DE1=as.data.frame(MSWT3_MSGN2)
as.matrix(names(MSWT3_MSGN2_DE1))
MSWT3_MSGN2_DE1=merge(MSWT3_MSGN2_DE1,gff_anno,by.x="row.names",by.y="gene_id",all.x=T)
colnames(MSWT3_MSGN2_DE1) <- sub(pattern = "Row.names",replacement = "Id", colnames(MSWT3_MSGN2_DE1))

MSLT2_MS202=results(dds2, contrast=c("condition2","MS202_3","MSLT2_1"))
MSLT2_MS202_DE2=as.data.frame(MSLT2_MS202)
as.matrix(names(MSLT2_MS202_DE2))
MSLT2_MS202_DE2=merge(MSLT2_MS202_DE2,gff_anno,by.x="row.names",by.y="gene_id",all.x=T)
colnames(MSLT2_MS202_DE2) <- sub(pattern = "Row.names",replacement = "Id", colnames(MSLT2_MS202_DE2))

MSLT2_MS204=results(dds3, contrast=c("condition3","MS204_3","MSLT2_1"))
MSLT2_MS204_DE3=as.data.frame(MSLT2_MS204)
as.matrix(names(MSLT2_MS204))
MSLT2_MS204_DE3=merge(MSLT2_MS204_DE3,gff_anno,by.x="row.names",by.y="gene_id",all.x=T)
colnames(MSLT2_MS204_DE3) <- sub(pattern = "Row.names",replacement = "Id", colnames(MSLT2_MS204_DE3))

MS202_MS204=results(dds4, contrast=c("condition4","MS204_3","MS202_3"))
MS202_MS204_DE4=as.data.frame(MS202_MS204)
as.matrix(names(MS202_MS204_DE4))
MS202_MS204_DE4=merge(MS202_MS204_DE4,gff_anno,by.x="row.names",by.y="gene_id",all.x=T)
colnames(MS202_MS204_DE4) <- sub(pattern = "Row.names",replacement = "Id", colnames(MS202_MS204_DE4))

# Control_vs_Test5=results(dds5, contrast=c("condition5","SAMPLE4","SAMPLE3"))
# Control_vs_Test_DE5=as.data.frame(Control_vs_Test5)
# as.matrix(names(Control_vs_Test_DE5))
# Control_vs_Test_DE5=merge(Control_vs_Test_DE5,gff_anno,by.x="row.names",by.y="Id",all.x=T)
# colnames(Control_vs_Test_DE5) <- sub(pattern = "Row.names",replacement = "Id", colnames(Control_vs_Test_DE5))
# 
# Control_vs_Test6=results(dds6, contrast=c("condition6","SAMPLE7","SAMPLE3"))
# Control_vs_Test_DE6=as.data.frame(Control_vs_Test6)
# as.matrix(names(Control_vs_Test_DE6))
# Control_vs_Test_DE6=merge(Control_vs_Test_DE6,gff_anno,by.x="row.names",by.y="Id",all.x=T)
# colnames(Control_vs_Test_DE6) <- sub(pattern = "Row.names",replacement = "Id", colnames(Control_vs_Test_DE6))
# 
# Control_vs_Test7=results(dds7, contrast=c("condition7","SAMPLE8","SAMPLE4"))
# Control_vs_Test_DE7=as.data.frame(Control_vs_Test7)
# as.matrix(names(Control_vs_Test_DE7))
# Control_vs_Test_DE7=merge(Control_vs_Test_DE7,gff_anno,by.x="row.names",by.y="Id",all.x=T)
# colnames(Control_vs_Test_DE7) <- sub(pattern = "Row.names",replacement = "Id", colnames(Control_vs_Test_DE7))
# 
# Control_vs_Test8=results(dds8, contrast=c("condition8","SAMPLE8","SAMPLE7"))
# Control_vs_Test_DE8=as.data.frame(Control_vs_Test8)
# as.matrix(names(Control_vs_Test_DE8))
# Control_vs_Test_DE8=merge(Control_vs_Test_DE8,gff_anno,by.x="row.names",by.y="Id",all.x=T)
# colnames(Control_vs_Test_DE8) <- sub(pattern = "Row.names",replacement = "Id", colnames(Control_vs_Test_DE8))

##add fold change
as.matrix(names(MSWT3_MSGN2_DE1))
MSWT3_MSGN2_DE1$FoldChange=2^(MSWT3_MSGN2_DE1$log2FoldChange)
MSWT3_MSGN2_DE1$FoldChange=format(round(MSWT3_MSGN2_DE1$FoldChange, 5), nsmall = 3)
MSWT3_MSGN2_DE1=MSWT3_MSGN2_DE1[,c(1,2,ncol(MSWT3_MSGN2_DE1),3:(ncol(MSWT3_MSGN2_DE1)-1))]
as.matrix(names(MSWT3_MSGN2_DE1))

as.matrix(names(MSLT2_MS202_DE2))
MSLT2_MS202_DE2$FoldChange=2^(MSLT2_MS202_DE2$log2FoldChange)
MSLT2_MS202_DE2$FoldChange=format(round(MSLT2_MS202_DE2$FoldChange, 5), nsmall = 3)
MSLT2_MS202_DE2=MSLT2_MS202_DE2[,c(1,2,ncol(MSLT2_MS202_DE2),3:(ncol(MSLT2_MS202_DE2)-1))]
as.matrix(names(MSLT2_MS202_DE2))

as.matrix(names(MSLT2_MS204_DE3))
MSLT2_MS204_DE3$FoldChange=2^(MSLT2_MS204_DE3$log2FoldChange)
MSLT2_MS204_DE3$FoldChange=format(round(MSLT2_MS204_DE3$FoldChange, 5), nsmall = 3)
MSLT2_MS204_DE3=MSLT2_MS204_DE3[,c(1,2,ncol(MSLT2_MS204_DE3),3:(ncol(MSLT2_MS204_DE3)-1))]
as.matrix(names(MSLT2_MS204_DE3))

as.matrix(names(MS202_MS204_DE4))
MS202_MS204_DE4$FoldChange=2^(MS202_MS204_DE4$log2FoldChange)
MS202_MS204_DE4$FoldChange=format(round(MS202_MS204_DE4$FoldChange, 5), nsmall = 3)
MS202_MS204_DE4=MS202_MS204_DE4[,c(1,2,ncol(MS202_MS204_DE4),3:(ncol(MS202_MS204_DE4)-1))]
as.matrix(names(MS202_MS204_DE4))

# as.matrix(names(Control_vs_Test_DE5))
# Control_vs_Test_DE5$FoldChange=2^(Control_vs_Test_DE5$log2FoldChange)
# Control_vs_Test_DE5$FoldChange=format(round(Control_vs_Test_DE5$FoldChange, 5), nsmall = 3)
# Control_vs_Test_DE5=Control_vs_Test_DE5[,c(1,2,ncol(Control_vs_Test_DE5),3:(ncol(Control_vs_Test_DE5)-1))]
# as.matrix(names(Control_vs_Test_DE5))
# 
# as.matrix(names(Control_vs_Test_DE6))
# Control_vs_Test_DE6$FoldChange=2^(Control_vs_Test_DE6$log2FoldChange)
# Control_vs_Test_DE6$FoldChange=format(round(Control_vs_Test_DE6$FoldChange, 5), nsmall = 3)
# Control_vs_Test_DE6=Control_vs_Test_DE6[,c(1,2,ncol(Control_vs_Test_DE6),3:(ncol(Control_vs_Test_DE6)-1))]
# as.matrix(names(Control_vs_Test_DE6))
# 
# as.matrix(names(Control_vs_Test_DE7))
# Control_vs_Test_DE7$FoldChange=2^(Control_vs_Test_DE7$log2FoldChange)
# Control_vs_Test_DE7$FoldChange=format(round(Control_vs_Test_DE7$FoldChange, 5), nsmall = 3)
# Control_vs_Test_DE7=Control_vs_Test_DE7[,c(1,2,ncol(Control_vs_Test_DE7),3:(ncol(Control_vs_Test_DE7)-1))]
# as.matrix(names(Control_vs_Test_DE7))
# 
# as.matrix(names(Control_vs_Test_DE8))
# Control_vs_Test_DE8$FoldChange=2^(Control_vs_Test_DE8$log2FoldChange)
# Control_vs_Test_DE8$FoldChange=format(round(Control_vs_Test_DE8$FoldChange, 5), nsmall = 3)
# Control_vs_Test_DE8=Control_vs_Test_DE8[,c(1,2,ncol(Control_vs_Test_DE8),3:(ncol(Control_vs_Test_DE8)-1))]
# as.matrix(names(Control_vs_Test_DE8))


####add normalize counts and mean of replicates###
MSWT3_MSGN2_DE1=merge(MSWT3_MSGN2_DE1,norm_counts1,by = "Id",all.x= "T")
as.matrix(names(MSWT3_MSGN2_DE1))

MSLT2_MS202_DE2=merge(MSLT2_MS202_DE2,norm_counts2,by = "Id",all.x= "T")
as.matrix(names(MSLT2_MS202_DE2))

MSLT2_MS204_DE3=merge(MSLT2_MS204_DE3,norm_counts3,by = "Id",all.x= "T")
as.matrix(names(MSLT2_MS204_DE3))

MS202_MS204_DE4=merge(MS202_MS204_DE4,norm_counts4,by = "Id",all.x= "T")
as.matrix(names(MS202_MS204_DE4))

# Control_vs_Test_DE5=merge(Control_vs_Test_DE5,norm_counts5,by = "Id",all.x= "T")
# as.matrix(names(Control_vs_Test_DE5))
# 
# Control_vs_Test_DE6=merge(Control_vs_Test_DE6,norm_counts6,by = "Id",all.x= "T")
# as.matrix(names(Control_vs_Test_DE6))
# 
# Control_vs_Test_DE7=merge(Control_vs_Test_DE7,norm_counts7,by = "Id",all.x= "T")
# as.matrix(names(Control_vs_Test_DE7))
# 
# Control_vs_Test_DE8=merge(Control_vs_Test_DE8,norm_counts8,by = "Id",all.x= "T")
# as.matrix(names(Control_vs_Test_DE8))
#write.table(Control_vs_Test_DE1, file = "DE_Test_Control_Control_All", col.names = TRUE, sep = "\t")
MSWT3_MSGN2_DE1=MSWT3_MSGN2_DE1[,c(1,12,13,2:11)]
as.matrix(names(MSWT3_MSGN2_DE1))

MSLT2_MS202_DE2=MSLT2_MS202_DE2[,c(1,12,13,2:11)]
as.matrix(names(MSLT2_MS202_DE2))

MSLT2_MS204_DE3=MSLT2_MS204_DE3[,c(1,12,13,2:11)]
as.matrix(names(MSLT2_MS204_DE3))

MS202_MS204_DE4=MS202_MS204_DE4[,c(1,12,13,2:11)]
as.matrix(names(MS202_MS204_DE4))

# Control_vs_Test_DE5=Control_vs_Test_DE5[,c(1,13,14,2:12)]
# as.matrix(names(Control_vs_Test_DE5))
# 
# Control_vs_Test_DE6=Control_vs_Test_DE6[,c(1,13,14,2:12)]
# as.matrix(names(Control_vs_Test_DE6))
# 
# Control_vs_Test_DE7=Control_vs_Test_DE7[,c(1,13,14,2:12)]
# as.matrix(names(Control_vs_Test_DE7))
# 
# Control_vs_Test_DE8=Control_vs_Test_DE8[,c(1,13,14,2:12)]
# as.matrix(names(Control_vs_Test_DE8))


#####################################################################################################
####################################################################################################
dir.create("MS202_MS204")
openxlsx::write.xlsx(MSWT3_MSGN2_DE1, file = "MSWT3_MSGN2/MSWT3_vs_MSGN2.Differential_Expression.xlsx")
openxlsx::write.xlsx(MSLT2_MS202_DE2, file = "MSLT2_MS202/MSLT2_vs_MS202.Differential_Expression.xlsx")
openxlsx::write.xlsx(MSLT2_MS204_DE3, file = "MSLT2_MS204/MSLT2_vs_MS204.Differential_Expression.xlsx")
openxlsx::write.xlsx(MS202_MS204_DE4, file = "MS202_MS204/MS202_vs_MS204.Differential_Expression.xlsx")
# openxlsx::write.xlsx(Control_vs_Test_DE5, file = "Control_vs_Test/Sample3_vs_Sample4.Differential_Expression.xlsx")
# openxlsx::write.xlsx(Control_vs_Test_DE6, file = "Control_vs_Test/Sample3_vs_Sample7.Differential_Expression.xlsx")
# openxlsx::write.xlsx(Control_vs_Test_DE7, file = "Control_vs_Test/Sample4_vs_Sample8.Differential_Expression.xlsx")
# openxlsx::write.xlsx(Control_vs_Test_DE8, file = "Control_vs_Test/Sample7_vs_Sample8.Differential_Expression.xlsx")
#Significant DiffExprs
MSWT3_MSGN2_DE.sig1=subset(MSWT3_MSGN2_DE1, abs(log2FoldChange)>=1 & pvalue <= 0.05)
dim(MSWT3_MSGN2_DE.sig1)
openxlsx::write.xlsx(MSWT3_MSGN2_DE.sig1, file = "MSWT3_MSGN2/MSWT3_vs_MSGN2.Significant_1.xlsx")
as.matrix(names(MSWT3_MSGN2_DE.sig1))

MSLT2_MS202_DE.sig2=subset(MSLT2_MS202_DE2, abs(log2FoldChange)>=1 & pvalue <= 0.05)
dim(MSLT2_MS202_DE.sig2)
openxlsx::write.xlsx(MSLT2_MS202_DE.sig2, file = "MSLT2_MS202/MSLT2_vs_M202.Significant_1.xlsx")
as.matrix(names(MSLT2_MS202_DE.sig2))

MSLT2_MS204_DE.sig3=subset(MSLT2_MS204_DE3, abs(log2FoldChange)>=1 & pvalue <= 0.05)
dim(MSLT2_MS204_DE.sig3)
openxlsx::write.xlsx(MSLT2_MS204_DE.sig3, file = "MSLT2_MS204/MSLT2_vs_MS204.Significant_1.xlsx")
as.matrix(names(MSLT2_MS204_DE.sig3))

MS202_MS204_DE.sig4=subset(MS202_MS204_DE4, abs(log2FoldChange)>=1 & pvalue <= 0.05)
dim(MS202_MS204_DE.sig4)
openxlsx::write.xlsx(MS202_MS204_DE.sig4, file = "MS202_MS204/MS202_vs_MS204.Significant_1.xlsx")
as.matrix(names(MS202_MS204_DE.sig4))

# Control_vs_Test_DE.sig5=subset(Control_vs_Test_DE5, abs(log2FoldChange)>=1 & pvalue <= 0.05)
# dim(Control_vs_Test_DE.sig5)
# openxlsx::write.xlsx(Control_vs_Test_DE.sig5, file = "Control_vs_Test/Sample3_vs_Sample4.Significant_1.xlsx")
# as.matrix(names(Control_vs_Test_DE.sig5))
# 
# Control_vs_Test_DE.sig6=subset(Control_vs_Test_DE6, abs(log2FoldChange)>=1 & pvalue <= 0.05)
# dim(Control_vs_Test_DE.sig6)
# openxlsx::write.xlsx(Control_vs_Test_DE.sig6, file = "Control_vs_Test/Sample3_vs_Sample7.Significant_1.xlsx")
# as.matrix(names(Control_vs_Test_DE.sig6))
# 
# Control_vs_Test_DE.sig7=subset(Control_vs_Test_DE7, abs(log2FoldChange)>=1 & pvalue <= 0.05)
# dim(Control_vs_Test_DE.sig7)
# openxlsx::write.xlsx(Control_vs_Test_DE.sig7, file = "Control_vs_Test/Sample4_vs_Sample8.Significant_1.xlsx")
# as.matrix(names(Control_vs_Test_DE.sig7))
# 
# Control_vs_Test_DE.sig8=subset(Control_vs_Test_DE8, abs(log2FoldChange)>=1 & pvalue <= 0.05)
# dim(Control_vs_Test_DE.sig8)
# openxlsx::write.xlsx(Control_vs_Test_DE.sig8, file = "Control_vs_Test/Sample7_vs_Sample8.Significant_1.xlsx")
# as.matrix(names(Control_vs_Test_DE.sig8))

#Over and Under expressed transcripts
MSWT3_MSGN2_up1=subset(MSWT3_MSGN2_DE1, MSWT3_MSGN2_DE1$log2FoldChange >= 1 & MSWT3_MSGN2_DE1$pvalue<= 0.05)
dim(MSWT3_MSGN2_up1)
openxlsx::write.xlsx(MSWT3_MSGN2_up1, file = "MSWT3_MSGN2/MSWT3_vs_MSGN2.Significant_Up_1.xlsx")

MSLT2_MS202_up2=subset(MSLT2_MS202_DE2, MSLT2_MS202_DE2$log2FoldChange >= 1 & MSLT2_MS202_DE2$pvalue<= 0.05)
dim(MSLT2_MS202_up2)
openxlsx::write.xlsx(MSLT2_MS202_up2, file = "MSLT2_MS202/MSLT2_vs_MS202.Significant_Up_1.xlsx")

MSLT2_MS204_up3=subset(MSLT2_MS204_DE3, MSLT2_MS204_DE3$log2FoldChange >= 1 & MSLT2_MS204_DE3$pvalue<= 0.05)
dim(MSLT2_MS204_up3)
openxlsx::write.xlsx(MSLT2_MS204_up3, file = "MSLT2_MS204/MSLT2_vs_MS204.Significant_Up_1.xlsx")

MS202_MS204_up4=subset(MS202_MS204_DE4, MS202_MS204_DE4$log2FoldChange >= 1 & MS202_MS204_DE4$pvalue<= 0.05)
dim(MS202_MS204_up4)
openxlsx::write.xlsx(MS202_MS204_up4, file = "MS202_MS204/MS202_vs_MS204.Significant_Up_1.xlsx")

# Control_vs_Test_up5=subset(Control_vs_Test_DE5, Control_vs_Test_DE5$log2FoldChange >= 1 & Control_vs_Test_DE5$pvalue<= 0.05)
# dim(Control_vs_Test_up5)
# openxlsx::write.xlsx(Control_vs_Test_up5, file = "Control_vs_Test/Sample3_vs_Sample4.Significant_Up_1.xlsx")
# 
# Control_vs_Test_up6=subset(Control_vs_Test_DE6, Control_vs_Test_DE6$log2FoldChange >= 1 & Control_vs_Test_DE6$pvalue<= 0.05)
# dim(Control_vs_Test_up6)
# openxlsx::write.xlsx(Control_vs_Test_up6, file = "Control_vs_Test/Sample3_vs_Sample7.Significant_Up_1.xlsx")
# 
# Control_vs_Test_up7=subset(Control_vs_Test_DE7, Control_vs_Test_DE7$log2FoldChange >= 1 & Control_vs_Test_DE7$pvalue<= 0.05)
# dim(Control_vs_Test_up7)
# openxlsx::write.xlsx(Control_vs_Test_up7, file = "Control_vs_Test/Sample4_vs_Sample8.Significant_Up_1.xlsx")
# 
# Control_vs_Test_up8=subset(Control_vs_Test_DE8, Control_vs_Test_DE8$log2FoldChange >= 1 & Control_vs_Test_DE8$pvalue<= 0.05)
# dim(Control_vs_Test_up8)
# openxlsx::write.xlsx(Control_vs_Test_up8, file = "Control_vs_Test/Sample7_vs_Sample8.Significant_Up_1.xlsx")
# 

MSWT3_MSGN2_down1=subset(MSWT3_MSGN2_DE1, MSWT3_MSGN2_DE1$log2FoldChange <= -1 & MSWT3_MSGN2_DE1$pvalue<= 0.05)
dim(MSWT3_MSGN2_down1)
openxlsx::write.xlsx(MSWT3_MSGN2_down1, file = "MSWT3_MSGN2/MSWT3_vs_MSGN2.Significant_Down_.xlsx")

MSLT2_MS202_down2=subset(MSLT2_MS202_DE2, MSLT2_MS202_DE2$log2FoldChange <= -1 & MSLT2_MS202_DE2$pvalue<= 0.05)
dim(MSLT2_MS202_down2)
openxlsx::write.xlsx(MSLT2_MS202_down2, file = "MSLT2_MS202/MSLT2_vs_MS202.Significant_Down_.xlsx")

MSLT2_MS204_down3=subset(MSLT2_MS204_DE3, MSLT2_MS204_DE3$log2FoldChange <= -1 & MSLT2_MS204_DE3$pvalue<= 0.05)
dim(MSLT2_MS204_down3)
openxlsx::write.xlsx(MSLT2_MS204_down3, file = "MSLT2_MS204/MSLT2_vs_MS204.Significant_Down_.xlsx")

MS202_MS204_down4=subset(MS202_MS204_DE4, MS202_MS204_DE4$log2FoldChange <= -1 & MS202_MS204_DE4$pvalue<= 0.05)
dim(MS202_MS204_down4)
openxlsx::write.xlsx(MS202_MS204_down4, file = "MS202_MS204/MS202_vs_MS204.Significant_Down_.xlsx")

# Control_vs_Test_down5=subset(Control_vs_Test_DE5, Control_vs_Test_DE5$log2FoldChange <= -1 & Control_vs_Test_DE5$pvalue<= 0.05)
# dim(Control_vs_Test_down5)
# openxlsx::write.xlsx(Control_vs_Test_down5, file = "Control_vs_Test/Sample3_vs_Sample4.Significant_Down_.xlsx")
# 
# Control_vs_Test_down6=subset(Control_vs_Test_DE6, Control_vs_Test_DE6$log2FoldChange <= -1 & Control_vs_Test_DE6$pvalue<= 0.05)
# dim(Control_vs_Test_down6)
# openxlsx::write.xlsx(Control_vs_Test_down6, file = "Control_vs_Test/Sample3_vs_Sample7.Significant_Down_.xlsx")
# 
# Control_vs_Test_down7=subset(Control_vs_Test_DE7, Control_vs_Test_DE7$log2FoldChange <= -1 & Control_vs_Test_DE7$pvalue<= 0.05)
# dim(Control_vs_Test_down7)
# openxlsx::write.xlsx(Control_vs_Test_down7, file = "Control_vs_Test/Sample4_vs_Sample8.Significant_Down_.xlsx")
# 
# Control_vs_Test_down8=subset(Control_vs_Test_DE8, Control_vs_Test_DE8$log2FoldChange <= -1 & Control_vs_Test_DE8$pvalue<= 0.05)
# dim(Control_vs_Test_down8)
# openxlsx::write.xlsx(Control_vs_Test_down8, file = "Control_vs_Test/Sample7_vs_Sample8.Significant_Down_.xlsx")

# Make a basic volcano plot
#par(bg = 'lightgrey')
png("MSWT3_MSGN2/MSWT3_vs_MSGN2.Volcano_plot.png", type="cairo", res = 150, units="in", width=6, height=7)
plot(MSWT3_MSGN2_DE1$log2FoldChange,-log10(MSWT3_MSGN2_DE1$pvalue), main="Differential Expression", pch=20, col="darkgrey", xlim=c(-8,8), xlab="log2 Fold change", ylab="-log10(p value)")
abline(v=c(-1,1), h=1.3, col="hotpink", lty="dotdash")
grid (NULL,NULL, lty = 3, col = "cornsilk2") 
#color the points as per p-value
with(subset(MSWT3_MSGN2_DE1, pvalue<=0.05 & abs(log2FoldChange)>=1), points(log2FoldChange, -log10(pvalue), pch=20, col="turquoise3"))
with(subset(MSWT3_MSGN2_DE1, padj<=0.05 & abs(log2FoldChange)>=1), points(log2FoldChange, -log10(pvalue), pch=20, col="tomato2"))
#legend(-7.9,0.2, c("adj.p<=0.05","p value<=0.05"), xjust=0, box.lty = "blank", text.col = c("tomato2", "turquoise3"), cex = 0.5)
legend(-7.9,0.3, c("p value<=0.05"), xjust=0, yjust = 0.5, box.lty = "blank", text.col = c("turquoise3","tomato2"), cex = 0.5)

dev.off()

png("MSLT2_MS202/MSLT2_vs_MS202.Volcano_plot.png", type="cairo", res = 150, units="in", width=6, height=7)
plot(MSLT2_MS202_DE2$log2FoldChange,-log10(MSLT2_MS202_DE2$pvalue), main="Differential Expression", pch=20, col="darkgrey", xlim=c(-8,8), xlab="log2 Fold change", ylab="-log10(p value)")
abline(v=c(-1,1), h=1.3, col="hotpink", lty="dotdash")
grid (NULL,NULL, lty = 3, col = "cornsilk2") 
#color the points as per p-value
with(subset(MSLT2_MS202_DE2, pvalue<=0.05 & abs(log2FoldChange)>=1), points(log2FoldChange, -log10(pvalue), pch=20, col="turquoise3"))
with(subset(MSLT2_MS202_DE2, padj<=0.05 & abs(log2FoldChange)>=1), points(log2FoldChange, -log10(pvalue), pch=20, col="tomato2"))
#legend(-7.9,0.2, c("adj.p<=0.05","p value<=0.05"), xjust=0, box.lty = "blank", text.col = c("tomato2", "turquoise3"), cex = 0.5)
legend(-7.9,0.3, c("p value<=0.05"), xjust=0, yjust = 0.5, box.lty = "blank", text.col = c("turquoise3","tomato2"), cex = 0.5)

dev.off()

png("MSLT2_MS204/MSLT2_vs_MS204.Volcano_plot.png", type="cairo", res = 150, units="in", width=6, height=7)
plot(MSLT2_MS204_DE3$log2FoldChange,-log10(MSLT2_MS204_DE3$pvalue), main="Differential Expression", pch=20, col="darkgrey", xlim=c(-8,8), xlab="log2 Fold change", ylab="-log10(p value)")
abline(v=c(-1,1), h=1.3, col="hotpink", lty="dotdash")
grid (NULL,NULL, lty = 3, col = "cornsilk2") 
#color the points as per p-value
with(subset(MSLT2_MS204_DE3, pvalue<=0.05 & abs(log2FoldChange)>=1), points(log2FoldChange, -log10(pvalue), pch=20, col="turquoise3"))
with(subset(MSLT2_MS204_DE3, padj<=0.05 & abs(log2FoldChange)>=1), points(log2FoldChange, -log10(pvalue), pch=20, col="tomato2"))
#legend(-7.9,0.2, c("adj.p<=0.05","p value<=0.05"), xjust=0, box.lty = "blank", text.col = c("tomato2", "turquoise3"), cex = 0.5)
legend(-7.9,0.3, c("p value<=0.05"), xjust=0, yjust = 0.5, box.lty = "blank", text.col = c("turquoise3","tomato2"), cex = 0.5)

dev.off()

png("MS202_MS204/MS202_vs_MS204.Volcano_plot.png", type="cairo", res = 150, units="in", width=6, height=7)
plot(MS202_MS204_DE4$log2FoldChange,-log10(MS202_MS204_DE4$pvalue), main="Differential Expression", pch=20, col="darkgrey", xlim=c(-8,8), xlab="log2 Fold change", ylab="-log10(p value)")
abline(v=c(-1,1), h=1.3, col="hotpink", lty="dotdash")
grid (NULL,NULL, lty = 3, col = "cornsilk2") 
#color the points as per p-value
with(subset(MS202_MS204_DE4, pvalue<=0.05 & abs(log2FoldChange)>=1), points(log2FoldChange, -log10(pvalue), pch=20, col="turquoise3"))
with(subset(MS202_MS204_DE4, padj<=0.05 & abs(log2FoldChange)>=1), points(log2FoldChange, -log10(pvalue), pch=20, col="tomato2"))
#legend(-7.9,0.2, c("adj.p<=0.05","p value<=0.05"), xjust=0, box.lty = "blank", text.col = c("tomato2", "turquoise3"), cex = 0.5)
legend(-7.9,0.3, c("p value<=0.05"), xjust=0, yjust = 0.5, box.lty = "blank", text.col = c("turquoise3","tomato2"), cex = 0.5)

dev.off()

#require("ggplot2") #Best plots
require("ggrepel") #Avoid overlapping labels
require("plyr")
as.matrix(names(MSWT3_MSGN2_DE1))
#mutateddf <- mutate(Control_vs_Test_Control_DE,Legend=ifelse(padj<=0.05 & abs(log2FoldChange)>=2, "log2FoldChange>=2",ifelse(padj>=0.05 & abs(log2FoldChange)>=2, "log2FoldChange>=2" , "Not Significant")))
mutateddf <- mutate(MSWT3_MSGN2_DE1,Legend=ifelse(padj<=0.05 & abs(log2FoldChange)>=1 & pvalue<=0.05 , "adj.p<=0.05",ifelse(padj>=0.05 & pvalue<=0.05 & abs(log2FoldChange)>=1, "p value<=0.05" , "Not Significant")))

input <- cbind(gene=rownames(mutateddf ), mutateddf ) #convert the rownames to a column
xmini=min(MSWT3_MSGN2_DE1$log2FoldChange)-4
xmaxi=max(MSWT3_MSGN2_DE1$log2FoldChange)
volc = ggplot(na.omit(input), aes(log2FoldChange, -log10(pvalue))) + #volcanoplot with log2Foldchange versus pvalue
  geom_point(aes(col=Legend)) + #add points colored by significance
  scale_color_manual(values=c("grey","turquoise3","tomato2")) +
  #scale_color_manual(values=c("grey","turquoise3")) +
  scale_x_continuous(limits = c(-10,10)) +
  # scale_y_continuous(limits = c(0,20)) + 
  geom_hline(yintercept = 1.3,size=0.3, color="hotpink", linetype="dashed") +
  geom_vline(xintercept = c(-1,1),size=0.3, color="hotpink", linetype="dashed") +
  ggtitle("Differential Expression") +
  theme(plot.title = element_text(hjust = 0.5))
volc=volc+geom_text_repel(data=head(input.neg[order(input.neg$pvalue,decreasing = F),],20), aes(label=gene_name),size = 3, color = 'black', box.padding = 0.1, point.padding = 0.1, segment.color = 'grey50', max.overlaps = 20) #adding text for the top 10  genes
volc=volc+geom_text_repel(data=head(input.pos[order(input.pos$pvalue,decreasing = F),],20), aes(label=gene_name),size = 3, color = 'black', box.padding = 0.1, point.padding = 0.1, segment.color = 'grey50',max.overlaps = 20) #adding text for the top 10 genes
input.pos=subset(input, log2FoldChange>=1 & pvalue <= 0.05)
#input.neg=subset(input, log2FoldChange>=-1 & pvalue <= 0.05)
input.neg=subset(input, log2FoldChange<=-1 & pvalue <= 0.05)
#volc=volc+geom_text_repel(data=head(input.pos[order(input.pos$pvalue,decreasing = F),],25), aes(label=Id),size = 3, color = 'black', box.padding = 0.1, point.padding = 0.1, segment.color = 'grey50') #adding text for the top 10  genes
#volc=volc+geom_text_repel(data=head(input.neg[order(input.neg$pvalue,decreasing = F),],25), aes(label=Id),size = 3, color = 'black', box.padding = 0.1, point.padding = 0.1, segment.color = 'grey50') #adding text for the top 10 genes
ggsave(volc, filename = "MSWT3_MSGN2/MSWT3_vs_MSGN2.Volcano_plot.genename_1.png",height = 10,width = 10,device = "png",units = "in",dpi = 300)
dev.off()
#ggsave(volc, filename = "Control_vs_ RMTP10_1/Control_vs_ RMTP10.Volcano_plot.genename.png")
rm(list="xmini", "xmaxi" , "input.pos", "input.neg","mutateddf","input")

as.matrix(names(MSLT2_MS202_DE2))
#mutateddf <- mutate(Control_vs_Test_Control_DE,Legend=ifelse(padj<=0.05 & abs(log2FoldChange)>=2, "log2FoldChange>=2",ifelse(padj>=0.05 & abs(log2FoldChange)>=2, "log2FoldChange>=2" , "Not Significant")))
mutateddf <- mutate(MSLT2_MS202_DE2,Legend=ifelse(padj<=0.05 & abs(log2FoldChange)>=1 & pvalue<=0.05 , "adj.p<=0.05",ifelse(padj>=0.05 & pvalue<=0.05 & abs(log2FoldChange)>=1, "p value<=0.05" , "Not Significant")))

input <- cbind(gene=rownames(mutateddf ), mutateddf ) #convert the rownames to a column
xmini=min(MSLT2_MS202_DE2$log2FoldChange)-4
xmaxi=max(MSLT2_MS202_DE2$log2FoldChange)
volc = ggplot(na.omit(input), aes(log2FoldChange, -log10(pvalue))) + #volcanoplot with log2Foldchange versus pvalue
  geom_point(aes(col=Legend)) + #add points colored by significance
  scale_color_manual(values=c("grey","turquoise3","tomato2")) +
  #scale_color_manual(values=c("grey","turquoise3")) +
  scale_x_continuous(limits = c(-10,10)) +
  # scale_y_continuous(limits = c(0,20)) + 
  geom_hline(yintercept = 1.3,size=0.3, color="hotpink", linetype="dashed") +
  geom_vline(xintercept = c(-1,1),size=0.3, color="hotpink", linetype="dashed") +
  ggtitle("Differential Expression") +
  theme(plot.title = element_text(hjust = 0.5))
volc=volc+geom_text_repel(data=head(input.neg[order(input.neg$pvalue,decreasing = F),],20), aes(label=gene_name),size = 3, color = 'black', box.padding = 0.1, point.padding = 0.1, segment.color = 'grey50') #adding text for the top 10  genes
volc=volc+geom_text_repel(data=head(input.pos[order(input.pos$pvalue,decreasing = F),],20), aes(label=gene_name),size = 3, color = 'black', box.padding = 0.1, point.padding = 0.1, segment.color = 'grey50') #adding text for the top 10 genes
input.pos=subset(input, log2FoldChange>=1 & pvalue <= 0.05)
#input.neg=subset(input, log2FoldChange>=-1 & pvalue <= 0.05)
input.neg=subset(input, log2FoldChange<=-1 & pvalue <= 0.05)
#volc=volc+geom_text_repel(data=head(input.pos[order(input.pos$pvalue,decreasing = F),],25), aes(label=Id),size = 3, color = 'black', box.padding = 0.1, point.padding = 0.1, segment.color = 'grey50') #adding text for the top 10  genes
#volc=volc+geom_text_repel(data=head(input.neg[order(input.neg$pvalue,decreasing = F),],25), aes(label=Id),size = 3, color = 'black', box.padding = 0.1, point.padding = 0.1, segment.color = 'grey50') #adding text for the top 10 genes
ggsave(volc, filename = "MSLT2_MS202/MSLT2_vs_MS202.Volcano_plot.genename.png",height = 10,width = 10,device = "png",units = "in",dpi = 300)
dev.off()
#ggsave(volc, filename = "Control_vs_ RMTP10_1/Control_vs_ RMTP10.Volcano_plot.genename.png")
rm(list="xmini", "xmaxi" , "input.pos", "input.neg","mutateddf","input")

as.matrix(names(MSLT2_MS204_DE3))
#mutateddf <- mutate(Control_vs_Test_Control_DE,Legend=ifelse(padj<=0.05 & abs(log2FoldChange)>=2, "log2FoldChange>=2",ifelse(padj>=0.05 & abs(log2FoldChange)>=2, "log2FoldChange>=2" , "Not Significant")))
mutateddf <- mutate(MSLT2_MS204_DE3,Legend=ifelse(ifelse(padj>=0.05 & pvalue<=0.05 & abs(log2FoldChange)>=1, "p value<=0.05" , "Not Significant")))

input <- cbind(gene=rownames(mutateddf ), mutateddf ) #convert the rownames to a column
xmini=min(MSLT2_MS204_DE3$log2FoldChange)-4
xmaxi=max(MSLT2_MS204_DE3$log2FoldChange)
volc = ggplot(na.omit(input), aes(log2FoldChange, -log10(pvalue))) + #volcanoplot with log2Foldchange versus pvalue
  geom_point(aes(col=Legend)) + #add points colored by significance
  scale_color_manual(values=c("grey","turquoise3","tomato2")) +
  #scale_color_manual(values=c("grey","turquoise3")) +
  scale_x_continuous(limits = c(-10,10)) +
  # scale_y_continuous(limits = c(0,20)) + 
  geom_hline(yintercept = 1.3,size=0.3, color="hotpink", linetype="dashed") +
  geom_vline(xintercept = c(-1,1),size=0.3, color="hotpink", linetype="dashed") +
  ggtitle("Differential Expression") +
  theme(plot.title = element_text(hjust = 0.5))
volc=volc+geom_text_repel(data=head(input.neg[order(input.neg$pvalue,decreasing = F),],20), aes(label=gene_name),size = 3, color = 'black', box.padding = 0.1, point.padding = 0.1, segment.color = 'grey50', max.overlaps = 20) #adding text for the top 10  genes
volc=volc+geom_text_repel(data=head(input.pos[order(input.pos$pvalue,decreasing = F),],20), aes(label=gene_name),size = 3, color = 'black', box.padding = 0.1, point.padding = 0.1, segment.color = 'grey50',max.overlaps = 20) #adding text for the top 10 genes
input.pos=subset(input, log2FoldChange>=1 & pvalue <= 0.05)
#input.neg=subset(input, log2FoldChange>=-1 & pvalue <= 0.05)
input.neg=subset(input, log2FoldChange<=-1 & pvalue <= 0.05)
#volc=volc+geom_text_repel(data=head(input.pos[order(input.pos$pvalue,decreasing = F),],25), aes(label=Id),size = 3, color = 'black', box.padding = 0.1, point.padding = 0.1, segment.color = 'grey50') #adding text for the top 10  genes
#volc=volc+geom_text_repel(data=head(input.neg[order(input.neg$pvalue,decreasing = F),],25), aes(label=Id),size = 3, color = 'black', box.padding = 0.1, point.padding = 0.1, segment.color = 'grey50') #adding text for the top 10 genes
ggsave(volc, filename = "MSLT2_MS204/MSLT2_vs_MS204.Volcano_plot.genename_1.png",height = 10,width = 10,device = "png",units = "in",dpi = 300)
dev.off()

as.matrix(names(MS202_MS204_DE4))
#mutateddf <- mutate(Control_vs_Test_Control_DE,Legend=ifelse(padj<=0.05 & abs(log2FoldChange)>=2, "log2FoldChange>=2",ifelse(padj>=0.05 & abs(log2FoldChange)>=2, "log2FoldChange>=2" , "Not Significant")))
mutateddf <- mutate(MS202_MS204_DE4,Legend=ifelse(padj<=0.05 & abs(log2FoldChange)>=1 & pvalue<=0.05 , "adj.p<=0.05",ifelse(padj>=0.05 & pvalue<=0.05 & abs(log2FoldChange)>=1, "p value<=0.05" , "Not Significant")))

input <- cbind(gene=rownames(mutateddf ), mutateddf ) #convert the rownames to a column
xmini=min(MS202_MS204_DE4$log2FoldChange)-4
xmaxi=max(MS202_MS204_DE4$log2FoldChange)
volc = ggplot(na.omit(input), aes(log2FoldChange, -log10(pvalue))) + #volcanoplot with log2Foldchange versus pvalue
  geom_point(aes(col=Legend)) + #add points colored by significance
  scale_color_manual(values=c("grey","turquoise3","tomato2")) +
  #scale_color_manual(values=c("grey","turquoise3")) +
  scale_x_continuous(limits = c(-10,10)) +
  # scale_y_continuous(limits = c(0,20)) + 
  geom_hline(yintercept = 1.3,size=0.3, color="hotpink", linetype="dashed") +
  geom_vline(xintercept = c(-1,1),size=0.3, color="hotpink", linetype="dashed") +
  ggtitle("Differential Expression") +
  theme(plot.title = element_text(hjust = 0.5))
volc=volc+geom_text_repel(data=head(input.neg[order(input.neg$pvalue,decreasing = F),],20), aes(label=gene_name),size = 3, color = 'black', box.padding = 0.1, point.padding = 0.1, segment.color = 'grey50', max.overlaps = 20) #adding text for the top 10  genes
volc=volc+geom_text_repel(data=head(input.pos[order(input.pos$pvalue,decreasing = F),],20), aes(label=gene_name),size = 3, color = 'black', box.padding = 0.1, point.padding = 0.1, segment.color = 'grey50',max.overlaps = 20) #adding text for the top 10 genes
input.pos=subset(input, log2FoldChange>=1 & pvalue <= 0.05)
#input.neg=subset(input, log2FoldChange>=-1 & pvalue <= 0.05)
input.neg=subset(input, log2FoldChange<=-1 & pvalue <= 0.05)
#volc=volc+geom_text_repel(data=head(input.pos[order(input.pos$pvalue,decreasing = F),],25), aes(label=Id),size = 3, color = 'black', box.padding = 0.1, point.padding = 0.1, segment.color = 'grey50') #adding text for the top 10  genes
#volc=volc+geom_text_repel(data=head(input.neg[order(input.neg$pvalue,decreasing = F),],25), aes(label=Id),size = 3, color = 'black', box.padding = 0.1, point.padding = 0.1, segment.color = 'grey50') #adding text for the top 10 genes
ggsave(volc, filename = "MS202_MS204/MS202_vs_MS204.Volcano_plot.genename_1.png",height = 10,width = 10,device = "png",units = "in",dpi = 300)
dev.off()
#ggsave(volc, filename = "Control_vs_ RMTP10_1/Control_vs_ RMTP10.Volcano_plot.genename.png")
rm(list="xmini", "xmaxi" , "input.pos", "input.neg","mutateddf","input")





###heatmap####
MSWT3_MSGN2_DE.sig.norm1= normcounts1[rownames(normcounts1) %in% MSWT3_MSGN2_DE.sig1$Id,]
dim(MSWT3_MSGN2_DE.sig.norm1)
MSWT3_MSGN2_DE.sig.norm1=data.frame(MSWT3_MSGN2_DE.sig.norm1)
as.matrix(names(MSWT3_MSGN2_DE.sig.norm1))
#colnames(Test_vs_CAPE_DE.sig.norm)=c("Test_Control1" , "Test_Control2" ,"Control", "Control1" , "Control2")
#Control_vs_Test_Control_DE.sig.norm=Control_vs_Test_Control_DE.sig.norm[,c(1,3)]

require(pheatmap)
png("MSWT3_MSGN2/MSWT3_vs_MSGN2.Heatmap_R.png", type="cairo", res = 200, units="in", width=8, height=25)
#anno=data.frame(factor(c(rep("Test_Control1",1),rep("Test_Control2",1),rep("Control",1),rep("Control1",1),rep("Control2",1))))
anno=data.frame(factor(c(rep("MSWT3",1),rep("MSGN2",1))))
rownames(anno)=colnames(data.frame(MSWT3_MSGN2_DE.sig.norm1))
colnames(anno)="Group"
newCols <- colorRampPalette(colorspace::heat_hcl(length(unique(anno$Group))))
mycolors <- newCols(length(unique(anno$Group)))
names(mycolors) <- unique(anno$Group)
mycolors <- list(Group = mycolors)
pheatmap(data.frame(MSWT3_MSGN2_DE.sig.norm1),color = colorRampPalette(c("green4", "yellow", "firebrick2"))(50),show_rownames = F,show_colnames = T,cluster_cols = T,scale = "row",annotation= anno,annotation_colors = mycolors,border_color = F,labels_col = colnames(MSWT3_MSGN2_DE.sig.norm1),fontsize_col = 8,fontsize = 6,angle_col = 0)
dev.off()

rm(MSWT3_MSGN2_DE.sig.norm1)

MSLT2_MS202_DE.sig.norm2= normcounts2[rownames(normcounts2) %in% MSLT2_MS202_DE.sig2$Id,]
dim(MSLT2_MS202_DE.sig.norm2)
MSLT2_MS202_DE.sig.norm2=data.frame(MSLT2_MS202_DE.sig.norm2)
as.matrix(names(MSLT2_MS202_DE.sig.norm2))
#colnames(Test_vs_CAPE_DE.sig.norm)=c("Test_Control1" , "Test_Control2" ,"Control", "Control1" , "Control2")
#Control_vs_Test_Control_DE.sig.norm=Control_vs_Test_Control_DE.sig.norm[,c(1,3)]

require(pheatmap)
png("MSLT2_MS202/MSLT2_vs_MS202.Heatmap_R.png", type="cairo", res = 200, units="in", width=8, height=25)
#anno=data.frame(factor(c(rep("Test_Control1",1),rep("Test_Control2",1),rep("Control",1),rep("Control1",1),rep("Control2",1))))
anno=data.frame(factor(c(rep("MSLT2_1",1),rep("MS202_3",1))))
rownames(anno)=colnames(data.frame(MSLT2_MS202_DE.sig.norm2))
colnames(anno)="Group"
newCols <- colorRampPalette(colorspace::heat_hcl(length(unique(anno$Group))))
mycolors <- newCols(length(unique(anno$Group)))
names(mycolors) <- unique(anno$Group)
mycolors <- list(Group = mycolors)
pheatmap(data.frame(MSLT2_MS202_DE.sig.norm2),color = colorRampPalette(c("green4", "yellow", "firebrick2"))(50),show_rownames = F,show_colnames = T,cluster_cols = T,scale = "row",annotation= anno,annotation_colors = mycolors,border_color = F,labels_col = colnames(MSLT2_MS202_DE.sig.norm2),fontsize_col = 8,fontsize = 6,angle_col = 0)
dev.off()

rm(MSLT2_MS202_DE.sig.norm2)

MSLT2_MS202_DE.sig.norm2= normcounts2[rownames(normcounts2) %in% MSLT2_MS202_DE.sig2$Id,]
dim(MSLT2_MS202_DE.sig.norm2)
MSLT2_MS202_DE.sig.norm2=data.frame(MSLT2_MS202_DE.sig.norm2)
as.matrix(names(MSLT2_MS202_DE.sig.norm2))
#colnames(Test_vs_CAPE_DE.sig.norm)=c("Test_Control1" , "Test_Control2" ,"Control", "Control1" , "Control2")
#Control_vs_Test_Control_DE.sig.norm=Control_vs_Test_Control_DE.sig.norm[,c(1,3)]

require(pheatmap)
png("MSLT2_MS202/MSLT2_vs_MS202.Heatmap_R.png", type="cairo", res = 200, units="in", width=8, height=25)
#anno=data.frame(factor(c(rep("Test_Control1",1),rep("Test_Control2",1),rep("Control",1),rep("Control1",1),rep("Control2",1))))
anno=data.frame(factor(c(rep("MSLT2_1",1),rep("MS202_3",1))))
rownames(anno)=colnames(data.frame(MSLT2_MS202_DE.sig.norm2))
colnames(anno)="Group"
newCols <- colorRampPalette(colorspace::heat_hcl(length(unique(anno$Group))))
mycolors <- newCols(length(unique(anno$Group)))
names(mycolors) <- unique(anno$Group)
mycolors <- list(Group = mycolors)
pheatmap(data.frame(MSLT2_MS202_DE.sig.norm2),color = colorRampPalette(c("green4", "yellow", "firebrick2"))(50),show_rownames = F,show_colnames = T,cluster_cols = T,scale = "row",annotation= anno,annotation_colors = mycolors,border_color = F,labels_col = colnames(MSLT2_MS202_DE.sig.norm2),fontsize_col = 8,fontsize = 6,angle_col = 0)
dev.off()

rm(MSLT2_MS202_DE.sig.norm2)

MSLT2_MS204_DE.sig.norm3= normcounts3[rownames(normcounts3) %in% MSLT2_MS204_DE.sig3$Id,]
dim(MSLT2_MS204_DE.sig.norm3)
MSLT2_MS204_DE.sig.norm3=data.frame(MSLT2_MS204_DE.sig.norm3)
as.matrix(names(MSLT2_MS204_DE.sig.norm3))
#colnames(Test_vs_CAPE_DE.sig.norm)=c("Test_Control1" , "Test_Control2" ,"Control", "Control1" , "Control2")
#Control_vs_Test_Control_DE.sig.norm=Control_vs_Test_Control_DE.sig.norm[,c(1,3)]

require(pheatmap)
png("MSLT2_MS204/MSLT2_vs_MS204.Heatmap_R.png", type="cairo", res = 200, units="in", width=8, height=25)
#anno=data.frame(factor(c(rep("Test_Control1",1),rep("Test_Control2",1),rep("Control",1),rep("Control1",1),rep("Control2",1))))
anno=data.frame(factor(c(rep("MSLT2_1",1),rep("MS204_3",1))))
rownames(anno)=colnames(data.frame(MSLT2_MS204_DE.sig.norm3))
colnames(anno)="Group"
newCols <- colorRampPalette(colorspace::heat_hcl(length(unique(anno$Group))))
mycolors <- newCols(length(unique(anno$Group)))
names(mycolors) <- unique(anno$Group)
mycolors <- list(Group = mycolors)
pheatmap(data.frame(MSLT2_MS204_DE.sig.norm3),color = colorRampPalette(c("green4", "yellow", "firebrick2"))(50),show_rownames = F,show_colnames = T,cluster_cols = T,scale = "row",annotation= anno,annotation_colors = mycolors,border_color = F,labels_col = colnames(MSLT2_MS204_DE.sig.norm3),fontsize_col = 8,fontsize = 6,angle_col = 0)
dev.off()

rm(MSLT2_MS204_DE.sig.norm3)

###############################################################################################################
###############################################################################################################
#load stringer package
require(tidyr)
require(readr)
require(stringr)
require(devtools)
require(ggfortify)

#MAKE LIST OF de FILES
diff.merged=mget(ls(pattern = "_DE"))
diff.merged$MSWT3_MSGN2_DE.sig1=NULL

#make lists of sig, up, doen
diff.list.sig=mget(ls(pattern = "_DE.sig"))
diff.list.Up=mget(ls(pattern = "_up"))

diff.list.Down=mget(ls(pattern = "_down"))


#######################GENE ONTOLOGY Analysis##############################
require(clusterProfiler)
require(GOplot)
require(grid)
#load organism database

require(org.Hs.eg.db)
library("org.Hs.eg.db")
keytypes(org.Hs.eg.db)
require(tidyr)

universelist <- unique(separate(data = diff.merged[[13]][1], col = Id, into = c("Id", "tmp"), sep = "\\."))
universelist <- universelist[1]
head(universelist)
universelist=as.character(unlist(unique(mapIds(org.Hs.eg.db, keys = as.character(unlist(universelist[1])),column = "ENTREZID", keytype = "ENSEMBL",multiVals="first"))))
head(universelist)

genelist <- unique(separate(data = diff.list.sig[[5]][1], col = Id, into = c("Id", "tmp"), sep = "\\."))
genelist <- genelist[1]
head(genelist)
genelist=as.character(unlist(unique(mapIds(org.Hs.eg.db, keys = as.character(unlist(genelist[1])),column = "ENTREZID", keytype = "ENSEMBL",multiVals="first"))))
head(genelist)

enrichGO(gene=genelist, universe = universelist, OrgDb = org.Hs.eg.db, pvalueCutoff = 0.05,
         ont = "BP", pAdjustMethod = "BH", minGSSize = 1, keyType = "ENTREZID", qvalueCutoff = 0.05,
         readable=T) -> tmp

if(class(tmp)=="NULL"){next}else{
  assign(paste(names(diff.list.sig)[5], "BP", sep="_"),
         slot(simplify(tmp, cutoff=0.7, by="p.adjust", select_fun=min),"result"))}

rm(tmp)

enrichGO(gene=genelist, universe = universelist, OrgDb = org.Hs.eg.db, pvalueCutoff = 0.05,
         ont = "CC", pAdjustMethod = "BH", minGSSize = 1, keyType = "ENTREZID", qvalueCutoff = 0.05,
         readable=T) -> tmp

if(class(tmp)=="NULL"){next}else{
  assign(paste(names(diff.list.sig)[5], "CC", sep="_"),
         slot(simplify(tmp, cutoff=0.7, by="p.adjust", select_fun=min),"result"))
}

rm(tmp)

enrichGO(gene=genelist, universe = universelist, OrgDb = org.Hs.eg.db, pvalueCutoff = 0.05,
         ont = "MF", pAdjustMethod = "BH", minGSSize = 1, keyType = "ENTREZID", qvalueCutoff = 0.05,
         readable=T) -> tmp

if(class(tmp)=="NULL"){next}else{
  assign(paste(names(diff.list.sig)[5], "MF", sep="_"),
         slot(simplify(tmp, cutoff=0.7, by="p.adjust", select_fun=min),"result"))
}
rm(tmp)

allGO=mget(ls(pattern = "BP|MF|CC"))
rm(list = ls(pattern = "BP|MF|CC"))

for(i in 1:length(allGO))
{
  if(grepl("BP",names(allGO[i]))==T)
  {
    allGO[[i]]$Category=c(rep("BP", dim(allGO[[i]])[1]))
  }
  else if(grepl("MF",names(allGO[i]))==T)
  {
    allGO[[i]]$Category=c(rep("MF", dim(allGO[[i]])[1]))
  }
  else if(grepl("CC",names(allGO[i]))==T)
  {
    allGO[[i]]$Category=c(rep("CC", dim(allGO[[i]])[1]))
  }
}
rm(i)
compariosions=unique(do.call(rbind,str_split(names(allGO),"_DE."))[,2])
for(i in 1:length(compariosions))
{
  assign(paste(compariosions[i],"GO_export",sep="_"),
         do.call(rbind,allGO[grep(compariosions[i], names(allGO))])[,c(10,1:6,8,9)])
}
rm(i)

GO.export=mget(ls(pattern = "_GO_export"))
rm(list=ls(pattern = "_GO_export"))

for(i in 1:length(GO.export))
{
  openxlsx::write.xlsx(x = GO.export[[i]],file = paste0(str_split(names(GO.export)[i], pattern = "__")[[1]][1],
                                                        ".Gene_Ontology_Enrichment5.xlsx"))
}
rm(i)

compariosions=unique(do.call(rbind,str_split(names(allGO),"_DE."))[,2])
for(i in 1:length(compariosions)){
  assign(paste(compariosions[i],"combined_GO",sep="_"),
         do.call(rbind,allGO[grep(compariosions[i], names(allGO))])[,c(10,1,2,6,8)])
}
rm(i)
combi.GO=mget(ls(pattern = "_combined_GO"))
rm(list=ls(pattern = "_combined_GO"))

for(i in 1:length(combi.GO)){
  names(combi.GO[[i]])=paste(c("Category", "ID", "term", "adj_pval","Genes"))
  combi.GO[[i]]$Genes=gsub("/", ",", combi.GO[[i]]$Genes)
}
rm(i)
for(i in 1:length(combi.GO)){
  names(combi.GO[[i]])=paste(c("Category", "ID", "term", "adj_pval","genes"))
}
genelist.GO=diff.list.sig[5]

for(i in 1:length(genelist.GO)){
  genelist.GO[[i]]=genelist.GO[[i]][,c(11,6,9,10,1)]
  names(genelist.GO[[i]])=paste(c("ID", "logFC", "pvalue","padj","ENSEMBL"))
}
rm(i)

for(i in 1:length(compariosions))
{
  assign(paste(compariosions[2], "circ",sep="_"),
         circle_dat(combi.GO[[grep(compariosions[2], names(combi.GO))]], genelist.GO$Control_vs_Test_DE.sig5))
}
rm(i)
GO.circ=mget(ls(pattern = "_circ"))

#merge(authors, books, by.x="surname", by.y="name")
#l1= merge.data.frame(sig_BP_circ, sig_CC_circ, all=T)
#l2= merge.data.frame(l1, sig_MF_circ, all=T)
l1= merge.data.frame(sig5_BP_circ, sig5_CC_circ, all=T)
l2=merge(l1, sig5_MF_circ, all=T)
#head(l2)
#names(l2)
data.circ=l1[!duplicated(l1$term),]
data.circ2=subset(data.circ, data.circ$count>1)
xmini=min(data.circ2$zscore)-2
xmaxi=max(data.circ2$zscore)+2
sub2.data.circ=subset(data.circ2, data.circ2$adj_pval<=0.1)
rang_data.circ <- c(1, 10)
cols_data.circ <- c("chartreuse4", "brown2", "cornflowerblue")
bubplot=ggplot(data.circ2, aes(zscore, -log10(adj_pval), (alpha=0.1), fill=category, size = count)) +
  labs(title = "GO bubble plot", x = "z-score", y = "-log10(adj p-value)") +
  theme(plot.title = element_text(hjust = 0.5))+
  geom_hline(yintercept = -log10(0.1), col = "orange") + scale_size(range = c(1,10)) +
  geom_point(shape = 21, colour="black", stroke=0.1)+
  geom_text(data = sub2.data.circ, aes(zscore, -log10(adj_pval), label=term), size = 3, check_overlap = T)+
  theme(legend.position = "bottom", legend.direction = "horizontal")+
  theme(legend.box = "horizontal",legend.key.size = unit(1, "cm"))+
  labs(x = "z-score", y = "-log10(adj p-value)")
ggsave(filename = "Comparison5_GO_ALL.Bubble_plot5.png",
       plot = bubplot,device = "png", dpi = 300, height = 15, width = 30, units = "in")

rm(allGO, bubplot,combi.GO,data.circ,data.circ2,genelist.GO,sig5_BP_circ,sig5_MF_circ,sig5_CC_circ, GO.circ,GO.export,l1,l2,sub2.data.circ,universelist1,universelist2,compariosions,compariosions2, pasAnno,pasCts,rang_data.circ,xmaxi,xmini,cols_data.circ)



#KEGG pathway analysis
#require(BiocInstaller)
#biocLite("pathview")

BiocManager::install("pathview")
library("pathview")
setwd()
#diff.list.sig$Id <- sub("\\.\\d", "", diff.list.sig$Id)

keys1=data.frame(lapply(diff.list.sig[[5]][1], as.character),stringsAsFactors = F)
keys1= unique(separate(data = keys1, col = Id, into = c("Id", "tmp"), sep = "\\."))[1]
entrez.list=AnnotationDbi::select(org.Hs.eg.db, keys1$Id,"ENTREZID","ENSEMBL")

keys2=data.frame(lapply(diff.merged[[13]][1], as.character),stringsAsFactors = F)
keys2= unique(separate(data = keys2, col = Id, into = c("Id", "tmp"), sep = "\\."))[1]
entrez.list.universe=AnnotationDbi::select(org.Hs.eg.db, keys2$Id,"ENTREZID","ENSEMBL")
Kegg_gene.list=as.character(entrez.list$ENTREZID)
kegg_universal.list=as.character(entrez.list.universe$ENTREZID)
kegg.results=enrichKEGG(gene=Kegg_gene.list, organism = "hsa", keyType = "kegg", pvalueCutoff = 0.05,
                        pAdjustMethod = "BH", universe=kegg_universal.list, minGSSize = 2,
                        qvalueCutoff = 0.05, use_internal_data = FALSE)
if(isTRUE(class(kegg.results)=="NULL")) { next } else kegg.out=data.frame(kegg.results@result)
openxlsx::write.xlsx(kegg.out[,c(1:6,8,9)], file=paste0(str_split(names(diff.list.sig[5]),"_DE.")[[1]][2],".KEGG_Pathway_Enrichment5.xlsx"))
tmp.sig.list=diff.list.sig[[5]]
#edit gene_id column(remove "." part from ids (ENSG00000001626.10 as ENSG00000001626))
#tmp.sig.list <- separate(data = tmp.sig.list, col = Id, into = c("Id", "tmp"), sep = "\\.")
tmp.sig.list.with.entrezid=merge(tmp.sig.list, entrez.list, by.x="ENTREZID", by.y="ENTREZID")
genelist=tmp.sig.list.with.entrezid$log2FoldChange
names(genelist)=c(tmp.sig.list.with.entrezid$gene_name)
#top10kegg=kegg.out$ID[1:10]

#for (j in 1:length(kegg.out$ID))
#{
pathview(gene.data=genelist, pathway.id = kegg.out$ID, species = "hsa",
         out.suffix=paste0(str_split(names(diff.list.sig)[5], pattern = "_DE.")[[1]][2]),
         limit = list(gene=max(round(abs(genelist)))), key.align = "y")
#}

# ###################################################
# ###merge gene name and ensembl ids in kegg sheet###
# ###################################################
# ####For Control vs Diabetic
# ###merge by genename
DAVID_LT=readxl::read_xlsx("sig1.KEGG_Pathway_Enrichment1.xlsx")
kegg.split_S1Controldown_vs_S1Mir21pulldown=read.delim(file = "S1-Controldown_vs_S1-Mir21-pulldown/KEGG/S1-Controldown_vs_S1-Mir21-pulldown.KEGG_Pathway_Enrichment.splitted.txt",header = T,sep = "\t",comment.char = "#",stringsAsFactors = F)
as.matrix(names(gff_anno))
gff_anno.kegg=gff_anno[,c(5,1,2)]
kegg.split_S1Controldown_vs_S1Mir21pulldown.geneanno=merge(kegg.split_S1Controldown_vs_S1Mir21pulldown,gff_anno.kegg,by.x = "geneID",by.y = "ENTREZID",all.x = T)
as.matrix(names(kegg.split_S1Controldown_vs_S1Mir21pulldown.geneanno))
kegg.split_S1Controldown_vs_S1Mir21pulldown.geneanno=kegg.split_S1Controldown_vs_S1Mir21pulldown.geneanno[,c(2:8,1,9:10)]
colnames(kegg.split_S1Controldown_vs_S1Mir21pulldown.geneanno) <- sub(pattern = "geneID",replacement = "ENTREZID", colnames(kegg.split_S1Controldown_vs_S1Mir21pulldown.geneanno))
colnames(kegg.split_S1Controldown_vs_S1Mir21pulldown.geneanno) <- sub(pattern = "Id",replacement = "ENSEMBL_ID", colnames(kegg.split_S1Controldown_vs_S1Mir21pulldown.geneanno))
# 
# 
# #collapse miRNAs in all known targets
require(dplyr)
as.matrix(names(kegg.split_S1Controldown_vs_S1Mir21pulldown.geneanno))
kegg.split_S1Controldown_vs_S1Mir21pulldown.geneanno.collapse <- kegg.split_S1Controldown_vs_S1Mir21pulldown.geneanno %>% group_by(ID) %>%  summarize(Description = paste(unique(Description),collapse="; "),Count = unique(Count), GeneRatio = unique(GeneRatio),  BgRatio = unique(BgRatio), pvalue = unique(pvalue), p.adjust = unique(p.adjust), ENTREZID = paste(unique(ENTREZID),collapse="; "), ENSEMBL_ID = paste(unique(ENSEMBL_ID),collapse="; "),  gene_name = paste(unique(gene_name),collapse="; "))# as.matrix(names(kegg.split_S1Controldown_vs_S1Mir21pulldown.geneanno.collapse))
openxlsx::write.xlsx(kegg.split_S1Controldown_vs_S1Mir21pulldown.geneanno.collapse,file = "S1-Controldown_vs_S1-Mir21-pulldown/KEGG/S1-Controldown_vs_S1-Mir21-pulldown.KEGG_Pathway_Enrichment.anno.xlsx")
# 

##############
##UpSetR Plots
setwd(workdir)
dir.create("UpSet_Plot")
dir()
setwd("UpSet_Plot/")

for(i in 1:length(diff.list.Down)){
  write.table(file=paste0(names(diff.list.Down[i]),".txt"), x=diff.list.Down[[i]][1],
              row.names = F, col.names = F, quote = F)
}
for(i in 1:length(diff.list.Up)){
  write.table(file=paste0(names(diff.list.Up[i]),".txt"), x=diff.list.Up[[i]][1],
              row.names = F, col.names = F, quote = F)
}
require(UpSetR)
#use VennPainter in windows to make a matrix file (.mf) [to this give all down and up regulated file of group1]
Venn=read.delim(file = "venn.mf",header = T,sep = "\t",comment.char = "#",stringsAsFactors = F)
Venn$X=as.factor(Venn$X)
Venn=data.frame(Venn)
as.matrix(names(Venn))
colnames(Venn)=c("X" , "Control_vs_Test_down","Control_vs_Test_up")
png("Venn_Diagram.png", type="cairo", res = 200, units="in", width=18, height=16)
upset(Venn,nsets = 14,nintersects = 50,keep.order = F, main.bar.color = "red",
      matrix.color = "blue", order.by = "freq")
dev.off()
setwd("../")


#######Count and make stats for DE,Sig,up and down
setwd(workdir)
#DE counts
DE_counts=NULL
for( i in 1:length(diff.merged))
{
  tmp=cbind(names(diff.merged)[i],dim(diff.merged[[i]])[1])
  DE_counts=rbind(DE_counts,tmp)
}
rm(i)
DE_counts=data.frame(DE_counts)
names(DE_counts)=c("condition" , "DE")
DE_counts <- separate(data = DE_counts, col = condition, into = c("condition", "tmp"), sep = "\\_DE")
DE_counts=DE_counts[,c(1,3)]

#Sig Counts
sig_counts=NULL
for( i in 1:length(diff.list.sig))
{
  tmp=cbind(names(diff.list.sig)[i],dim(diff.list.sig[[i]])[1])
  sig_counts=rbind(sig_counts,tmp)
}
rm(i)
sig_counts=data.frame(sig_counts)
names(sig_counts)=c("condition" , "sig")
sig_counts <- separate(data = sig_counts, col = condition, into = c("condition", "tmp"), sep = "\\_DE")
sig_counts=sig_counts[,c(1,3)]

#Up Counts
up_counts=NULL
for( i in 1:length(diff.list.Up))
{
  tmp=cbind(names(diff.list.Up)[i],dim(diff.list.Up[[i]])[1])
  up_counts=rbind(up_counts,tmp)
}
rm(i)
up_counts=data.frame(up_counts)
names(up_counts)=c("condition" , "up")
up_counts <- separate(data = up_counts, col = condition, into = c("condition", "tmp"), sep = "\\_up")
up_counts=up_counts[,c(1,3)]

##Down Counts
down_counts=NULL
for( i in 1:length(diff.list.Down))
{
  tmp=cbind(names(diff.list.Down)[i],dim(diff.list.Down[[i]])[1])
  down_counts=rbind(down_counts,tmp)
}
rm(i)
down_counts=data.frame(down_counts)
names(down_counts)=c("condition" , "down")
down_counts <- separate(data = down_counts, col = condition, into = c("condition", "tmp"), sep = "\\_down")
down_counts=down_counts[,c(1,3)]

#merge all
DE_sig_counts <- merge(DE_counts,sig_counts,by = "condition")
DE_sig_up_counts <- merge(DE_sig_counts,up_counts,by = "condition")
DE_sig_up_down_counts <- merge(DE_sig_up_counts,down_counts,by = "condition")
#In case of genes
names(DE_sig_up_down_counts)=c("Condition" , "Total Expressed Genes" , "Significantly Expressed Genes" , "Up Regulated Genes" , "Down Regulated Genes")
#In case of transcripts
#names(DE_sig_up_down_counts)=c("Condition" , "Total Expressed Transcripts" , "Significantly Expressed Transcripts" , "Over Expressed Transcripts" , "Under Expressed Transcripts")
#write.table(DE_sig_up_down_counts,file = "DE_analysis_statistics.txt", quote = F,row.names = F,sep = "\t")
openxlsx::write.xlsx(DE_sig_up_down_counts,file = "DE_analysis_statistics.xlsx")

rm(down_counts,sig_counts,tmp,up_counts,DE_counts,DE_sig_counts,DE_sig_up_counts)
