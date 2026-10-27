gene_count=read.delim(file = "gene_expression.txt",header = T,sep = "\t",comment.char = "#",stringsAsFactors = F)
rownames(gene_count)=gene_count$Geneid
as.matrix(names(gene_count))
gene_count=gene_count[,c(1,11:13]
as.matrix(names(gene_count))
colnames(gene_count)=gsub(colnames(gene_count),pattern = ".Aligned.sortedByCoord.out.bam",replacement = "")
as.matrix(names(gene_count))
TMP=gene_count
TMP=TMP[,-1]

require(devtools)
lapply(names(TMP),function(x){write.table(cbind(rownames(TMP),TMP[[x]]), file = paste("output", x, sep = ""),row.names = F, quote = F, col.names = F,sep = "\t")})
getwd()

##For DESeq analysis, we need a matrix as input so,transform countdata data frame into matrix of counts with transcripts ids as its row names.##
y=as.matrix(gene_count[,c(2:4)])
rownames(y)=gene_count$Geneid
as.matrix(names(gene_count))
condition=factor(c("Control","Patient","Patient"))


##############create DESeqDataSet###########
library(DESeq2)
require(DESeq2)
dde=DESeqDataSetFromMatrix(countData = y, colData=DataFrame(condition),design = ~condition)


##Retain genes which have at least 1 read in any of the samples##
dde=estimateSizeFactors(dde)
sizeFactors(dde)

##Retain genes which have at least 5 reads in any of thesamples##
dds=dde[rowSums(counts(dde)>3)>=3,]

dde
#dim: 5678 10  
dds
#dim: 5651 10 

##Now specify the condition that should be considered as reference for differential expression analysis.##
#dds$condition<- relevel(dds$condition, ref="Control")

##Perform Differential expression (DE) analysis##
dds=DESeq(dds)
resultsNames(dds)

normcounts=counts(dde, normalized=TRUE)

#calulate corrilation matrix
corr_matrix=cor(normcounts, method = "spearman")
openxlsx::write.xlsx(corr_matrix, file = "03_Correlation_matrix.xlsx")

require(PerformanceAnalytics)
png("Correlation_plot.png", type="cairo", res = 600, units="in", width=ncol(corr_matrix), height=nrow(corr_matrix))
chart.Correlation(corr_matrix, histogram=TRUE, pch=19,method = "spearman")
dev.off()
detach("package:PerformanceAnalytics")

#Import target file with grouping information
require(readr)
target= read.table("target.txt",sep="\t",header=T)               
require(ggfortify)
require(ggplot2)
require(ggrepel)
#Draw PCA plot
tpca=t(normcounts)
tpca=data.frame(tpca)
tpca$group=target$group
png("PCA_plot.png", type="cairo", res = 300, units="in", width=14, height=10)
autoplot(prcomp(tpca[,c(1:(ncol(tpca)-1))]), data=tpca, colour='group')+
  geom_text_repel(aes(label=rownames(tpca),color=factor(group)),size=6)
dev.off()


# ##### Annotation to add with files ########
###Annotation
gff_anno=read.table("annotation.txt",sep="\t",header=T)

library(dplyr)
#RAW COUNTS AND NORMALIZED COUNTS
norm_counts=as.data.frame(normcounts)
norm_counts$Id <- rownames(norm_counts)
norm_counts=norm_counts[,c(ncol(norm_counts),1:(ncol(norm_counts)-1))]
norm_counts.anno=merge(norm_counts,gff_anno,by.x="Id",by.y="gene_id",all.x=T)
openxlsx::write.xlsx(norm_counts.anno, file = "01_Normalized_Expression_Counts.xlsx")
gene_count.anno=merge(gene_count,gff_anno,by.x="Geneid",by.y="gene_id",all.x=T)
openxlsx::write.xlsx(gene_count.anno, file = "02_Expression_Counts.xlsx")



###############################################
######## mild vs moderate. #######
###############################################
#PUT mild IN SECOND PLACE like AB-mild
Control_vs_Patient=results(dds, contrast=c("condition","Patient","Control"))
Control_vs_Patient_DE=as.data.frame(Control_vs_Patient)
as.matrix(names(Control_vs_Patient_DE))

Control_vs_Patient_DE=merge(Control_vs_Patient_DE,gff_anno,by.x="row.names",by.y="gene_id",all.x=T)
colnames(Control_vs_Patient_DE) <- sub(pattern = "Row.names",replacement = "Id", colnames(Control_vs_Patient_DE))

##add fold change
as.matrix(names(Control_vs_Patient_DE))
Control_vs_Patient_DE$FoldChange=2^(Control_vs_Patient_DE$log2FoldChange)
Control_vs_Patient_DE$FoldChange=format(round(Control_vs_Patient_DE$FoldChange, 5), nsmall = 3)
Control_vs_Patient_DE=Control_vs_Patient_DE[,c(1,2,ncol(Control_vs_Patient_DE),3:(ncol(Control_vs_Patient_DE)-1))]
as.matrix(names(Control_vs_Patient_DE))

####add normalize counts and mean of replicates###
Control_vs_Patient_DE=merge(Control_vs_Patient_DE,norm_counts,by = "Id",all.x= "T")
as.matrix(names(Control_vs_Patient_DE))
as.matrix(names(Control_vs_Patient_DE))

dir.create("Control_vs_Patient")
openxlsx::write.xlsx(Control_vs_Patient_DE, file = "Control_vs_Patient/Control_vs_Patient.Differential_Expression.xlsx")

#Significant DiffExprs
Control_vs_Patient_DE.sig=subset(Control_vs_Patient_DE, abs(log2FoldChange)>=1 & pvalue <= 0.01)
dim(Control_vs_Patient_DE.sig)
openxlsx::write.xlsx(Control_vs_Patient_DE.sig, file = "Control_vs_Patient/Control_vs_Patient.Significant.xlsx")
as.matrix(names(Control_vs_Patient_DE.sig))
#Over and Under expressed transcripts
Control_vs_Patient_up=subset(Control_vs_Patient_DE, Control_vs_Patient_DE$log2FoldChange >= 1 & Control_vs_Patient_DE$pvalue<= 0.01)
dim(Control_vs_Patient_up)
openxlsx::write.xlsx(Control_vs_Patient_up, file = "Control_vs_Patient/Control_vs_Patient.Significant_Up.xlsx")

Control_vs_Patient_down=subset(Control_vs_Patient_DE, Control_vs_Patient_DE$log2FoldChange <= -1 & Control_vs_Patient_DE$pvalue<= 0.01)
dim(Control_vs_Patient_down)
openxlsx::write.xlsx(Control_vs_Patient_down, file = "Control_vs_Patient/Control_vs_Patient.Significant_Down.xlsx")

# Make a basic volcano plot
png("Control_vs_Patient/Control_vs_Patient.Volcano_plot.png", type="cairo", res = 150, units="in", width=6, height=7)
plot(Control_vs_Patient_DE$log2FoldChange,-log10(Control_vs_Patient_DE$pvalue), main="Differential Expression",
     pch=20, col="darkgrey", xlim=c(-10,10), xlab="log2 Fold change", ylab="-log10(p value)")
abline(v=c(-1,1), h=-log10(0.01), col="hotpink", lty="dotdash")
grid (NULL,NULL, lty = 3, col = "cornsilk2") 
#color the points as per p-value
with(subset(Control_vs_Patient_DE, pvalue<=0.01 & abs(log2FoldChange)>=1), points(log2FoldChange, -log10(pvalue), pch=20, col="turquoise3"))
with(subset(Control_vs_Patient_DE, padj<=0.01 & abs(log2FoldChange)>=1), points(log2FoldChange, -log10(pvalue), pch=20, col="tomato2"))
legend(-10,3, c("adj.p<=0.01","p value<=0.01"), xjust=0, box.lty = "blank", text.col = c("tomato2", "turquoise3"), cex = 0.7)
dev.off()


require("ggrepel") #Avoid overlapping labels
require("plyr")
as.matrix(names(Control_vs_Patient_DE))
mutateddf <- mutate(Control_vs_Patient_DE,Legend=ifelse(padj<=0.01 & abs(log2FoldChange)>=1 & pvalue<=0.01 , "adj.p<=0.01",
                            ifelse(padj>=0.01 & pvalue<=0.01 & abs(log2FoldChange)>=1, "p value<=0.01" , "Not Significant")))

input <- cbind(gene=rownames(mutateddf ), mutateddf ) #convert the rownames to a column
xmini=min(Control_vs_Patient_DE$log2FoldChange)-4
xmaxi=max(Control_vs_Patient_DE$log2FoldChange)
volc = ggplot(na.omit(input), aes(log2FoldChange, -log10(pvalue))) + #volcanoplot with log2Foldchange versus pvalue
  geom_point(aes(col=Legend)) + #add points colored by significance
  scale_color_manual(values=c("grey","turquoise3")) +
  scale_x_continuous(limits = c(-10,10)) +
  geom_hline(yintercept = -log10(0.01),size=0.3, color="hotpink", linetype="dashed") +
  geom_vline(xintercept = c(-1,1),size=0.3, color="hotpink", linetype="dashed") +
  ggtitle("Differential Expression") +
  theme(plot.title = element_text(hjust = 0.5))
input.pos=subset(input, log2FoldChange>=1 & padj <= 0.01)
input.neg=subset(input, log2FoldChange<=-1 & padj <= 0.01)
volc=volc+geom_text_repel(data=head(input.pos[order(input.pos$pvalue,decreasing = F),],25), aes(label=Id),size = 3, color = 'black', box.padding = 0.1, point.padding = 0.1, segment.color = 'grey50') #adding text for the top 10  genes
volc=volc+geom_text_repel(data=head(input.neg[order(input.neg$pvalue,decreasing = F),],25), aes(label=Id),size = 3, color = 'black', box.padding = 0.1, point.padding = 0.1, segment.color = 'grey50') #adding text for the top 10 genes
ggsave(volc, filename = "Control_vs_Patient/Control_vs_Patient.Volcano_plot.genename.png",height = 10,width = 10,device = "png",units = "in",dpi = 300)
#ggsave(volc, filename = "mild_vs_ RMTP10_1/mild_vs_ RMTP10.Volcano_plot.genename.png")
rm(list="xmini", "xmaxi" , "input.pos", "input.neg","mutateddf","input")


###heatmap####
Control_vs_Patient_DE.sig.norm= normcounts[rownames(normcounts) %in% Control_vs_Patient_DE.sig$Id,]
dim(Control_vs_Patient_DE.sig.norm)
Control_vs_Patient_DE.sig.norm=data.frame(Control_vs_Patient_DE.sig.norm)
as.matrix(names(Control_vs_Patient_DE.sig.norm))

require(pheatmap)
png("Control_vs_Patient/Control_vs_Patient.Heatmap_R.png", type="cairo", res = 200, units="in", width=8, height=25)
anno=target[,2]
rownames(anno)=colnames(data.frame(Control_vs_Patient_DE.sig.norm))
colnames(anno)="Group"
newCols <- colorRampPalette(colorspace::heat_hcl(length(unique(anno$Group))))
mycolors <- newCols(length(unique(anno$Group)))
names(mycolors) <- unique(anno$Group)
mycolors <- list(Group = mycolors)
pheatmap(data.frame(Control_vs_Patient_DE.sig.norm),color = colorRampPalette(c("green4", "yellow", "firebrick2"))(50),
         show_rownames = F,show_colnames = T,cluster_cols = T,scale = "row",annotation= anno,annotation_colors = mycolors,
         border_color = F,labels_col = colnames(Control_vs_Patient_DE.sig.norm),fontsize_col = 8,fontsize = 6,angle_col = 0,cluster_rows = T)
dev.off()

rm(Control_vs_Patient_DE.sig.norm)

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
keytypes(org.Hs.eg.db)

#edit gene_id column(remove "." part from ids (ENSG00000001626.10 as ENSG00000001626))
require(tidyr)
universelist <- diff.merged[[1]]$ENTREZ

#universelist=unique(diff.merged[1])
universelist=as.character(unique(universelist))
for(i in 1:length(diff.list.sig)){
  genelist=unique(as.character(diff.list.sig[[i]]$ENTREZ))
  enrichGO(gene=genelist, universe = universelist, OrgDb = org.Hs.eg.db, pvalueCutoff = 0.05,
           ont = "BP", pAdjustMethod = "BH", minGSSize = 2, keyType = "ENTREZID", qvalueCutoff = 0.05,
           readable=T) -> tmp
  if(class(tmp)=="NULL"){next}else{
    assign(paste(names(diff.list.sig[i]), "BP", sep="_"),
           slot(simplify(tmp, cutoff=0.7, by="p.adjust", select_fun=min),"result"))
  }
  rm(tmp)
  rm("genelist", "genelist1")
}
rm(i)

for(i in 1:length(diff.list.sig)){
  genelist=as.character(unique(diff.list.sig[[i]]$ENTREZ))
  enrichGO(gene=genelist1, universe = universelist, OrgDb = org.Hs.eg.db, pvalueCutoff = 0.05,
           ont = "CC", pAdjustMethod = "BH", minGSSize = 2, keyType = "ENTREZID", qvalueCutoff = 0.05,
           readable=T) -> tmp
  if(class(tmp)=="NULL"){next}else{
    assign(paste(names(diff.list.sig[i]), "CC", sep="_"),
           slot(simplify(tmp, cutoff=0.7, by="p.adjust", select_fun=min),"result"))
  }
  rm(tmp)
  rm("genelist", "genelist1")
}
rm(i)

for(i in 1:length(diff.list.sig)){
  genelist=unique(as.character(diff.list.sig[[i]]$ENTREZ))
  enrichGO(gene=genelist, universe = universelist, OrgDb = org.Hs.eg.db, pvalueCutoff = 0.05,
           ont = "MF", pAdjustMethod = "BH", minGSSize = 2, keyType = "ENTREZID", qvalueCutoff = 0.05,
           readable=T) -> tmp
  if(class(tmp)=="NULL"){next}else{
    assign(paste(names(diff.list.sig[i]), "MF", sep="_"),
           slot(simplify(tmp, cutoff=0.7, by="p.adjust", select_fun=min),"result"))
  }
  rm(tmp)
  rm("genelist", "genelist1")
}
rm(i)

allGO=mget(ls(pattern = "BP|MF|CC"))
rm(list = ls(pattern = "BP|MF|CC"))
#insert Category in GO output
for(i in 1:length(allGO)){
  if(grepl("BP",names(allGO[i]))==T){
    allGO[[i]]$Category=c(rep("BP", dim(allGO[[i]])[1]))
  }else if(grepl("MF",names(allGO[i]))==T){
    allGO[[i]]$Category=c(rep("MF", dim(allGO[[i]])[1]))
  }else if(grepl("CC",names(allGO[i]))==T){
    allGO[[i]]$Category=c(rep("CC", dim(allGO[[i]])[1]))
  }
}
rm(i)

#combine BP,MF,CC of each comparions to export GO files
compariosions=unique(do.call(rbind,str_split(names(allGO),"_DE"))[,1])
for(i in 1:length(compariosions)){
  assign(paste(compariosions[i],"_GO_export",sep="_"),
         do.call(rbind,allGO[grep(compariosions[i], names(allGO))])[,c(10,1:6,8,9)])
}
rm(i)
GO.export=mget(ls(pattern = "_GO_export"))
rm(list=ls(pattern = "_GO_export"))
#export GO files
for(i in 1:length(GO.export)){
  openxlsx::write.xlsx(x = GO.export[[i]],file = paste0("./",str_split(names(GO.export)[i], pattern = "__")[[1]][1],
                                                        "/",str_split(names(GO.export)[i], pattern = "__")[[1]][1],
                                                        ".Gene_Ontology_Enrichment.xlsx"))
}
rm(i)


#combine BP,MF,CC of each comparions
compariosions=unique(do.call(rbind,str_split(names(allGO),"_DE"))[,1])
for(i in 1:length(compariosions)){
  assign(paste(compariosions[i],"combined_GO",sep="_"),
         do.call(rbind,allGO[grep(compariosions[i], names(allGO))])[,c(10,1,2,6,8)])
}
rm(i)
combi.GO=mget(ls(pattern = "_combined_GO"))
rm(list=ls(pattern = "_combined_GO"))
#rename colums in GO
for(i in 1:length(combi.GO)){
  names(combi.GO[[i]])=paste(c("category", "ID", "term", "adj_pval", "genes"))
  combi.GO[[i]]$genes=gsub("/", ",", combi.GO[[i]]$genes)
}

genelist.GO=diff.list.sig
as.matrix(names(genelist.GO[[1]]))
#generte gene list files
for(i in 1:length(genelist.GO)){
  genelist.GO[[i]]=genelist.GO[[i]][,c(1,4)]
  names(genelist.GO[[i]])=paste(c("ID", "logFC"))
}
rm(i)

#create Circ_dat
for(i in 1:length(compariosions)){
  assign(paste(compariosions[i], "circ",sep="_"),
         circle_dat(combi.GO[[grep(compariosions[i], names(combi.GO[i]))]], genelist.GO[[grep(compariosions[i], names(genelist.GO[i]))]])
         
  )
}
rm(i)
GO.circ=mget(ls(pattern = "_circ"))
rm(list = ls(pattern = "_circ"))


#Make bubble plot
for(i in 1:length(GO.circ)){
  data.circ=GO.circ[[i]][!duplicated(GO.circ[[i]]$term),]
  data.circ2=subset(data.circ, data.circ$count>1)
  sub2.data.circ=subset(data.circ2, data.circ2$adj_pval<=0.05)
  rang_data.circ <- c(1, 10)
  cols_data.circ <- c("chartreuse4", "brown2", "cornflowerblue")
  data.circ3=sub2.data.circ
  bubplot=ggplot(data.circ3, aes(zscore, -log10(adj_pval), (alpha=0.1), fill=category, size = count)) +
    labs(title = "GO bubble plot", x = "z-score", y = "-log10(adj p-value)") +
    theme(plot.title = element_text(hjust = 0.5))+
    geom_hline(yintercept = -log10(0.05), col = "orange") + 
    scale_size(range = c(min(data.circ3$zscore)-5,max(data.circ3$zscore)+5)) +
    geom_point(shape = 21, colour="black", stroke=0.1)+
    geom_text(data = data.circ3, aes(zscore, -log10(adj_pval), label=term), size = 5, check_overlap = T)+
    theme(legend.position = "bottom", legend.direction = "horizontal")+
    theme(legend.box = "horizontal",legend.key.size = unit(1, "cm"))+
    labs(x = "z-score", y = "-log10(adj p-value)")+
    xlim(c(min(data.circ3$zscore)-1,max(data.circ3$zscore)+1))
  ggsave(filename = paste0("./",str_split(names(GO.circ[i]),"_circ")[[1]][1],"/",str_split(names(GO.circ[i]),"_circ")[[1]][1],".Bubble_plot1.png"),
         plot = bubplot,device = "png", dpi = 300, height = 10, width = 15, units = "in"
  )
  rm(list=c("bubplot","data.circ", "data.circ2","sub2.data.circ","xmini","xmaxi"))
}
rm(i)

require(pathview)
setwd(workdir)
for(i in 1:length(diff.list.sig))
{
  entrez.list=as.character(unique(diff.list.sig[[i]]$ENTREZ))
  entrez.list.universe=as.character(unique(diff.merged[[i]]$ENTREZ))
  kegg.results=enrichKEGG(gene=entrez.list, organism = "hsa", keyType = "kegg", pvalueCutoff = 0.05,
                          pAdjustMethod = "BH", universe=entrez.list.universe, minGSSize = 2,
                          qvalueCutoff = 0.05, use_internal_data = FALSE)
  if(isTRUE(class(kegg.results)=="NULL")) { next } else kegg.out=data.frame(kegg.results@result)
  openxlsx::write.xlsx(kegg.out[,c(1:6,8,9)], file=paste0("./",str_split(names(diff.list.sig)[i],"_DE")[[1]][1],"/",str_split(names(diff.list.sig)[i],"_DE")[[1]][1],".KEGG_Pathway_Enrichment1.xlsx"))
  
  genelist=diff.list.sig[[i]]$log2FoldChange
  names(genelist)=diff.list.sig[[i]]$ENTREZ
  top10kegg=kegg.out$ID[1:10]
  dir.create(paste0(str_split(names(diff.list.sig)[i], pattern = "_DE")[[1]][1],"/KEGG_pathways"))
  setwd(paste0(str_split(names(diff.list.sig)[i], pattern = "_DE")[[1]][1],"/KEGG_pathways","/"))
  for (j in 1:10){
    pathview(gene.data=genelist, pathway.id = top10kegg[j], species = "hsa",gene.idtype="KEGG",
             out.suffix=paste0(str_split(names(diff.list.sig)[i], pattern = "_DE")[[1]][1]),
             limit = 10)
  }
  rm(list = c("keys1","entrez.list","keys2","entrez.list.universe","Kegg_gene.list","kegg_universal.list",
              "kegg.results","kegg.out","tmp.sig.list","tmp.sig.list.with.entrezid","genelist","top10kegg"))
}
