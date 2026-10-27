#test of 16srRna pipeline
setwd(getwd())

library(dplyr)
#Reaquired file: 1.final.tx.1.cons.taxonomy(Got from mothur analysis) 2.final.tx.shared
#Rename the file "final.tx.1.cons.taxonomy"
dir.create("Files")
file.copy("final.tx.1.cons.taxonomy.NEW","Files/final.tx.1.cons.taxonomy")
file.copy("final.tx.shared.NEW","Files/final.tx.shared")
setwd("Files")
file.rename("final.tx.1.cons.taxonomy","final.tx.1.cons.taxonomy.NEW")

final.tx.1.cons.taxonomy.NEW=read.csv("final.tx.1.cons.taxonomy.NEW",sep = "\t")
#Find the "unknown" Otu in "final.tx.1.cons.taxonomy" file
#Note: note down the unknown Otu ID same ID u have to remove from "final.tx.shared.NEW.transposed" file
final.tx.1.cons.taxonomy.NEW1<-final.tx.1.cons.taxonomy.NEW[!grepl("unknown",final.tx.1.cons.taxonomy.NEW$Taxonomy), ]
#removing the selected rows and assigning to new data frame
write.table(final.tx.1.cons.taxonomy.NEW1,file = "final.tx.1.cons.taxonomy.NEW",sep = "\t",row.names = FALSE)


##########
#######
#####
###

# Alpha diversity analysis
require(phyloseq)
require(ggplot2)
require(pheatmap)
library(ape)
library(dplyr)
library(ggplot2)
library(gplots)
library(lme4)
library(phangorn)
library(phyloseq)
library(plotly)
library(tidyr)
library(VennDiagram)
library(phyloseq)
#library(venneuler)
require(vegan)
require(reshape)


# Load the final.tx.shared file convert into table form

shared=read.delim("final.tx.shared",header = T,stringsAsFactors = F, sep = "\t")
selectcols<-final.tx.1.cons.taxonomy.NEW1$OTU
a<-select(shared,c(1:3,selectcols))
shared_new=data.frame(t(a))
write.table(shared_new, file = "final.tx.shared.NEW.transposed",sep = "\t",row.names = T,col.names = F,quote = F)
#rm(shared)

#file for merging with the taxonomy file
b<-select(shared,c(2,selectcols))
b_transposed <- cbind(rownames(t(b)),(data.frame(t(b), row.names=NULL)))
# changing the first row as header
names(b_transposed) <- b_transposed[1,]
#deleting that row
b_transposed <- b_transposed[-1,]
#Renaming the 1 st col
names(b_transposed)[1]<-"OTU"

#############
###########
########
######
####
#
#Merge Taxonomy file and sample otu's into new dataframe, on the basis of same OTU Id's

merged_df<-merge(final.tx.1.cons.taxonomy.NEW,b_transposed,by = "OTU")
#Store the _mothur__OTU__Table.tab
merged_df$Taxonomy=gsub("\\(100)", "", merged_df$Taxonomy)
#otus_norm$Taxonomy=gsub("\\;", "|", otus_norm$Taxonomy)
merged_df$Taxonomy=gsub("k__", "", merged_df$Taxonomy)
merged_df$Taxonomy=gsub("\\p__", "", merged_df$Taxonomy)
merged_df$Taxonomy=gsub("\\c__", "", merged_df$Taxonomy)
merged_df$Taxonomy=gsub("\\o__", "", merged_df$Taxonomy)
merged_df$Taxonomy=gsub("f__", "", merged_df$Taxonomy)
merged_df$Taxonomy=gsub("\\g__", "", merged_df$Taxonomy)
merged_df$Taxonomy=gsub("s__", "", merged_df$Taxonomy)

ProjectID="CG_DN_500"
write.table(merged_df,paste(ProjectID,'_mothur_OTU_Table.tab',sep = ""), sep = "\t",row.names = FALSE)

############
## Normalized the OTUs on the basis of all samples

otus=read.delim(file = paste(ProjectID,'_mothur_OTU_Table.tab',sep = ""),header = T,sep = "\t",comment.char = "#",stringsAsFactors = F)
as.matrix(names(otus))
otus_norm=otus
for( i in 4:length(otus_norm))
{
  otus_norm[,i]=round((otus_norm[,i]/sum(otus_norm[,i]))*100000,digits = 1)
}
write.table(otus_norm, file = paste(ProjectID,"_OTU_Table.normalized.tab",sep=""),sep = "\t",row.names = F,col.names = T,quote = F)

##############
#Relative abundance
otus=read.delim(file = paste(ProjectID,"_mothur_OTU_Table.tab",sep=""),header = T,sep = "\t",comment.char = "#",stringsAsFactors = F)
as.matrix(names(otus))
otus_norm=otus
for( i in 4:length(otus_norm))
{
  otus_norm[,i]=round((otus_norm[,i]/sum(otus_norm[,i]))*100,digits = 2)
}
write.table(otus_norm, file = paste(ProjectID,"_mothur_OTU_Table.relative-abundance.tab",sep=""),sep = "\t",row.names = F,col.names = T,quote = F)


##relative abundance for LEfSe input file
otus=read.delim(file = paste(ProjectID,"_mothur_OTU_Table.tab",sep=""),header = T,sep = "\t",comment.char = "#",stringsAsFactors = F)
as.matrix(names(otus))


otus_norm=otus
for( i in 4:length(otus_norm))
{
  otus_norm[,i]=round((otus_norm[,i]/sum(otus_norm[,i]))*100,digits = 2)
}
write.table(otus_norm, file = paste(ProjectID,"_mothur_OTU_Table.relative-abundance.LEfSe-file.tab",sep=""),sep = "\t",row.names = F,col.names = T,quote = F)


##########
#####Normalized OTUs
#Input for STAMP Analysis
otus=read.delim(file = paste(ProjectID,"_mothur_OTU_Table.tab",sep=""),header = T,sep = "\t",comment.char = "#",stringsAsFactors = F)
as.matrix(names(otus))
otus_norm=otus
for( i in 4:length(otus_norm))
{
  otus_norm[,i]=round((otus_norm[,i]/sum(otus_norm[,i]))*100000,digits = 1)
}

normalized_OTUs=otus_norm
rownames(normalized_OTUs)=normalized_OTUs$OTU
#For STAMP input
as.matrix(names(normalized_OTUs))
normalized_OTUs1=normalized_OTUs[,c(3:ncol(normalized_OTUs))]
colnames(normalized_OTUs1)=sub(pattern = "Taxonomy",replacement = "Feature", colnames(normalized_OTUs1))
normalized_OTUs1$Feature=gsub("\\(100)", "", normalized_OTUs1$Feature)
write.table(normalized_OTUs1, file = "Normalized.OTU.spf",sep = "\t",row.names = F,col.names = T,quote = F)


#############
######
###
#Splitting OTU_table file taxonomy column into all features
splitted_otu=separate(otus,Taxonomy,c("Kingdom","Phylum","Class","Order","Family","Genus","Species"),sep = "\\;")

write.table(splitted_otu,paste(ProjectID,"_mothur_OTU_Table.splitted.tab",sep=""),sep="\t",row.names = FALSE)



################
#############
########
######
###
#
###make Krona input files
#Extract all taxonomy features and saving into krona file format to visualize in pie-chart form
rawOTUs=read.delim(file = paste(ProjectID,"_mothur_OTU_Table.splitted.tab",sep=""),header = T,sep = "\t",comment.char = "#",stringsAsFactors = F)
as.matrix(names(rawOTUs))
column_names=colnames(rawOTUs)
dir.create("Krona_pie-chart")
for( i in 9:length(rawOTUs))
{
  filename=paste("Krona_pie-chart/",column_names[i],"_krona.in",sep = "")
  write.table(rawOTUs[,c(i,3:9)], file = filename,sep = "\t",row.names = F,col.names = F,quote = F)
}
rm(rawOTUs)

#Perl scripts to run krona chart
setwd("Krona_pie-chart")
getwd()
cmd=paste("perl /home/clevergene/anaconda3/opt/krona/scripts/ImportText.pl *in")
system(cmd)

setwd("..")
#########
#####
##
#
# ALPHA DIVERSITY
APS_mothur=import_mothur(mothur_shared_file = "final.tx.shared",mothur_constaxonomy_file = "final.tx.1.cons.taxonomy.NEW")
colnames(tax_table(APS_mothur))=c("Kingdom","Phylum", "Class", "Order", "Family", "Genus","Species")
sampleID=colnames(APS_mothur@otu_table)
sampleID


sampleType=c(sampleID)
sampletable=data.frame(cbind(sampleID,sampleType),stringsAsFactors = FALSE)
colnames(sampletable)=c("SampleID","SampleType")
rownames(sampletable)=sampleID
sample_data(APS_mothur)=sampletable


############
########
#####
#PCoA PLOT New with OTUs samples
png("PCoA_plots.png", type="cairo", res = 250, units="in", width=7, height=6)
APS_mothur_ordinates_pca=ordinate(APS_mothur,"PCoA")
#PCA=plot_ordination(physeq = APS_mothur,ordination = res_mothur_ordinates_pca,type = "samples",shape = "SampleType",color = "SampleType",title = "PCoA")
PCA=plot_ordination(physeq = APS_mothur,ordination = APS_mothur_ordinates_pca,type = "samples",shape = NULL, color = "SampleType",title = "PCoA",label = "SampleID")
PCA +geom_point(size=1.5)
PCA +geom_label(label.size=2.5)
dev.off()




####################################################
############ ALPHA DIVERSITY #######################
####################################################
#shared file
APS_mothur=import_mothur(mothur_shared_file = "final.tx.shared",mothur_constaxonomy_file = "final.tx.1.cons.taxonomy.NEW")
colnames(tax_table(APS_mothur))=c("Kingdom","Phylum", "Class", "Order", "Family", "Genus", "species")
sampleID=colnames(APS_mothur@otu_table)
#sampleType=c("Sample_2669_18A" , "Sample_2669_18B" , "Sample_2670_18A" , "Sample_2670_18B" , "Sample_2672_18A" , "Sample_2672_18B" , "Sample_2673_18A" , "Sample_2673_18B" , "Sample_2984_18A" , "Sample_2984_18B" , "Sample_3773_18A" , "Sample_3773_18B" , "Sample_3776_18A" , "Sample_3776_18B" , "Sample_3777_18A" , "Sample_3777_18B" , "Sample_3778_18A" , "Sample_3778_18B" , "Sample_3779_18A" , "Sample_3779_18B")
#sampleType=c('P1R1','P1R2','P1R3','P2R1','P2R2','P2R3')
sampleType=sampleID
sampletable=data.frame(cbind(sampleID,sampleType),stringsAsFactors = FALSE)
colnames(sampletable)=c("SampleID","SampleType")
rownames(sampletable)=sampleID
sample_data(APS_mothur)=sampletable

#PCoA PLOT
#PCoA PLOT New without OTUs only samples
png("PCoA_plot.png", type="cairo", res = 250, units="in", width=7, height=6)
APS_mothur_ordinates_pca=ordinate(APS_mothur,"PCoA")
#PCA=plot_ordination(physeq = APS_mothur,ordination = res_mothur_ordinates_pca,type = "samples",shape = "SampleType",color = "SampleType",title = "PCoA")
PCA=plot_ordination(physeq = APS_mothur,ordination = APS_mothur_ordinates_pca,type = "samples",shape = NULL,color = "SampleID",title = "PCoA")
PCA +geom_point(size=1.5)
###PCA +theme_bw() +geom_text(size = 10, vjust = 1.5) +theme(text = element_text(size = 16))
PCA +geom_text_repel(aes(label=PCA$data$SampleID),size=4)
dev.off()


#### AlphaDiversity plot####
png("AlphaDiversity_plot.Group-wise.png", type="cairo", res = 250, units="in", width=15, height=10)
AlphaDiversity=plot_richness(APS_mothur,color = "SampleType",x="SampleType")
AlphaDiversity +theme(axis.text.x = element_text(size = 6.5, angle = 90, hjust = 1))
dev.off()

png("AlphaDiversity_plot.Chao1_ACE_Shannon.sample-wise.png", type="cairo", res = 250, units="in", width=14, height=10)
AlphaDiversity=plot_richness(APS_mothur,color = "SampleType",measures=c("Chao1", "ACE", "Shannon"))
AlphaDiversity +theme(axis.text.x = element_text(size = 7, angle = 90, hjust = 1))
dev.off()

png("AlphaDiversity_plot.Simpson_InvSimpson_Fisher.sample-wise.png", type="cairo", res = 250, units="in", width=12, height=10)
AlphaDiversity=plot_richness(APS_mothur,color = "SampleType",measures=c("Simpson", "InvSimpson", "Fisher"))
AlphaDiversity +theme(axis.text.x = element_text(size = 5, angle = 90, hjust = 1))
dev.off()

#table
AlphaDiversity_table=estimate_richness(APS_mothur)
write.table(AlphaDiversity_table, file = "AlphaDiversity_table.txt",sep = "\t",row.names = T,col.names = T,quote = F)

### NEW rarefraction plot
png("Rarefraction_plot.png", type="cairo", res = 200, units="in", width=14, height=10)
rarecurve(t(otu_table(APS_mothur)))
dev.off()



##############################################
###table,barplot and heatmap phylum level#####
##############################################
#Table
APS_mothur_phylum=tax_glom(APS_mothur,taxrank = "Phylum")
#phylum=cbind(data.frame(APS_mothur_phylum@tax_table),data.frame(APS_mothur_phylum@otu_table))
phylum=cbind(APS_mothur_phylum@tax_table,APS_mothur_phylum@otu_table)
phylum=data.frame(phylum,stringsAsFactors = F)
phylum$OTU=rownames(phylum)
as.matrix(names(phylum))
phylum=phylum[,c(ncol(phylum),1,2,8:(ncol(phylum)-1))]
as.matrix(names(phylum))
#removing unclassified 
phylum_filtered=filter(phylum,!grepl("*_unclassified",Phylum))
write.table( APS_mothur_phylum10.merge, file="Phylum10_level.txt", sep = "\t",row.names = F,col.names = T, quote = F)
APS_mothur_phylum10.merge=read.delim("Phylum_level10.txt",header = T,stringsAsFactors = F, sep = "\t")

#Bar Plot
#APS_mothur_phylum10 = names(sort(taxa_sums(APS_mothur_phylum), TRUE)[1:10]))
#install.packages("reshape")
library(reshape)
#APS_mothur_phylum filter by removal of unclassified Phylum
APS_mothur_phylum_Filtered = subset_taxa(APS_mothur_phylum,!grepl("*_unclassified",Phylum))
APS_mothur_phylum10 = data.frame(sort(taxa_sums(APS_mothur_phylum_Filtered), TRUE)[1:10])
APS_mothur_phylum10$OTU=rownames(APS_mothur_phylum10)
APS_mothur_phylum10=APS_mothur_phylum10[,c(2,1)]
colnames(APS_mothur_phylum10)=c("OTU","Size")
APS_mothur_phylum10.merge=merge(APS_mothur_phylum10,phylum,by = "OTU",all.x = T)
APS_mothur_phylum10.merge=APS_mothur_phylum10.merge[rev(order(APS_mothur_phylum10.merge$Size)),]
as.matrix(names(APS_mothur_phylum10.merge))
names(APS_mothur_phylum10.merge)[c(5:ncol(APS_mothur_phylum10.merge))]<-c("GR-6","HR-6","KR-1","MR-6","HDR-6","KR-6","VC-1","VS-6")
APS_mothur_phylum10.merge=APS_mothur_phylum10.merge[,c(4:ncol(APS_mothur_phylum10.merge))]
write.table( APS_mothur_phylum10.merge, file="Phylum10_level.txt", sep = "\t",row.names = F,col.names = T, quote = F)
read.table("Phylum10_level.txt", header = T,stringsAsFactors = F, sep = "\t" )
APS_mothur_phylum10.merge=read.delim("Phylum10_level.txt",header = T,stringsAsFactors = F, sep = "\t")
#APS_mothur_phylum10.merge$Athul_MetG1 <- as.numeric(as.character(APS_mothur_phylum10.merge$Athul_MetG1))
cols.num <- names(APS_mothur_phylum10.merge[,c(2:ncol(APS_mothur_phylum10.merge))])
APS_mothur_phylum10.merge[cols.num] <- sapply(APS_mothur_phylum10.merge[cols.num],as.numeric)

for( i in 2:length(APS_mothur_phylum10.merge))
{
  APS_mothur_phylum10.merge[,i]=round((APS_mothur_phylum10.merge[,i]/sum(APS_mothur_phylum10.merge[,i]))*100,digits = 1)
}
rm(i)
APS_mothur_phylum10.merge_melt=melt(APS_mothur_phylum10.merge)
colnames(APS_mothur_phylum10.merge_melt)=c("Phylum","Sample","Percentage")
png("Top10_normalized_OTUs_Phylum_Barplot.png", type="cairo", res = 200, units="in", width=10, height=6)
phylum_colors=rainbow(length(unique(APS_mothur_phylum10.merge$Phylum)), s=.6, v=.9)[sample(1:length(unique(APS_mothur_phylum10.merge$Phylum)),length(unique(APS_mothur_phylum10.merge$Phylum)))]
ggplot(APS_mothur_phylum10.merge_melt,aes(x = Sample,y = Percentage,fill=Phylum))+
  geom_bar(stat = "identity") + 
  scale_fill_manual(values = phylum_colors)+
  theme(axis.text.x = element_text(angle = 0, hjust = 1))
dev.off()


#Heatmap
APS_mothur_phylum20 = data.frame(sort(taxa_sums(APS_mothur_phylum_Filtered), TRUE)[1:20])
APS_mothur_phylum20$OTU=rownames(APS_mothur_phylum20)
APS_mothur_phylum20=APS_mothur_phylum20[,c(2,1)]
colnames(APS_mothur_phylum20)=c("OTU","Size")
APS_mothur_phylum20.merge=merge(APS_mothur_phylum20,phylum,by = "OTU",all.x = T)
APS_mothur_phylum20.merge=APS_mothur_phylum20.merge[rev(order(APS_mothur_phylum20.merge$Size)),]
as.matrix(names(APS_mothur_phylum20.merge))
APS_mothur_phylum20.merge=APS_mothur_phylum20.merge[,c(4:ncol(APS_mothur_phylum20.merge))]
#APS_mothur_phylum20.merge$Athul_MetG1 <- as.numeric(as.character(APS_mothur_phylum20.merge$Athul_MetG1))
cols.num <- names(APS_mothur_phylum20.merge[,c(2:ncol(APS_mothur_phylum20.merge))])
APS_mothur_phylum20.merge[cols.num] <- sapply(APS_mothur_phylum20.merge[cols.num],as.numeric)

for( i in 2:length(APS_mothur_phylum20.merge))
{
  APS_mothur_phylum20.merge[,i]=round((APS_mothur_phylum20.merge[,i]/sum(APS_mothur_phylum20.merge[,i]))*100,digits = 1)
}
rm(i)
write.table( APS_mothur_phylum20.merge, file="Phylum20_level.txt", sep = "\t",row.names = F,col.names = T, quote = F)
read.table("Phylum20_level.txt", header = T,stringsAsFactors = F, sep = "\t" )
APS_mothur_phylum20.merge=read.delim("Phylum20_level.txt",header = T,stringsAsFactors = F, sep = "\t")
names(APS_mothur_phylum20.merge)[c(2:ncol(APS_mothur_phylum20.merge))]<-c("GR-6","HR-6","KR-1","MR-6","HDR-6","KR-6","VC-1","VS-6")
rownames(APS_mothur_phylum20.merge)=APS_mothur_phylum20.merge$Phylum
APS_mothur_phylum20.merge=APS_mothur_phylum20.merge[,-1]
anno=data.frame(factor(c(rep("A",25),rep("B",12),rep("C",6),rep("D",3))))
#rownames(anno)=colnames(APS_mothur_phylum20.merge)
#rownames(anno)=colnames("P1R1", "P1R2", "P1R3", "P2R1", "P2R2", "P2R3")
colnames(anno)="Group"

png("Top20_normalized_OTUs_Phylum_heatmap.png", type="cairo", res = 200, units="in", width=7, height=4)
# pheatmap(APS_mothur_phylum20.merge,show_rownames = T,show_colnames = T,cluster_cols = F,scale = "none",annotation_col = anno,border_color = F,fontsize_col = 7,fontsize_row = 8,fontsize = 7,cellheight = 10)
pheatmap(APS_mothur_phylum20.merge,show_rownames = T,show_colnames = T,cluster_cols = T,scale = "none",border_color = F,fontsize_col = 7,fontsize_row = 8,fontsize = 7,angle_col = 0)
dev.off()




# ##############################################
# ###table,barplot and heatmap species level#####
# ##############################################
# #Table
APS_mothur_species=tax_glom(APS_mothur,taxrank = "species")
species=cbind(data.frame(APS_mothur_species@tax_table),data.frame(APS_mothur_species@otu_table))
species=cbind(APS_mothur_species@tax_table,APS_mothur_species@otu_table)
species=data.frame(species,stringsAsFactors = F)
species$OTU=rownames(species)
as.matrix(names(species))
species=species[,c(ncol(species),1:(ncol(species)-1))]
#removing unclassified 
species_filtered=filter(species,!grepl("*_unclassified",species))
write.table(species_filtered, file="Species_level.txt", sep = "\t",row.names = F,col.names = T, quote = F)


# 
# #Bar Plot
APS_mothur_species_Filtered = subset_taxa(APS_mothur_species,!grepl("*_unclassified",species))
APS_mothur_species10 = data.frame(sort(taxa_sums(APS_mothur_species_Filtered), TRUE)[1:10])
APS_mothur_species10$OTU=rownames(APS_mothur_species10)
APS_mothur_species10=APS_mothur_species10[,c(2,1)]
colnames(APS_mothur_species10)=c("OTU","Size")
APS_mothur_species10.merge=merge(APS_mothur_species10,species_filtered,by = "OTU",all.x = T)
APS_mothur_species10.merge=APS_mothur_species10.merge[rev(order(APS_mothur_species10.merge$Size)),]
as.matrix(names(APS_mothur_species10.merge))
APS_mothur_species10.merge=APS_mothur_species10.merge[,c(7:ncol(APS_mothur_species10.merge))]
cols.num <- names(APS_mothur_species10.merge[,c(4:ncol(APS_mothur_species10.merge))])
APS_mothur_species10.merge[cols.num] <- sapply(APS_mothur_species10.merge[cols.num],as.numeric)

for( i in 4:length(APS_mothur_species10.merge))
{
  APS_mothur_species10.merge[,i]=round((APS_mothur_species10.merge[,i]/sum(APS_mothur_species10.merge[,i]))*100,digits = 1)
}
rm(i)

#concatination of two columns and making in to 1 new column
APS_mothur_species10.merge$Genus=do.call(rbind,strsplit(APS_mothur_species10.merge$Genus,split = "__"))[,2]
APS_mothur_species10.merge$species=do.call(rbind,strsplit(APS_mothur_species10.merge$species,split = "__"))[,2]
APS_mothur_species10.merge=cbind(newcol=paste(APS_mothur_species10.merge$Genus,APS_mothur_species10.merge$species),APS_mothur_species10.merge)
library(ggplot2)
library(gplots)
names(APS_mothur_species10.merge)[c(5:ncol(APS_mothur_species10.merge))]<-c("GR-6","HR-6","KR-1","MR-6","HDR-6","KR-6","VC-1","VS-6")
APS_mothur_species10.merge_melt=melt(APS_mothur_species10.merge[c(1,5:ncol(APS_mothur_species10.merge))])
colnames(APS_mothur_species10.merge_melt)=c("Species","Sample","Percentage")

png("Top10_normalized_OTUs_species_Barplot.png", type="cairo", res = 200, units="in", width=9, height=7)
species_colors=rainbow(length(unique(APS_mothur_species10.merge_melt$Species)), s=.6, v=.9)[sample(1:length(unique(APS_mothur_species10.merge_melt$Species)),length(unique(APS_mothur_species10.merge_melt$Species)))]
ggplot(APS_mothur_species10.merge_melt,aes(x = Sample,y = Percentage,fill=Species))+
  geom_bar(stat = "identity") + 
  scale_fill_manual(values = species_colors)+
  theme(axis.text.x = element_text(angle = 0, hjust = 1))
dev.off()

# #Heatmap
APS_mothur_species20 = data.frame(sort(taxa_sums(APS_mothur_species_Filtered), TRUE)[1:20])
APS_mothur_species20$OTU=rownames(APS_mothur_species20)
APS_mothur_species20=APS_mothur_species20[,c(2,1)]
colnames(APS_mothur_species20)=c("OTU","Size")
APS_mothur_species20.merge=merge(APS_mothur_species20,species,by = "OTU",all.x = T)
APS_mothur_species20.merge=APS_mothur_species20.merge[rev(order(APS_mothur_species20.merge$Size)),]
as.matrix(names(APS_mothur_species20.merge))
APS_mothur_species20.merge=APS_mothur_species20.merge[,c(7:ncol(APS_mothur_species20.merge))]
#cols.num <- names(APS_mothur_species20.merge[,c(4:ncol(APS_mothur_species20.merge))])
APS_mothur_species20.merge[cols.num] <- sapply(APS_mothur_species20.merge[cols.num],as.numeric)

for( i in 4:length(APS_mothur_species20.merge))
{
  APS_mothur_species20.merge[,i]=round((APS_mothur_species20.merge[,i]/sum(APS_mothur_species20.merge[,i]))*100,digits = 1)
}
rm(i)

write.table( APS_mothur_species20.merge, file="species20_level.txt", sep = "\t",row.names = F,col.names = T, quote = F)
read.table("Phylum20_level.txt", header = T,stringsAsFactors = F, sep = "\t" )
APS_mothur_phylum20.merge=read.delim("Phylum20_level.txt",header = T,stringsAsFactors = F, sep = "\t")
names(APS_mothur_species20.merge)[c(4:ncol(APS_mothur_species20.merge))]<-c("GR-6","HR-6","KR-1","MR-6","HDR-6","KR-6","VC-1","VS-6")

APS_mothur_species20.merge$Genus=do.call(rbind,strsplit(APS_mothur_species20.merge$Genus,split = "__"))[,2]
APS_mothur_species20.merge$species=do.call(rbind,strsplit(APS_mothur_species20.merge$species,split = "__"))[,2]
APS_mothur_species20.merge=cbind(newcol=paste(APS_mothur_species20.merge$Genus,APS_mothur_species20.merge$species),APS_mothur_species20.merge)
names(APS_mothur_species20.merge)[1]<-"Species"
rownames(APS_mothur_species20.merge)=APS_mothur_species20.merge$Species
APS_mothur_species20.merge=APS_mothur_species20.merge[,c(-1:-4)]
names(APS_mothur_species20.merge)<-c("GR-6","HR-6","KR-1","MR-6","HDR-6","KR-6","VC-1","VS-6")
png("Top20_normalized_OTUs_species_heatmap.png", type="cairo", res = 200, units="in", width=6, height=4)
pheatmap(APS_mothur_species20.merge,show_rownames = T,show_colnames = T,cluster_cols = T,scale = "none",border_color = F,fontsize_col = 7,fontsize_row = 8,fontsize = 7,angle_col = 0)
dev.off()







##############################################
###table,barplot and heatmap Genus level#####
##############################################
#Table
APS_mothur_Genus=tax_glom(APS_mothur,taxrank = "Genus")
APS_mothur_Genus_Filtered = subset_taxa(APS_mothur_Genus,!grepl("*_unclassified",Genus))
Genus=cbind(APS_mothur_Genus_Filtered@tax_table,APS_mothur_Genus_Filtered@otu_table)
Genus=data.frame(Genus,stringsAsFactors = F)
Genus$OTU=rownames(Genus)
as.matrix(names(Genus))
Genus=Genus[,c(ncol(Genus),1:6,8:(ncol(Genus)-1))]
write.table(Genus,file="Genus_level.txt", sep = "\t",row.names = F,col.names = T, quote = F)

#Bar Plot

APS_mothur_Genus10 = data.frame(sort(taxa_sums(APS_mothur_Genus_Filtered), TRUE)[1:10])
APS_mothur_Genus10$OTU=rownames(APS_mothur_Genus10)
APS_mothur_Genus10=APS_mothur_Genus10[,c(2,1)]
colnames(APS_mothur_Genus10)=c("OTU","Size")
APS_mothur_Genus10.merge=merge(APS_mothur_Genus10,Genus,by = "OTU",all.x = T)
APS_mothur_Genus10.merge=APS_mothur_Genus10.merge[rev(order(APS_mothur_Genus10.merge$Size)),]
as.matrix(names(APS_mothur_Genus10.merge))
names(APS_mothur_Genus10.merge)[c(9:ncol(APS_mothur_Genus10.merge))]<-c("GR-6","HR-6","KR-1","MR-6","HDR-6","KR-6","VC-1","VS-6")
##names(APS_mothur_Genus10.merge)<-c("GR-6","HR-6","KR-1","MR-6","HDR-6","KR-6","VC-1","VS-6")
APS_mothur_Genus10.merge=APS_mothur_Genus10.merge[,c(8:ncol(APS_mothur_Genus10.merge))]
as.matrix(names(APS_mothur_Genus10.merge))
cols.num <- names(APS_mothur_Genus10.merge[,c(2:ncol(APS_mothur_Genus10.merge))])
write.table(APS_mothur_Genus10.merge, file="Genus20_level.txt", sep = "\t",row.names = F,col.names = T, quote = F)
APS_mothur_Genus10.merge=read.delim("Genus10_level.txt",header = T,stringsAsFactors = F, sep = "\t")
APS_mothur_Genus10.merge[cols.num] <- sapply(APS_mothur_Genus10.merge[cols.num],as.numeric)

for( i in 2:length(APS_mothur_Genus10.merge))
{
  APS_mothur_Genus10.merge[,i]=round((APS_mothur_Genus10.merge[,i]/sum(APS_mothur_Genus10.merge[,i]))*100,digits = 1)
}
rm(i)
APS_mothur_Genus10.merge_melt=melt(APS_mothur_Genus10.merge)
APS_mothur_Genus10.merge=read.delim("Genus10_level.txt",header = T,stringsAsFactors = F, sep = "\t")
names(APS_mothur_Genus10.merge)[c(2:ncol(APS_mothur_Genus10.merge))]<-c("GR-6","HR-6","KR-1","MR-6","HDR-6","KR-6","VC-1","VS-6")
colnames(APS_mothur_Genus10.merge_melt)=c("Genus","Sample","Percentage")
png("Top10_normalized_OTUs_Genus_Barplot.png", type="cairo", res = 200, units="in", width=10, height=6)
Genus_colors=rainbow(length(unique(APS_mothur_Genus10.merge$Genus)), s=.6, v=.9)[sample(1:length(unique(APS_mothur_Genus10.merge$Genus)),length(unique(APS_mothur_Genus10.merge$Genus)))]
ggplot(APS_mothur_Genus10.merge_melt,aes(x = Sample,y = Percentage,fill=Genus))+
  geom_bar(stat = "identity") + 
  scale_fill_manual(values = Genus_colors)+
  theme(axis.text.x = element_text(angle = 0, hjust = 1))
dev.off()

#Heatmap

APS_mothur_Genus20 = data.frame(sort(taxa_sums(APS_mothur_Genus_Filtered), TRUE)[1:20])
APS_mothur_Genus20$OTU=rownames(APS_mothur_Genus20)
APS_mothur_Genus20=APS_mothur_Genus20[,c(2,1)]
colnames(APS_mothur_Genus20)=c("OTU","Size")
APS_mothur_Genus20.merge=merge(APS_mothur_Genus20,Genus,by = "OTU",all.x = T)
APS_mothur_Genus20.merge=APS_mothur_Genus20.merge[rev(order(APS_mothur_Genus20.merge$Size)),]
as.matrix(names(APS_mothur_Genus20.merge))
names(APS_mothur_Genus20.merge)[c(9:ncol(APS_mothur_Genus20.merge))]<-c("GR-6","HR-6","KR-1","MR-6","HDR-6","KR-6","VC-1","VS-6")
APS_mothur_Genus20.merge=APS_mothur_Genus20.merge[,c(8:ncol(APS_mothur_Genus20.merge))]
as.matrix(names(APS_mothur_Genus20.merge))
write.table(APS_mothur_Genus20.merge, file="Genus20_level.txt", sep = "\t",row.names = F,col.names = T, quote = F)
read.table("Genus20_level.txt",header = T,stringsAsFactors = F, sep = "\t")
APS_mothur_Genus20.merge=read.delim("Genus20_level.txt",header = T,stringsAsFactors = F, sep = "\t")
cols.num <- names(APS_mothur_Genus20.merge[,c(2:ncol(APS_mothur_Genus20.merge))])
APS_mothur_Genus20.merge[cols.num] <- sapply(APS_mothur_Genus20.merge[cols.num],as.numeric)

for( i in 2:length(APS_mothur_Genus20.merge))
{
  APS_mothur_Genus20.merge[,i]=round((APS_mothur_Genus20.merge[,i]/sum(APS_mothur_Genus20.merge[,i]))*100,digits = 1)
}
rm(i)
rownames(APS_mothur_Genus20.merge)=APS_mothur_Genus20.merge$Genus
names(APS_mothur_Genus20.merge)[c(2:ncol(APS_mothur_Genus20.merge))]<-c("GR-6","HR-6","KR-1","MR-6","HDR-6","KR-6","VC-1","VS-6")
APS_mothur_Genus20.merge=APS_mothur_Genus20.merge[,-1]
png("Top20_normalized_OTUs_Genus_heatmap.png", type="cairo", res = 200, units="in", width=8, height=4)
# pheatmap(APS_mothur_Genus20.merge,show_rownames = T,show_colnames = T,cluster_cols = F,scale = "none",annotation_col = anno,bGenus_color = F,fontsize_col = 7,fontsize_row = 8,fontsize = 7,cellheight = 10)
pheatmap(APS_mothur_Genus20.merge,show_rownames = T,show_colnames = T,cluster_cols = T,scale = "none",border_color = F,fontsize_col = 7,fontsize_row = 8,fontsize = 7,angle_col = 0)
dev.off()



##################################################
###### contig length distribution histogram ######
##################################################

length_distribution=read.delim(file = "input.trim.contigs.good.unique.good.summary",header = T,sep = "\t",comment.char = "#",stringsAsFactors = F)
png("Contig_length_distribution_histogram.png", type="cairo", res = 200, units="in", width=6, height=5)
hist(length_distribution$nbases,breaks = 100,xlim = c(350,550),main="Contig length Distribution",xlab = "Contig Length (Bases)",ylab = "Frequency")
dev.off()

