##assuming you have a file containing gene expression named as gexprs.txt
##reading a text file 
getwd()
workdir=getwd()
##to set the workdir
setwd(workdir)

##read the delim file
countdata <- read.delim("C:/Users/Eshwari/Desktop/Bioinfo_scripts/GeneExpression.txt",header=T,sep = "\t",stringsAsFactors = F,comment.char = '#') 
head(countdata)
##library(DESeq2)
##retain all rows, extract only the first column and the columns containing the read counts
as.matrix(countdata)
##assigning rownames of the countdata
rownames(countdata)=countdata$Geneid
##check if the rownames are assigned 
rownames(countdata)
##
names(countdata)
##subsetting only the gene expression count matrix 
gene_count <- as.matrix[,c(1,7:ncol(countdata))]
gene_count <- gene_count[,c(1,7:ncol(countdata))]


##create data frame, assign columns rows and names, add delete columns rows, subset columns names
## read files in R
## lappaly, sapply, tidyr functions 
## bioconductor pacakges 
## plot heatmaps, PCA plots, volcano plots
## perform normalization in R 


read.delim(filename="geneexpression.txt", sep = "\t", header = T, stringsAsFactors = F)
read.csv
read.table 

## to view data 
head(data)
summary(data)
hist(data$columnname)

##data manipulation
filtered_data <- as.matrix(countdata[,c(1,7:ncol(countdata))])
##lapply, apply(), sapply()
""" this family of functions help us to apply certain function to a data frame, list, or vector, and return result as list or vector depending on the function we use
##1. apply()
apply(inputmatrix, 1, sum) ##1 is for across rows
apply(inputmatrix, 2, mean) ##2 is for across columns

##2. lapply()
##applies functions to list objects and returns list object of the same length
input = list, vector, data frame
output = list

#3. sapply()
input = list, vector, data frame
output = array or matrix


##tidyr - this package is used to tidy up the rows and columns in R
##dplyr, tidyr, ggplot2
%>% - pipe operator 


##rbind and cbind

##merging

##printing output to a file

##printing a report

##creating pca plots, heatmaps, volcano plots





