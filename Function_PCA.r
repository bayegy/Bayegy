library("ggbiplot")
library("stringr")
library(optparse)
require(reshape)
require(ggplot2)
require(ggpubr)
library(dplyr)
library(devtools)
#install_github("vqv/ggbiplot", ref = "experimental") 

option_list <- list( 
    make_option(c("-i", "--input"),metavar="path", dest="input",help="Specify the path of the PICRUSt output file",default=NULL),
    make_option(c("-m", "--map"),metavar="path",dest="map", help="Specify the path of mapping file",default=NULL),
    make_option(c("-g", "--category"),metavar="string",dest="group", help="Specify category name in mapping file",default="none")
    )

opt <- parse_args(OptionParser(option_list=option_list,description = "R script for generating Dunn test output"))

args <- commandArgs(trailingOnly = TRUE)
KEGG_function_txt = opt$input
meta_txt = opt$map
category_1= opt$group

this.dir <- getwd()
setwd(this.dir)

#setwd("~/Desktop")
#KEGG_function_txt = "feature-table.metagenome.L2.PCA.txt"
#meta_txt = "sample-metadata.PCA.txt"
#category_1= "Group1"

design = read.table(meta_txt, header=T, row.names= 1, sep="\t",check.names = F, na.strings = "",  fill = TRUE) 
#head(design)

#design<-design[!is.na(design[category_1]),]


table = read.table(KEGG_function_txt, row.names= 1,  header=T, sep="\t",check.names=F, na.strings = "",  fill = TRUE)
#head(table)

#remove the last column and transpose
table2<-t(table[,-(ncol(table))])
#summary(table2)

#Sort the table
table3<-table2[order(rownames(table2)), ]

#Sort the design file
design2<-design[order(rownames(design)), ]
#dim(table3)
#dim(design2)

#Find the common samples in both design and table
idx = rownames(design2) %in% rownames(table3)
design3 = design2[idx,]

design4<-design3[category_1]
colnames(design4)<-"Group"
#head(map3)

design5<-design3[!is.na(design4["Group"]),]

#may not necessary for table4, as the number should match already but just put here in case.
table4 = table3[rownames(design4), ]
#dim(table4)

#Convert to proportion
table5<-prop.table(as.matrix(table4), margin = 1)

#Remove the ones with zero variance
table6<-table5[ , apply(table5, 2, var) != 0]
rowSums(table6,na=T)
dim(table6)

joinedtab<-data.frame(design4,table6,check.rows = T,check.names = F)
data<-joinedtab[!is.na(joinedtab["Group"]),]
table6_2<-data[,-1]

#Transpose
table7<-t(table6_2)

# 筛选mad值：按mad值排序取前10波动最大的OTUs
table8 = head(table7[order(apply(table7,1,mad), decreasing=T),],n=10)

#Transpose back
table9<-t(table8)

pc <- prcomp(table9, center=TRUE, scale. = TRUE)

print(paste("Making PCA plots for", KEGG_function_txt, sep=" "))
KEGG_function_txt<-str_replace(KEGG_function_txt,"PCA.txt","")
###########The name is bad here use tools::file_path_sans_ext("ABCD.csv") to obtain the basename in the future.
PCA_plot_outputpdfname <- paste(KEGG_function_txt, category_1,".PCA.pdf", sep="")
p<-ggbiplot(pc, obs.scale = 1, var.scale = 1, groups = design5[[category_1]], ellipse = TRUE, varname.adjust = 1.2, varname.abbrev=T, varname.size=3) +
  scale_color_discrete(name = '') +
  theme_bw()
  

ggsave(PCA_plot_outputpdfname, plot=p, height=6, width=8, dpi = 300)
