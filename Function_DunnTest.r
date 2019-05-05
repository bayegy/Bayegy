#R script for Dunn test
library(optparse)
require(reshape)
require(ggplot2)
require(ggpubr)
library(dplyr)
library(ggsignif)
library(dunn.test)

option_list <- list( 
    make_option(c("-i", "--input"),metavar="path", dest="input",help="Specify the path of the PICRUSt output file",default=NULL),
    make_option(c("-m", "--map"),metavar="path",dest="map", help="Specify the path of mapping file",default=NULL),
    make_option(c("-g", "--category"),metavar="string",dest="group", help="Specify category name in mapping file",default="none")
    )

opt <- parse_args(OptionParser(option_list=option_list,description = "R script for generating Dunn test output"))

#setwd("~/Desktop")
#map = "sample-metadata.PCA.txt"
#group = "Group1"
#KEGG_table= "feature-table.metagenome.L2.PCA.txt"
filename_temp<-opt$input
KEGG_table_temp<-gsub(pattern = "\\.PCA.txt$", "", filename_temp)
#output_file="KEGG.txt"
output_file=paste(KEGG_table_temp,"_",opt$group,"_DunnTest.txt",sep="")

#print(opt$input)
#print(opt$map)
#print(opt$group)
#print(opt$out)

map<-read.table(opt$map, header = T, row.names = 1, check.names = F,stringsAsFactors = F,sep = "\t",comment.char = "", na.strings="")
table<-read.table(opt$input, header = T, row.names = 1, check.names = F,stringsAsFactors = F,sep = "\t",comment.char = "")
#head(map)

map2<-map[order(rownames(map)), ]
#head(map2)

map3<-map2[opt$group]
colnames(map3)<-"Group"
#head(map3)

table1<-t(table[,-ncol(table)])
#head(table1)

table2<-table1[order(rownames(table1)), ]
#head(table2)

table3<-prop.table(as.matrix(table2), margin = 1)
table4<-table3[ , apply(table3, 2, var) != 0]
rowSums(table4,na=T)

joinedtab<-data.frame(map3,table4,check.rows = T,check.names = F)

data<-joinedtab[!is.na(joinedtab["Group"]),]

data_colname<-colnames(data)
output<-data.frame("KEGG_pathway"=character(),"KW_pvalue"=character(),"DunnTest_comparison"=character(),"DunnTest_Z"=character(),"DunnTest_PValueAdjusted"=character())
#colnames(out)<-c("KEGG_pathway","KW_pvalue","DunnTest_comparison","DunnTest_Z","DunnTest_PValueAdjusted")

for (i in 2:ncol(data)){
    data$Group <- as.factor(data$Group)
    kw_test_results<-kruskal.test(data[,i]~data$Group)
    pvals<-kw_test_results$p.value
    post_hoc<-dunn.test(data[,i],g=data$Group, method="bonferroni")
    #print(data_colname[i])
    #print(post_hoc)
    dunntest_results = cbind.data.frame(data_colname[i],kw_test_results$p.value, post_hoc$comparisons,post_hoc$Z,post_hoc$P.adjusted)
    output <- rbind(as.matrix(output), as.matrix(dunntest_results))
}

write.table(as.matrix(output), output_file, quote=FALSE, col.names=NA, sep="\t")


