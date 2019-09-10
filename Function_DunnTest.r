#!/usr/bin/env Rscript 
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
    make_option(c("-c", "--category"),metavar="string",dest="group", help="Specify category name in mapping file",default="none"),
    make_option(c("-p", "--prefix"),metavar="str", dest="prefix",help="The prefix of output files, default if ''. Diabled where output is infer",default=""),
    make_option(c("-o", "--output"),metavar="string",dest="out", help="Where to store the results",default="infer")
    # make_option(c("-l", "--clean"),metavar="logical",dest="clean", help="if T, clean output file name, make it compatable to more piplines",default=FALSE)
    )

opt <- parse_args(OptionParser(option_list=option_list,description = "R script for generating Dunn test output"))





#setwd("~/Desktop")
#map = "sample-metadata.PCA.txt"
#group = "Group1"
#KEGG_table= "feature-table.metagenome.L2.PCA.txt"
filename_temp<-opt$input
KEGG_table_temp<-gsub(pattern = "\\.PCA.txt$", "", filename_temp)
#output_file="KEGG.txt"
if(opt$out=="infer"){
    output_file=paste(KEGG_table_temp,"_",opt$group,"_DunnTest.txt",sep="")
}else{
    if(!dir.exists(opt$out)){dir.create(opt$out,recursive = T)}
    opt$out<-paste(opt$out,"/",opt$prefix,sep="")
    output_file=paste(opt$out, opt$group,"_DunnTest.txt",sep="")
}


#print(opt$input)
#print(opt$map)
#print(opt$group)
#print(opt$out)

map<-read.table(opt$map, header = T, row.names = 1, check.names = F,stringsAsFactors = F,sep = "\t",comment.char = "", na.strings="", quote="")
table<-read.table(opt$input, header = T, check.names = F,stringsAsFactors = F,sep = "\t",comment.char = "", quote = "")
#head(map)
table <-  table[!duplicated(table[,1]),]
rownames(table) <- table[,1]
table <- table[, -1]

map2<-map[order(rownames(map)), ]
#head(map2)

map3<-map2[opt$group]
colnames(map3)<-"Group"
#head(map3)
if(!is.numeric(table[,ncol(table)])){
    table1<-t(table[,-ncol(table)])
}else{
    table1<-t(table)
}
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
    tryCatch({post_hoc<-dunn.test(data[,i],g=data$Group, method="bonferroni");
        dunntest_results = cbind.data.frame(data_colname[i],kw_test_results$p.value, post_hoc$comparisons,post_hoc$Z,post_hoc$P.adjusted);
        output <- rbind(as.matrix(output), as.matrix(dunntest_results))}, 
        error=function(e){print("Case can not be dunn-tested")})
    #print(post_hoc)
    #print(data_colname[i]) 
}

write.table(as.matrix(output), output_file, quote=FALSE, col.names=NA, sep="\t")


