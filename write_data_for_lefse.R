library(optparse)
library(stringr)

option_list <- list( 
    make_option(c("-i", "--input"),metavar="file", dest="otu",help="Specify the path of collapsed bacteria taxonomy file.The first row must not contain any comments",default=NULL),
    make_option(c("-m", "--map"),metavar="file",dest="map", help="Specify the path of mapping file",default=NULL),
    make_option(c("-c", "--category"),metavar="string",dest="group", help="Specify category name in mapping file. You may specify more than one category seprated by commas.",default="none"),
    make_option(c("-o", "--output"),metavar="path",dest="out", help="Specify the path of output files",default="./")
    )

opt <- parse_args(OptionParser(option_list=option_list,description = "This script is used to write the input file of Lefse."))
if(!dir.exists(opt$out)){dir.create(opt$out,recursive = T)}
#ag<-commandArgs(T)
#if(length(ag)<3){
#  print("please specify: 
#1.collapsed bacteria taxonomy file
#2.mapping file
#3.column name in mapping file of group seprated by ','
#4.full path of out file")
#}else{
meta<-read.table(opt$map,na.strings="",row.names=1,header = T,sep = "\t",comment.char = "",check.names = F,stringsAsFactors = F)
group<-str_split(opt$group,",")[[1]]
if(length(group)==1){
  meta<-na.omit(meta[group])
}else{
  meta<-na.omit(meta[,colnames(meta)%in%group])
}
meta<-data.frame(Subject=rownames(meta),meta)
data<-read.table(opt$otu,header = T,sep = "\t",comment.char = "",stringsAsFactors = F,check.names = F)
data<-data[data[,1]!="Others"&data[,1]!="unclassified",]
data<-data[,c(length(data),2:(length(data)-1))]
data[,1]<-str_replace(data[,1],";$","")
data[,1]<-str_replace_all(data[,1],";","|")
#data<-data[-dim(data)[1],]

meta<-t(meta)
data<-data[,c(1,match(meta[1,],colnames(data)))]
write.table(meta,file = opt$out,row.names = T,col.names = F,quote = F,sep = "\t",append = F)
write.table(data,file = opt$out,row.names = F,col.names = F,quote = F,sep = "\t",append = T)

#}

