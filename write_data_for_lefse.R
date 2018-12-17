#library(optparse)
library(stringr)

#option_list <- list( 
#    make_option(c("-i", "--input"),metavar="path", dest="otu",help="Specify the path of collapsed bacteria taxonomy file.The first row must not contain any comments",default=NULL),
#    make_option(c("-m", "--map"),metavar="path",dest="map", help="Specify the path of mapping file",default=NULL),
#    make_option(c("-c", "--category"),metavar="string",dest="group", help="Specify category name in mapping file. You may specify more than one category seprated by commas.",default="none"),
#    make_option(c("-o", "--output"),metavar="path",dest="out", help="Specify the path of output file",default="./")
#    )

#opt <- parse_args(OptionParser(option_list=option_list,description = "This script is used to write the input file of Lefse."))
#if(!dir.exists(opt$out)){dir.create(opt$out,recursive = T)}
ag<-commandArgs(T)
if(length(ag)<5){
  write("please specify: 
1.Input file, with taxonomy at last column, id or short name at first column

2.Mapping file

3.Categories name seprated by ',' in mapping file

4.Path of out file

5.Use 'T' to Skip the first line(e.g. comment line) of input file, 'F' not

Sample Usage: 
Rscript write_data_for_lefse.R  otu_table_with_comment.txt  mapping_file.txt  Group1  Group1_table_for_lefse.txt T", stdout())
}else{
meta<-read.table(ag[2],na.strings="",row.names=1,header = T,sep = "\t",comment.char = "",check.names = F,stringsAsFactors = F)
group<-str_split(ag[3],",")[[1]]

meta<-na.omit(meta[group])

meta<-data.frame(Subject=rownames(meta),meta)
if(as.logical(ag[5])){
	data<-read.table(ag[1],skip=1,header = T,sep = "\t",comment.char = "",stringsAsFactors = F,check.names = F)
}else{
	data<-read.table(ag[1],header = T,sep = "\t",comment.char = "",stringsAsFactors = F,check.names = F)
}
#calculate the relative abundance
#data[,2:(ncol(data)-1)]=apply(data[,2:(ncol(data)-1)],2,function(x){x/sum(x)})


data<-data[data[,1]!="Others"&data[,1]!="unclassified",]
data<-data[,c(ncol(data),2:(ncol(data)-1))]
data[,1]<-str_replace(data[,1],"; *[a-z]__ *;.*$","")
data[,1]<-str_replace(data[,1],"; *[a-z]__ *$","")
data[,1]<-str_replace(data[,1],";$","")
data[,1]<-str_replace_all(data[,1],";","|")

meta<-t(meta)
data<-data[,c(1,match(meta[1,],colnames(data)))]

#Filter otus
#otu_sum<-colSums(t(data[,-1])>0)
#data<-data[otu_sum>0.25*(ncol(data)-1),]

write.table(meta,file = ag[4],row.names = T,col.names = F,quote = F,sep = "\t",append = F)
write.table(data,file = ag[4],row.names = F,col.names = F,quote = F,sep = "\t",append = T)

}

