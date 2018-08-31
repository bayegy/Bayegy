
library(stringr)
ag<-commandArgs(T)

if(length(ag)<3){
  print("please specify: 
1.collapsed bacteria taxonomy file
2.mapping file
3.column name in mapping file of group seprated by ','
4.full path of out file")
}else{
  meta<-read.table(ag[2],row.names=1,header = T,sep = "\t",comment.char = "",check.names = F,stringsAsFactors = F)
  group<-str_split(ag[3],",")[[1]]
  if(length(group)==1){
    meta<-meta[group]
  }else{
    meta<-meta[,colnames(meta)%in%group]
  }
  meta<-data.frame(Subject=rownames(meta),meta)
  data<-read.table(ag[1],header = T,sep = "\t",comment.char = "",stringsAsFactors = F,check.names = F)
  data<-data[data[,1]!="Others",]
  data<-data[,c(length(data),2:(length(data)-1))]
  data[,1]<-str_replace(data[,1],";$","")
  data[,1]<-str_replace_all(data[,1],";","|")
  #data<-data[-dim(data)[1],]
  
  meta<-t(meta)
  meta<-meta[,match(colnames(data)[-1],meta[1,])]
  write.table(meta,file = ag[4],row.names = T,col.names = F,quote = F,sep = "\t",append = F)
  write.table(data,file = ag[4],row.names = F,col.names = F,quote = F,sep = "\t",append = T)
  
}

