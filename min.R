ag<-commandArgs(T)
data<-read.table(ag[1],sep='\t',row.names=1,header=T)
min=floor(min(colSums(data[,-ncol(data)]))/1000)*1000
write(as.character(min), stdout())