ag<-commandArgs(T)
data<-read.table(ag[1],sep='\t',row.names=1,header=T)
# 20201102,tangmaomao modify the min code
min=floor(min(colSums(data[,-ncol(data)])))
write(as.character(min), stdout())
