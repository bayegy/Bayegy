ag<-commandArgs(T)
data<-read.table(ag[1],sep=',',header=F)
write(as.character(min(data[,2])), stdout())