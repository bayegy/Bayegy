library('optparse')
option_list <- list(
	make_option(c("-i", "--input"),metavar="path", dest="data",help="Specify the abundance table you want to debatch, you should not contain any comment in the file, taxonomy or otuid must in the firsr as well as last column",default=NULL),
	make_option(c('-f',"--feature-table-taxonomy"),metavar="path", dest="otu",help="Specify the otu table you want to debatch with comment at first row, if --input is supplied, this one will be ignored.",default=NULL),
    make_option(c("-m", "--map"),metavar="path",dest="map", help="Specify the path of mapping file",default=NULL),
    make_option(c("-c", "--category"),metavar="string",dest="group", help="Specify category name in mapping file, the effect of which will be eliminated. If not passed, no batch effect will be corrected",default=NULL),
    make_option(c("--min-samples"),metavar="int or float",dest="sample", help="Pass this to filter your otu table before debatch",default=NULL),
    make_option(c("--min-frequency"),metavar="int or float",dest="frequency", help="Pass this to filter your otu table before debatch",default=NULL),
    make_option(c("--only-mean"),metavar='logical',dest="onlymean", help="If TRUE, only mean will be corrected, otherwise both variance and mean will be corrected",default=TRUE),
	make_option(c("-o", "--output"),metavar="directory",dest="out", help="Specify the directory where the debatched abundance file will be placed",default="./")  
	)

opt <- parse_args(OptionParser(option_list=option_list))

library(stringr)
library(sva)

if(!is.null(opt$data)){
  data<- read.table(opt$data,comment.char="",check.names=F,stringsAsFactors=F, header = TRUE, sep = "\t",na.strings='')
}else{
  data<- read.table(opt$otu,skip = 1,comment.char="",check.names=F,stringsAsFactors=F, header = TRUE, sep = "\t",na.strings='')
}


backdata<-data[,c(1,length(data))]
colname_backdata<-colnames(backdata)
data<-data[,-c(1,length(data))]

if(!is.null(opt$sample)){
  ifelse(as.numeric(opt$sample)>=1,sel<-apply(data, 1, function(x){return(sum(x>0)>=as.numeric(opt$sample))}),
         sel<-apply(data, 1, function(x){return(sum(x>0)>=(as.numeric(opt$sample)*ncol(data)))}))
  data<-data[sel,]
  backdata<-backdata[sel,]
}

if(!is.null(opt$frequency)){
  sel<-(colSums(t(data))>=as.numeric(opt$frequency))
  data<-data[sel,]
  backdata<-backdata[sel,]
}
if (!is.null(opt$group)){
  map<-read.table(opt$map,row.names = 1,comment.char="",check.names=F,stringsAsFactors=F, header = TRUE, sep = "\t",na.strings='')
  map<-map[match(colnames(data),rownames(map)),]
  batch <- as.factor(map[opt$group][,1])
  modcombat = model.matrix(~1, data=map)
  data <- ComBat(dat=data, batch=batch, mod=modcombat, par.prior=F, prior.plots=F,mean.only = as.logical(opt$onlymean))
}

data<-data.frame(backdata[,1],data,backdata[,2],check.names = F)
colnames(data)[c(1,length(data))]<-colname_backdata

if(!is.null(opt$data)){
  write.table(data,opt$out,sep = '\t',quote = F,row.names = F)
}else{
  cm<-"#Feature table filterd or debatched"
  write.table(cm,opt$out,sep = '\t',quote = F,row.names = F,col.names=F,append = F)
  write.table(data,opt$out,sep = '\t',quote = F,row.names = F,append = T)
}



