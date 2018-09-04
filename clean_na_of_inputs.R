library(optparse)

#######arguments
option_list <- list( 
    make_option(c("-i", "--input"), dest="i",help="the files needed to clean na. Sample id must be the column name but not first and last column",default=NA),
    make_option(c("-m", "--map"),dest="m", help="Specify the path of mapping file",default=NULL),
    make_option(c("-g", "--group"),dest="g", help="Specify group name according to which the map and input files will be cleaned",default="none"),
    make_option(c("-o", "--output"),dest="out", help="Specify the path of output files",default="./")
    make_option(c("-a", "--na-string"),dest="a", help="Specify the character represent the missing values",default="NA")
    )

opt <- parse_args(OptionParser(option_list=option_list))

####check outdir
if(!dir.exists(opt$out)){dir.create(opt$out,recursive = T)}


####clean procedure
dat <- read.table(opt$m,comment.char="",check.names=F,stringsAsFactors=F, header = TRUE, sep = "\t",na.strings=opt$a)
otu <- read.table(opt$i,comment.char="",check.names=F,stringsAsFactors=F, header = TRUE, sep = "\t",na.strings=opt$a)
group<-dat[opt$g]
sel<-(!is.na(group))
dat<-dat[sel,]

id<-colnames(otu)[-c(1,length(otu))]
sel1<-c(T,sel[match(id,rownames(group))],T)
otu<-otu[,sel1]

write.table(dat,file = paste(opt$out,"/","cleaned_map.txt",sep = ""),row.names = F,col.names = T,quote = F,sep = "\t",append = F)
write.table(otu,file = paste(opt$out,"/","cleaned_table.txt",sep = ""),row.names = F,col.names = T,quote = F,sep = "\t",append = F)