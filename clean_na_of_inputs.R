library(optparse)

#######arguments
option_list <- list( 
    make_option(c("-i", "--input"), dest="i",help="the files needed to clean na. Sample id must be the column name. All the columns not in sample ID will be saved",default=NULL),
    make_option(c("-m", "--map"),dest="m", help="Specify the path of mapping file with the Sample ID at the first column",default=NULL),
    make_option(c("-g", "--group"),dest="g", help="Specify group name according to which the map and input files will be cleaned",default="none"),
    make_option(c("-o", "--output"),dest="out", help="Specify the path of output files",default="./"),
    make_option(c("-a", "--na-string"),dest="a", help="Specify the character representing the missing values",default="NA")
    )

opt <- parse_args(OptionParser(option_list=option_list))

####check outdir
if(!dir.exists(opt$out)){dir.create(opt$out,recursive = T)}


####clean procedure
dat <- read.table(opt$m,comment.char="",check.names=F,stringsAsFactors=F, header = TRUE, sep = "\t",na.strings=opt$a)

group<-dat[opt$g]
rownames(group)<-dat[,1]
sel<-(!is.na(group))
dat<-dat[sel,]
if(!is.null(opt$i)){
    otu <- read.table(opt$i,comment.char="",check.names=F,stringsAsFactors=F, header = TRUE, sep = "\t",na.strings=opt$a)
    sel1<-sel[match(colnames(otu),rownames(group))]
    sel1<-(is.na(sel1)|sel1)
    otu<-otu[,sel1]
    write.table(otu,file = paste(opt$out,"/","cleaned_table.txt",sep = ""),row.names = F,col.names = T,quote = F,sep = "\t",append = F)
}

write.table(dat,file = paste(opt$out,"/","cleaned_map.txt",sep = ""),row.names = F,col.names = T,quote = F,sep = "\t",append = F)
