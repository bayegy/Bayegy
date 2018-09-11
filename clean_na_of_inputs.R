library(optparse)

#######arguments
option_list <- list( 
    make_option(c("-m", "--map"),dest="m", help="Specify the path of mapping file with the Sample ID at the first column [Required]",default=NULL),
    make_option(c("-g", "--group"),dest="g", help="Specify group name according to which the map and input files will be cleaned [Required]",default=NULL),
    make_option(c("-t", "--table-qza"), dest="t",help="qza otu table needed to clean na.",default=NULL),
    make_option(c("-i", "--input"), dest="i",help="the otu table needed to clean na. Sample id must be the column name. All the columns not in sample ID will be saved",default=NULL),
    make_option(c("-d", "--distance-matrix"), dest="d",help="the .qza distance matrix needed to clean na.",default=NULL),
    make_option(c("-o", "--output"),dest="out", help="Specify the path of output files",default="./"),
    make_option(c("-a", "--na-string"),dest="a", help="Specify the character representing the missing values.Default is ''",default="")
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
write.table(dat,file = paste(opt$out,"/","cleaned_map.txt",sep = ""),row.names = F,col.names = T,quote = F,sep = "\t",append = F)

if(!is.null(opt$i)){
    otu <- read.table(opt$i,skip=1,comment.char="",check.names=F,stringsAsFactors=F, header = TRUE, sep = "\t",na.strings=opt$a)
    sel1<-sel[match(colnames(otu),rownames(group))]
    sel1<-(is.na(sel1)|sel1)
    otu<-otu[,sel1]
    write.table("# Constructed from biom file",file = paste(opt$out,"/","cleaned_feature_table.txt",sep = ""),row.names = F,col.names = F,quote = F,sep = "\t",append = F)
    write.table(otu,file = paste(opt$out,"/","cleaned_feature_table.txt",sep = ""),row.names = F,col.names = T,quote = F,sep = "\t",append = T)
}

if(!is.null(opt$d)){
    system(sprintf("qiime tools export --input-path %s --output-path %s",opt$d,opt$out))
    distance <- read.table(sprintf("%s/distance-matrix.tsv",opt$out),comment.char="",check.names=F,stringsAsFactors=F, header = TRUE, row.names=1,sep = "\t")
    selc<-sel[match(colnames(distance),rownames(group))]
    selr<-sel[match(rownames(distance),rownames(group))]
    distance<-distance[selr,selc]
    distance<-data.frame(s=rownames(distance),distance,check.names=F)
    colnames(distance)[1]<-""
    write.table(distance,file = paste(opt$out,"/","cleaned_distance.tsv",sep = ""),row.names = F,col.names = T,quote = F,sep = "\t",append = F)
    system(sprintf("qiime tools import --type DistanceMatrix --input-path %s/cleaned_distance.tsv  --output-path %s/cleaned_distance.qza",opt$out,opt$out))

}

if(!is.null(opt$t)){
    system(sprintf("qiime feature-table filter-samples --i-table %s --m-metadata-file %s/cleaned_map.txt --o-filtered-table %s/filtered_feature_table.qza",opt$t,opt$out,opt$out))
}
