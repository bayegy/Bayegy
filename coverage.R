library(optparse)

#######arguments
option_list <- list( 
    make_option(c("-i", "--input"),metavar="path", dest="otu",help="Specify the path of collapsed bacteria table",default=NULL),
    make_option(c("-a", "--append"),metavar="path", dest="append",help="alpha table need be merged",default=NULL),
    make_option(c("-p", "--prefix"),metavar="str", dest="prefix",help="The prefix of output files, default if null",default=""),
    make_option(c("-o", "--output"),metavar="directory",dest="out", help="Specify the directory of output files",default="./")
    )

opt <- parse_args(OptionParser(option_list=option_list,description = "This script is used to relate the bacteria species and enviromental factors(numeric), and a heatmap is used to visualize the rank correlation coefficent"))

library(pheatmap)

if(!dir.exists(opt$out)){dir.create(opt$out,recursive = T)}


#system(sprintf("qiime tools export --input-path %s --output-path %s&& biom convert --to-tsv -i %s/feature-table.biom -o %s/alpha_sampling_table.tsv",opt$rare,opt$out,opt$out,opt$out))


opt$out<-paste(opt$out,"/",opt$prefix,sep="")

otu <- read.table(opt$otu, row.names=1, skip=1, quote="",comment.char="",check.names=F,stringsAsFactors=F, header = TRUE, sep = "\t")

otu<-otu[,-ncol(otu)]

coverage<-apply(otu,2,function(x){1-(sum(x==1)/sum(x))})

coverage<-data.frame(coverage)

rownames(coverage)<-colnames(otu)


if(!is.null(opt$append)){
	alpha <- read.table(opt$append, row.names=1, quote="",comment.char="",check.names=F,stringsAsFactors=F, header = TRUE, sep = "\t")
	coverage<- data.frame(alpha,coverage=coverage[match(rownames(alpha),rownames(coverage)),],check.names = F,check.rows = T)
}


write.table(coverage, paste(opt$out, "alpha-summary.tsv", sep="") , quote=FALSE, col.names=NA, sep="\t")