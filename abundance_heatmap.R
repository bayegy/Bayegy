library(optparse)

#######arguments
option_list <- list( 
    make_option(c("-i", "--input"),metavar="path", dest="otu",help="Specify the path of collapsed bacteria table",default=NULL),
    make_option(c("-n", "--number"),metavar="int", dest="num",help="The number of most related species you want to plot, default is 20",default=15),
    make_option(c("-a", "--min-abundance"),metavar="float", dest="mina",help="Min correlation coefficent to label significance, default is 0",default=NULL),
    make_option(c("-p", "--prefix"),metavar="str", dest="prefix",help="The prefix of output files, default if null",default=""),
    make_option(c("-o", "--output"),metavar="directory",dest="out", help="Specify the directory of output files",default="./")
    )

opt <- parse_args(OptionParser(option_list=option_list,description = "This script is used to relate the bacteria species and enviromental factors(numeric), and a heatmap is used to visualize the rank correlation coefficent"))

library(pheatmap)

if(!dir.exists(opt$out)){dir.create(opt$out,recursive = T)}
opt$out<-paste(opt$out,"/",opt$prefix,sep="")

otu <- read.table(opt$otu,quote="",comment.char="",check.names=F,stringsAsFactors=F, header = TRUE, sep = "\t")

otu<-otu[otu[,1]!="Others"&otu[,1]!="unclassified",]

#otu<-otu[!duplicated(otu[,1]),]
rownames(otu)<-otu[,ncol(otu)]

otu<-t(otu[,-c(1,ncol(otu))])
sum_otu<-colSums(otu)

if(!is.null(opt$mina)){
	otu<-otu[,sum_otu>=as.numeric(opt$mina)]
}else{
	sel<-head(order(sum_otu,decreasing = T),as.numeric(opt$num))
	otu<-otu[,sel]
}

otu<-log(otu+1,base=10)
#otu<-log((otu+min(otu[otu!=0]))*10000)
#otu<-scale(otu,center=T,scale=T)

#apply(otu,2,mean)
#apply(otu,2,sd)

#print(otu)
pdf(paste(opt$out,"abundance_heatmap.pdf",sep = ""), height=7+0.5*dim(otu)[1],width=3+0.5*dim(otu)[2])
pheatmap(otu,fontsize=15,border_color = "black",
         fontsize_row =15,fontsize_col = 15,
         cluster_rows=T,clustering_distance_rows="euclidean",
         cluster_cols=T,clustering_distance_cols="euclidean")
dev.off()
