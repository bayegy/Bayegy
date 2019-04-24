library(optparse)

#######arguments
option_list <- list(
    make_option(c("-i", "--input"),metavar="path", dest="otu",help="Specify the path of collapsed bacteria table. Required",default=NULL),
    make_option(c("-m", "--map"),metavar="path",dest="map", help="Specify the path of mapping file. Optional",default=NULL),
    make_option(c("-c", "--category"),metavar="string",dest="group", help="Category to compare. Optional",default=NULL),
    make_option(c("-n", "--number"),metavar="int", dest="num",help="The number of most abundant species needed to be plotted, default is 20",default=20),
    make_option(c("-a", "--min-abundance"),metavar="float", dest="mina",help="The min abundance of species to be plotted, pass this will cause --number disabled",default=NULL),
    make_option(c("-p", "--prefix"),metavar="str", dest="prefix",help="The prefix of output files, default if ''",default=""),
    make_option(c("-b", "--by-groupMean"),metavar="logical", dest="bym",help="if T, to use the group mean to plot heatmap",default=FALSE),
    make_option(c("-l", "--last-column"),metavar="logical", dest="last",help="T(use the last column as the feature names) or F(use the first column as the feature names). The last column does not has to be feature name, could also be abundance.",default=TRUE),
    make_option(c("-t", "--transpose"),metavar="logical", dest="trans",help="if T, transpose the heatmap",default=TRUE),
    make_option(c("-u", "--cluster"),metavar="logical", dest="cluster",help="T(cluster the samples) or F",default=FALSE),
    make_option(c("-o", "--output"),metavar="directory",dest="out", help="Specify the directory of output files",default="./")
    )

opt <- parse_args(OptionParser(option_list=option_list,description = "This script is used to relate the bacteria species and enviromental factors(numeric), and a heatmap is used to visualize the rank correlation coefficent"))

library(pheatmap)
library(getopt)

if(!dir.exists(opt$out)){dir.create(opt$out,recursive = T)}
opt$out<-paste(opt$out,"/",opt$prefix,sep="")



no_group<-(is.null(opt$map)|is.null(opt$group))
if(!no_group){
    map<-read.table(opt$map,quote="",row.names = 1,na.strings = "",comment.char="",check.names=F,stringsAsFactors=F, header = TRUE, sep = "\t")
    group<-na.omit(map[c(opt$group,"Description")])
    group<-group[order(rownames(group)),]
    group<-group[order(group[,1]),]
}


otu <- read.table(opt$otu,quote="",comment.char="",check.names=F,stringsAsFactors=F, header = TRUE, sep = "\t")

otu<-otu[otu[,1]!="Others"&otu[,1]!="unclassified",]


if(opt$last){
    otu<-otu[!duplicated(otu[,ncol(otu)]),]
    rownames(otu)<-otu[,ncol(otu)]
}else{
    otu<-otu[!duplicated(otu[,1]),]
    rownames(otu)<-otu[,1]
}



p1<-max(nchar(rownames(otu)))

if(!is.numeric(otu[,ncol(otu)])){
  otu<-otu[,-ncol(otu)]
}

otu<-t(otu[,-1])
sum_otu<-colSums(otu)

if(!is.null(opt$mina)){
    otu<-otu[,sum_otu>=as.numeric(opt$mina)]
}else{
    sel<-head(order(sum_otu,decreasing = T),as.numeric(opt$num))
    otu<-otu[,sel]
}

otu<-log(otu+1,base=10)
p2<-3+(0.3*dim(otu)[1])+(0.05*p1)
p3<-dim(otu)[2]

if(!no_group){
    otu<-otu[match(rownames(group),rownames(otu)),]
}


if(opt$trans){
    otu<-t(otu)
}


ht<-ifelse(opt$trans,3+0.4*p3 ,ifelse(p2<50,p2,49.9))
wd<-ifelse(opt$trans, ifelse(p2<50,p2,49.9),3+0.4*p3)
cc<-ifelse(opt$trans,F,T)
cr<-ifelse(opt$trans,T,F)

if(opt$cluster){
    cc=TRUE
    cr=TRUE
}

annotation_row=NA
annotation_col=NA

if(!no_group&!opt$bym){
    if(opt$trans){
        annotation_row = NA
        annotation_col = group[opt$group]
    }else{
        annotation_row = group[opt$group]
        annotation_col = NA
    }
}


if(opt$bym){
    if(opt$trans){
        otu<-t(apply(otu,1,function(x){tapply(x,INDEX = group[,1],mean)}))
        p2<-3+(0.3*dim(otu)[2])+(0.05*p1)
        p3<-dim(otu)[1]
        wd<-ifelse(p2<50,p2,49.9)
        ht<-2+0.4*p3
    }else{
        otu<-apply(otu,2,function(x){tapply(x,INDEX = group[,1],mean)})
        p2<-3+(0.3*dim(otu)[1])+(0.05*p1)
        p3<-dim(otu)[2]
        ht<-ifelse(p2<50,p2,49.9)
        wd<-2+0.4*p3
    }
}

if(!no_group){
    base_dir<-normalizePath(dirname(get_Rscript_filename()))
    source(paste(base_dir,"/piputils/get_colors.R", sep = ""))
    groups_color<-get_colors(opt$group, opt$map)
    names(groups_color)<-sort(unique(group[,1]))
    rexp=sprintf("annotation_color=list(%s=groups_color)",opt$group)
    eval(parse(text = rexp))
}else{
    annotation_color=NA
}




write.table(t(otu),paste(opt$out,"table.txt",sep = ""),sep = "\t",quote=FALSE,col.names=NA)

pdf(paste(opt$out,"heatmap.pdf",sep = ""), height=ht,width=wd)
pheatmap(otu,annotation_row=annotation_row,
         annotation_col=annotation_col,fontsize=10,border_color = "black",
         color = colorRampPalette(colors = c("#FFCCCC","red","black"))(100),
         cluster_cols=cc,clustering_distance_cols="euclidean",
         cluster_rows=cr,clustering_distance_rows="euclidean",annotation_colors=annotation_color)
dev.off()

