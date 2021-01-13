library(optparse)

#######arguments
option_list <- list(
    make_option(c("-i", "--input"),metavar="path", dest="otu",help="Specify the path of collapsed bacteria table. Required",default=NULL),
    make_option(c("-m", "--map"),metavar="path",dest="map", help="Specify the path of mapping file. Optional",default=NULL),
    make_option(c("-c", "--category"),metavar="string",dest="group", help="Category to compare. Optional",default=NULL),
    make_option(c("-O", "--order"),metavar="string",dest="order", help="Sequence to reorder samples",default=NULL),
    make_option(c("-C", "--cellcolors"),metavar="string",dest="cellcolors", help="Comma seprated 3 colors",default="#FFCCCC,#FF0000,#222222"),
    make_option(c("-j", "--colors"),metavar="string",dest="colors", help="Comma seprated group colors.",default=NULL),
    make_option(c("-n", "--number"),metavar="int", dest="num",help="The number of most abundant species needed to be plotted, default is 20",default=20),
    make_option(c("-a", "--min-abundance"),metavar="float", dest="mina",help="The min abundance of species to be plotted, pass this will cause --number disabled",default=NULL),
    make_option(c("-p", "--prefix"),metavar="str", dest="prefix",help="The prefix of output files, default if ''",default=""),
    make_option(c("-b", "--by-groupMean"),metavar="logical", dest="bym",help="if T, to use the group mean to plot heatmap",default=FALSE),
    make_option(c("-l", "--last-column"),metavar="logical", dest="last",help="T(use the last column as the feature names) or F(use the first column as the feature names). The last column does not has to be feature name, could also be abundance.",default=TRUE),
    make_option(c("-t", "--transpose"),metavar="logical", dest="trans",help="if T, transpose the heatmap",default=TRUE),
    make_option(c("-M", "--method"),metavar="string", dest="method",help="method used to transform the values; log, zscore or none",default="log"),
    make_option(c("-u", "--cluster"),metavar="logical", dest="cluster",help="T(cluster the samples) or F",default=FALSE),
    make_option(c("-W", "--cellwidth"),metavar="int", dest="cellwidth",help="cell width in points",default=20),
    make_option(c("-H", "--cellheight"),metavar="int", dest="cellheight",help="cell height in points",default=20),
    make_option(c("-S", "--fontsize"),metavar="int", dest="fontsize",help="font size in pixels",default=10),
    make_option(c("-o", "--output"),metavar="directory",dest="out", help="Specify the directory of output files",default="./")
    )

opt <- parse_args(OptionParser(option_list=option_list,description = "This script is used to relate the bacteria species and enviromental factors(numeric), and a heatmap is used to visualize the rank correlation coefficent"))

library(pheatmap)
library(getopt)
library(stringr)

if(!dir.exists(opt$out)){dir.create(opt$out,recursive = T)}
opt$out<-paste(opt$out,"/",opt$prefix,sep="")



no_group<-(is.null(opt$map)|is.null(opt$group))

if(!is.null(opt$map)){
    map<-read.table(opt$map,quote="",row.names = 1,na.strings = "",comment.char="",check.names=F,stringsAsFactors=F, header = TRUE, sep = "\t")
}

if(!no_group){
    group<-map[opt$group]
    group$Description=0
    group<-na.omit(group)
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

if(!is.numeric(otu[,ncol(otu)])){
  otu<-otu[,-ncol(otu)]
}

otu<-t(otu[, -1])
sum_otu<-colSums(otu)

if(!is.null(opt$mina)){
    otu<-otu[,sum_otu>=as.numeric(opt$mina)]
}else{
    sel<-head(order(sum_otu,decreasing = T),as.numeric(opt$num))
    otu<-otu[,sel]
}


if(opt$method=="log"){
    otu<-log(otu+1,base=10)
}else if(opt$method=="zscore"){
    otu<-scale(otu, center = TRUE, scale = TRUE)
}

# p2<-1+(cellheight_inch*dim(otu)[1])+(fontsize_inch*p1)
# p3<-dim(otu)[2]

if(!no_group){
    otu<-otu[match(rownames(group),rownames(otu)),]
}

if(!is.null(opt$map)&&!is.null(opt$order)){
    sample_order <- map[opt$order]
    sample_order <- na.omit(sample_order)
    sample_order <- rownames(sample_order)[order(sample_order[, 1])]
    otu<-otu[match(sample_order,rownames(otu)),]
    if("group"%in%ls()){
        # groups <-unique(group[,1][rownames(group)%in%sample_order])
        groups <-unique(group[,1][match(sample_order, rownames(group))])
    }
}


if(opt$bym){
    otu<-apply(otu,2,function(x){tapply(x,INDEX = group[,1],mean)})
    if(!is.null(opt$order)){
        otu<- otu[match(groups, rownames(otu)), ]
    }
}

if(opt$trans){
    otu<-t(otu)
}

cellwidth<-as.numeric(opt$cellwidth) #points
cellheight<-as.numeric(opt$cellheight) #points
fontsize<-as.numeric(opt$fontsize) #pixels?

# cellwidth_inch<-cellwidth/72
# cellheight_inch<-cellheight/72
# fontsize_inch<-fontsize/96/2
# fontsize_inch<-fontsize/96
# max_ft_len<-max_ft_len*fontsize_inch

# n_row<-dim(otu)[1]
# n_col<-dim(otu)[2]

# colname_len<-max(nchar(colnames(otu)))*fontsize_inch
# rowname_len<-max(nchar(rownames(otu)))*fontsize_inch

# ht <- n_row*cellheight_inch+colname_len+3
# wd <- n_col*cellwidth_inch+rowname_len+4

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




if(!no_group){
    if(is.null(opt$colors)){
        base_dir<-normalizePath(dirname(get_Rscript_filename()))
        source(paste(base_dir,"/piputils/get_colors.R", sep = ""))
        groups_color<-get_colors(opt$group, opt$map)
    }else{
        groups_color<-str_split(opt$colors, ",")[[1]]
    }
    names(groups_color)<-sort(unique(group[,1]))
    if("groups"%in%ls()){
        groups_color <- groups_color[groups]
    }
    rexp=sprintf("annotation_color=list(%s=groups_color)",opt$group)
    eval(parse(text = rexp))
}else{
    annotation_color=NA
}


cellcolors<-str_split(opt$cellcolors, ",")[[1]]

write.table(t(otu),paste(opt$out,"table.xls",sep = ""),sep = "\t",quote=FALSE,col.names=NA)

# 1 inch equals to 72 points
# 1 pixels equals to 0.75 points
# 1 point equals to 4/3 pixels
# 1 inch equals to 96 pixels
# pdf(paste(opt$out,"heatmap.pdf",sep = ""), height=ht, width=wd)

pheatmap(otu, file=paste(opt$out,"heatmap.pdf",sep = ""), annotation_row=annotation_row,
         annotation_col=annotation_col, fontsize=fontsize, border_color = "black",
         color = colorRampPalette(colors = cellcolors)(100),
         cluster_cols=cc,clustering_distance_cols="euclidean",
         cellwidth=cellwidth, cellheight=cellheight,
         # cutree_cols=3,
         # height=ht,width=wd,
         cluster_rows=cr,clustering_distance_rows="euclidean",annotation_colors=annotation_color)
# dev.off()

