#!/usr/bin/env Rscript
library(optparse)


option_list <- list(
    make_option(c("-i", "--input"),metavar="path", dest="otu",help="Specify the path of collapsed bacteria table",default=NULL),
    make_option(c("-n", "--number"),metavar="int", dest="num",help="The number of most abundant species you want to plot",default=50,type="integer"),
    make_option(c("-c", "--cut"),metavar="float",dest="cut", help="Specify the threshold of correlation coefficent",default=0.6),
    make_option(c("-l", "--layout"),metavar="string",dest="layout", help="Network layout", default='layout_on_grid'),
    make_option(c("-r", "--coefficient"),metavar="string",dest="coefficient", help="type of the correlation coefficient", default='spearman'),
    make_option(c("-e", "--edge-colors"),metavar="string",dest="colors", help="edge colors in the formation of [pos color],[neg color]", default='#FF0000,#0000FF'),
    make_option(c("-d", "--delete"),metavar="logical",dest="delete", help="delete isolated vertices", default=FALSE),
    make_option(c("-C", "--curve"),metavar="logical",dest="curve", help="curve the edges", default=TRUE),
    make_option(c("-a", "--adjust"),metavar="logical",dest="adjust", help="adjust the false discovery rate", default=FALSE),
    make_option(c("-o", "--output"),metavar="directory",dest="out", help="Specify the directory of output files",default="./")
    )

opt <- parse_args(OptionParser(option_list=option_list))

library(igraph)
library(psych)
library(stringr)
asign_rainbow_cor<-function(x){
  unig<-unique(x)
  color_out<-x
  rainbow_colors<-rainbow(length(unig))
  for (i in 1:length(unig)) {
    color_out[color_out==unig[i]]<-rainbow_colors[i]
  }
  return(color_out)
}

up_to_down<-function(m){
  d2<-ncol(m)
  for(i in 1:d2){
    for (j in 1:d2) {
      if(i>j){
        m[i,j]<-m[j,i]
      }
    }
  }
  return(m)
}

if(!dir.exists(opt$out)){dir.create(opt$out,recursive = T)}
dat <- read.table(as.character(opt$otu),comment.char="",check.names=F,stringsAsFactors=F, header = TRUE, sep = "\t")
dat<-dat[dat[,1]!="Others"&dat[,1]!="unclassified",]
dat<-dat[!duplicated(dat[,1]),]
rownames(dat)=dat[,1]


annotation<-FALSE

if(!is.numeric(dat[,ncol(dat)])){
  annotation_row<-str_replace(str_extract(dat[, ncol(dat)],"p__[^;]+"),"p__","")
  unig<-unique(annotation_row)
  has_phylum <- TRUE
  if(sum(!is.na(unig))<1){
    annotation_row <- dat[, ncol(dat)]
    has_phylum <- FALSE
  }
  dat=dat[,-c(1,ncol(dat))]
  annotation<-TRUE
}else{
  dat=dat[,-1]
}
dat<-t(dat)


if(opt$num<dim(dat)[2]){
  sel<-(colSums(dat)>=sort(colSums(dat),T)[opt$num])
  dat<-dat[,sel]
  if(annotation){
    annotation_row<-annotation_row[sel]
  }
}

# correlation of top species
cor<-corr.test(dat, method = opt$coefficient, adjust=ifelse(as.logical(opt$adjust), "fdr", "none"))
p<-cor$p
r<-cor$r
p<-up_to_down(p)

write.table(r,paste(opt$out,"/", opt$coefficient, "_correlation_matrix.xls",sep = ""), sep="\t", col.names = NA)
write.table(p,paste(opt$out,"/","p_value_matrix.xls",sep = ""), sep="\t", col.names = NA)

r[p>=0.05|abs(r)<opt$cut] <- 0
igraph<-graph_from_adjacency_matrix(r,mode="undirected",weighted=TRUE,diag=FALSE)
ss<-scale(colSums(dat))
V(igraph)$size<-(ss-min(ss)+5)*2

#(log(ss1+min(ss1),10)+(-log(min(ss1))))*5
igraph.weight = E(igraph)$weight
E.color = igraph.weight

edge_colors <- str_split(opt$colors, ',')[[1]]
E.color = ifelse(E.color>0,edge_colors[1], ifelse(E.color<0, edge_colors[2], "grey"))


E(igraph)$color = as.character(E.color)
E(igraph)$width = abs(igraph.weight)*3

if(annotation){
  V(igraph)$color<-asign_rainbow_cor(annotation_row)
}else{
  V(igraph)$color<-'#98F5FF'
}

if(opt$delete){
  bad.vs<-V(igraph)[degree(igraph)==0]
  igraph<-delete.vertices(igraph,bad.vs)
}
#

pdf(paste(opt$out,"/","Correlation_network.pdf",sep = ""), height=10,width=14)

par(fig=c(0,0.83,0,1),xpd=F)
plot(igraph,vertex.label.color="black",layout=get(opt$layout),
     edge.curved=as.logical(opt$curve),margin=c(0,0,0,0)+0.06,vertex.label.cex=0.9)

if(annotation){
  unig<-unique(annotation_row)
  par(fig=c(0.7,1,0,1),xpd=F)
  legend('center',xpd = F,title=ifelse(has_phylum,"Phylum",""),bty="n",col = rainbow(length(unig)),pch = rep(16,length(unig)),pt.cex=2,legend = unig,horiz = FALSE)
}
dev.off()

