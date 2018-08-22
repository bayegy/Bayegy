library(optparse)
library(reshape)
library(igraph)
library(psych)
require(RColorBrewer)
library(stringr)

option_list <- list( 
    make_option(c("-i", "--input"), dest="otu",help="Specify the path of collapsed bacteria table",default=NULL),
    make_option(c("-n", "--number"), dest="num",help="The number of most abundant species you want to plot",default=20),
    make_option(c("-o", "--output"),dest="out", help="Specify the path of output files",default="./")
    )

opt <- parse_args(OptionParser(option_list=option_list))


dat <- read.table(as.character(opt$otu),comment.char="",check.names=F,stringsAsFactors=F, header = TRUE, sep = "\t")
dat<-dat[dat[,1]!="Others",]
dat<-dat[!duplicated(dat[,1]),]
rownames(dat)=dat[,1]
annotation_row<-str_extract(dat[,length(dat)],"p__[^;]+")
dat=dat[,-c(1,length(dat))]
dat<-t(dat)

if(opt$num<dim(cor_allft_r)[1]){
  sel<-(colSums(dat)>=sort(colSums(dat),T)[opt$num])
  dat<-dat[,sel]
  annotation_row<-annotation_row[sel]
}

cor<-corr.test(dat,method = "spearman",adjust="fdr")
p<-cor$p
r<-cor$r
r[cor$p>=0.05] <- 0 
igraph<-graph_from_adjacency_matrix(r,mode="undirected",weighted=TRUE,diag=FALSE)
ss1<-colSums(tb)
V(igraph)$size<-(log(ss1+min(ss1),10)+(-log(min(ss1))))*2
igraph.weight = E(igraph)$weight
E.color = igraph.weight
E.color = ifelse(E.color>0,rainbow(3)[1],ifelse(E.color<0, rainbow(3)[3],"grey"))
E(igraph)$color = as.character(E.color)
E(igraph)$width = abs(igraph.weight)*3
V(igraph)$color<-brewer.pal(12,"Set3")[1]

bad.vs<-V(igraph)[degree(igraph)==0]
igraph<-delete.vertices(igraph,bad.vs)


pdf(paste(opt$out,"/","Correlation_network.pdf",sep = ""), height=14,width=12)
plot(igraph,vertex.label.color="black",edge.width=1,
     edge.curved=TRUE,margin=c(0,0,0,0))
dev.off()

