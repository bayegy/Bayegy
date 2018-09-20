library(optparse)


option_list <- list( 
    make_option(c("-i", "--input"),metavar="path", dest="otu",help="Specify the path of collapsed bacteria table",default=NULL),
    make_option(c("-n", "--number"),metavar="int", dest="num",help="The number of most abundant species you want to plot",default=50,type="integer"),
    make_option(c("-c", "--cut"),metavar="float",dest="cut", help="Specify the threshold of correlation coefficent",default=0.6),
    make_option(c("-o", "--output"),metavar="directory",dest="out", help="Specify the directory of output files",default="./")  
    )

opt <- parse_args(OptionParser(option_list=option_list))

library(igraph)
library(psych)
library(stringr)


if(!dir.exists(opt$out)){dir.create(opt$out,recursive = T)}
dat <- read.table(as.character(opt$otu),comment.char="",check.names=F,stringsAsFactors=F, header = TRUE, sep = "\t")
dat<-dat[dat[,1]!="Others"&dat[,1]!="unclassified",]
dat<-dat[!duplicated(dat[,1]),]
rownames(dat)=dat[,1]
annotation_row<-str_replace(str_extract(dat[,length(dat)],"p__[^;]+"),"p__","")
dat=dat[,-c(1,length(dat))]
dat<-t(dat)

if(opt$num<dim(dat)[2]){
  sel<-(colSums(dat)>=sort(colSums(dat),T)[opt$num])
  dat<-dat[,sel]
  annotation_row<-annotation_row[sel]
}


unig<-unique(annotation_row)
asign_rainbow_cor<-function(x){
	unig<-unique(x)
	color_out<-x
	for (i in 1:length(unig)) {
		color_out[color_out==unig[i]]<-rainbow(length(unig))[i]
	}
	return(color_out)
}


cor<-corr.test(dat,method = "spearman",adjust="fdr")
p<-cor$p
r<-cor$r

write.table(r,paste(opt$out,"/","spearman_rank_correlation_matrix.txt",sep = ""),sep="\t")
write.table(p,paste(opt$out,"/","fdr_adjusted_p_value_matrix.txt",sep = ""),sep="\t")


r[cor$p>=0.05|abs(cor$r)<opt$cut] <- 0 
igraph<-graph_from_adjacency_matrix(r,mode="undirected",weighted=TRUE,diag=FALSE)
ss<-scale(colSums(dat))
V(igraph)$size<-(ss-min(ss)+5)*2

#(log(ss1+min(ss1),10)+(-log(min(ss1))))*5
igraph.weight = E(igraph)$weight
E.color = igraph.weight
E.color = ifelse(E.color>0,rainbow(3)[1],ifelse(E.color<0, rainbow(3)[3],"grey"))
E(igraph)$color = as.character(E.color)
E(igraph)$width = abs(igraph.weight)*3
V(igraph)$color<-asign_rainbow_cor(annotation_row)


bad.vs<-V(igraph)[degree(igraph)==0]
#igraph<-delete.vertices(igraph,bad.vs)




pdf(paste(opt$out,"/","Correlation_network.pdf",sep = ""), height=12,width=12)
plot(igraph,vertex.label.color="black",edge.width=1,layout=layout_on_grid,
     edge.curved=F,margin=c(0,0,0,0.4),vertex.label.cex=1)

legend(title="Phylum",x="right",y = "right",bty="n",col = rainbow(length(unig)),pt.bg = rainbow(length(unig)),pch = rep(16,length(unig)),pt.cex=2,legend = unig)
dev.off()

