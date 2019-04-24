library(optparse)

#######arguments
option_list <- list( 
    make_option(c("-i", "--input"),metavar="path", dest="otu",help="Specify the path of collapsed bacteria table",default=NULL),
    make_option(c("-m", "--map"),metavar="path",dest="map", help="Specify the path of mapping file. Required",default=NULL),
    make_option(c("-c", "--category"),metavar="string",dest="group", help="Category to compare. Required",default=NULL),
    make_option(c("-p", "--prefix"),metavar="str", dest="prefix",help="The prefix of output files, default if null",default=""),
    make_option(c("-n", "--number"),metavar="int", dest="num",help="The number of species needed to be plotted, default is 20",default=20),
    make_option(c("-a", "--add-postfix"),metavar="logical", dest="add",help="add postfix to the duplicated taxons",default=TRUE),
    make_option(c("-s", "--skip-first-line"),metavar="logical", dest="skip",help="If TRUE, skip the first line when read the input file.",default=FALSE),
    make_option(c("-b", "--by-groupMean"),metavar="logical", dest="bym",help="Pass this to use the group mean to plot barplot",default=FALSE),
    make_option(c("-o", "--output"),metavar="directory",dest="out", help="Specify the directory of output files",default="./")
    )

opt <- parse_args(OptionParser(option_list=option_list,description = "This script is used to draw the relative stack barplot of species"))

library(ggplot2)
library(stringr)
library(reshape)
library(RColorBrewer)
#display.brewer.all()


if(!dir.exists(opt$out)){dir.create(opt$out,recursive = T)}
opt$out<-paste(opt$out,"/",opt$prefix,sep="")

#otu <- read.table("otu_table.Genus.relative.txt",quote="",comment.char="",check.names=F,stringsAsFactors=F, header = TRUE, sep = "\t")
if(as.logical(opt$skip)){
  otu <-read.table(opt$otu,quote="",comment.char="",check.names=F,stringsAsFactors=F, header = TRUE, sep = "\t",skip = 1)
}else{
  otu <-read.table(opt$otu,quote="",comment.char="",check.names=F,stringsAsFactors=F, header = TRUE, sep = "\t")
}

if(!is.numeric(otu[,ncol(otu)])){
  otu<-otu[,-ncol(otu)]
}

if(!is.null(opt$map)&!is.null(opt$group)){
  map<-read.table(opt$map,quote="",row.names = 1,na.strings = "",comment.char="",check.names=F,stringsAsFactors=F, header = TRUE, sep = "\t")
#map<-read.table("mapping_file.txt",quote="",row.names = 1,na.strings = "",comment.char="",check.names=F,stringsAsFactors=F, header = TRUE, sep = "\t")
  map<-map[order(rownames(map)),]
  group<-na.omit(map[opt$group])
  label_order<-rownames(group)[order(group)]
}else{
  label_order<-ordered(colnames(otu)[-1])
}

add_postfix<-function(z){
  i<-1
  ad<-function(x){
    x[duplicated(x)]<-str_replace(x[duplicated(x)],"_[0-9]+$","")
    x[duplicated(x)]<-paste(x[duplicated(x)],"_",i,sep = "")
    return(x)
  }
  y<-z
  while(sum(duplicated(y))!=0){
    y<-ad(y)
    i<-i+1
  }
  return(y)
}


sum_abundance<-colSums(t(otu[,-1]))
otu<-otu[order(sum_abundance,decreasing = T),]

if(as.logical(opt$add)){
  otu[,1]<-add_postfix(otu[,1])
}else{
  otu<-otu[!duplicated(otu[,1]),]
}

#otu<-otu[otu[,1]!="Others"&otu[,1]!="unclassified",]
rownames(otu)<-otu[,1]

otu<-t(apply(otu[,-1],2,function(x){x/sum(x)}))*100

#otu<-t(otu[,-1])*100
if(!is.null(opt$map)&!is.null(opt$group)){
  otu<-otu[match(rownames(group),rownames(otu)),]
}

if(!is.null(opt$map)&!is.null(opt$group)&opt$bym){
  otu<-data.frame(apply(otu,2,function(x){tapply(x,INDEX = group[,1],mean)}),check.names=F)
  otu<-data.frame(t(apply(otu,1,function(x){x/sum(x)}))*100,check.names=F)
  label_order<-ordered(rownames(otu))
  #print(otu)
}
#deod<-order(colSums(otu),decreasing = T)
#sel<-head(deod,as.numeric(opt$num))
#other<-deod[(as.numeric(opt$num)+1):length(deod)]
#otu<-data.frame(otu[,sel],Other=apply(otu[,other],1,sum),check.names = F,check.rows = T)
num<-as.numeric(opt$num)
if(num<ncol(otu)-1){
  otu<-data.frame(otu[,1:num],Other=apply(otu[,(num+1):ncol(otu)],1,sum),check.names = F,check.rows = T)
}else{
  otu<-data.frame(otu,check.names = F,check.rows = T)
}


otu_out<-t(otu)/100
#otu<-data.frame(otu,group)
otu$id<-rownames(otu)
p1<-(max(nchar(colnames(otu)))*0.05+0.3)*ceiling(ncol(otu)/17)+2.5

#group<-map["Group1"]
#opt$group<-"Group1"
#otu<-melt(otu,id.vars = c("Group1","id"))
#reverse the order of species

otu<-otu[rev(colnames(otu))]

otu<-melt(otu,id.vars = "id")

pallet<-c(rev(brewer.pal(12,"Paired")),brewer.pal(8,"Set2")[-c(7,8)],brewer.pal(8,"Dark2"),brewer.pal(12,"Set3"),brewer.pal(8,"Accent"),brewer.pal(11,"Spectral"))

p<-ggplot(otu,aes(x=id,y=value,fill=variable))+geom_bar(stat = "identity",width = 0.7)+
  guides(fill=guide_legend(title = NULL))+
  scale_fill_manual(values = pallet)+
  scale_x_discrete(limits=label_order)+
  xlab("")+ylab("Sequence Number Percent(%)")+theme_bw()+
  theme(text = element_text(size = 10),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        axis.line = element_line(),panel.border =  element_blank(),
        axis.text.x = element_text(angle = 45,size = 10,hjust = 1))+
  scale_y_continuous(limits=c(0,101), expand = c(0, 0))


wd<-length(label_order)*0.2+p1
wd<-ifelse(wd<50,wd,49.9)
write.table(otu_out,paste(opt$out,"table.txt"),sep = "\t",quote=FALSE,row.names = TRUE,col.names = NA)
ggsave(plot = p,paste(opt$out,"barplot.pdf",sep = ""),width = wd,height = 7,dpi = 300)

