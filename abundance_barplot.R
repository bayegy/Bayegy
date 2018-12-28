library(optparse)

#######arguments
option_list <- list( 
    make_option(c("-i", "--input"),metavar="path", dest="otu",help="Specify the path of collapsed bacteria table",default=NULL),
    make_option(c("-m", "--map"),metavar="path",dest="map", help="Specify the path of mapping file. Required",default=NULL),
    make_option(c("-c", "--category"),metavar="string",dest="group", help="Category to compare. Required",default=NULL),
    make_option(c("-p", "--prefix"),metavar="str", dest="prefix",help="The prefix of output files, default if null",default=""),
    make_option(c("-n", "--number"),metavar="int", dest="num",help="The number of most related species you want to plot, default is 20",default=20),
    make_option(c("-a", "--add-posix"),action = "store_true", dest="add",help="add posix to the duplicated taxons"),
    make_option(c("-o", "--output"),metavar="directory",dest="out", help="Specify the directory of output files",default="./")
    )

opt <- parse_args(OptionParser(option_list=option_list,description = "This script is used to relate the bacteria species and enviromental factors(numeric), and a heatmap is used to visualize the rank correlation coefficent"))

library(ggplot2)
library(stringr)
library(reshape)
library(RColorBrewer)
display.brewer.all()


if(!dir.exists(opt$out)){dir.create(opt$out,recursive = T)}
opt$out<-paste(opt$out,"/",opt$prefix,sep="")

otu <- read.table(opt$otu,quote="",comment.char="",check.names=F,stringsAsFactors=F, header = TRUE, sep = "\t")
#otu <- read.table("otu_table.Genus.relative.txt",quote="",comment.char="",check.names=F,stringsAsFactors=F, header = TRUE, sep = "\t")

map<-read.table(opt$map,quote="",row.names = 1,na.strings = "",comment.char="",check.names=F,stringsAsFactors=F, header = TRUE, sep = "\t")
#map<-read.table("mapping_file.txt",quote="",row.names = 1,na.strings = "",comment.char="",check.names=F,stringsAsFactors=F, header = TRUE, sep = "\t")

group<-na.omit(map[opt$group])
label_order<-rownames(group)[order(group)]

add_posix<-function(z){
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


if(opt$add){
  otu[,1]<-add_posix(otu[,1])
}else{
  otu<-otu[!duplicated(otu[,1]),]
}

#otu<-otu[otu[,1]!="Others"&otu[,1]!="unclassified",]
rownames(otu)<-otu[,1]
otu<-t(otu[,-c(1,ncol(otu))])*100

otu<-otu[match(rownames(group),rownames(otu)),]

deod<-order(colSums(otu),decreasing = T)
sel<-head(deod,as.numeric(opt$num))
other<-deod[(as.numeric(opt$num)+1):length(deod)]
otu<-data.frame(otu[,sel],Other=apply(otu[,other],1,sum),check.names = F,check.rows = T)


#otu<-data.frame(otu,group)
otu$id<-rownames(otu)
#group<-map["Group1"]
#opt$group<-"Group1"
#otu<-melt(otu,id.vars = c("Group1","id"))
otu<-melt(otu,id.vars = "id")

pallet<-c(brewer.pal(12,"Paired"),brewer.pal(8,"Set2")[-c(7,8)],brewer.pal(8,"Dark2"),brewer.pal(12,"Set3"),brewer.pal(8,"Accent"),brewer.pal(11,"Spectral"))

p<-ggplot(otu,aes(x=id,y=value,fill=variable))+geom_bar(stat = "identity",width = 0.7)+
  guides(fill=guide_legend(title = NULL))+
  scale_fill_manual(values = pallet)+
  scale_x_discrete(limits=label_order)+
  xlab("")+ylab("Sequence Number Percent(%)")+theme_bw()+
  theme(text = element_text(size = 10),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        axis.line = element_line(),panel.border =  element_blank(),
        axis.text.x = element_text(angle = 90,size = 10,vjust = 0.5,hjust = 1))

wd<-length(label_order)*0.2+4

ggsave(plot = p,paste(opt$out,"barplot.pdf",sep = ""),width = wd,height = 7,dpi = 300)

