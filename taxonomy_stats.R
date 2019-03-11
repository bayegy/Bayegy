library(optparse)

option_list <- list( 
  make_option(c("-i", "--input"),metavar="path", dest="otu",help="Specify the path of otu taxonomy table, with Consensus Lineage at last column",default=NULL),
  make_option(c("-o", "--output"),metavar="directory",dest="out", help="Specify the directory of output files",default="./"),
  make_option(c("-m", "--map"),metavar="path",dest="map", help="Specify the path of mapping file. Required",default=NULL),
  make_option(c("-c", "--category"),metavar="string",dest="group", help="Category to compare. Required",default=NULL),
  make_option(c("-b", "--by-group"),metavar="logical", dest="by",help="Pass this to use the group to plot barplot",default=FALSE)
)

opt <- parse_args(OptionParser(option_list=option_list,description = "This script is used to draw the relative stack barplot of species"))
library(phyloseq)
library(reshape)
library(ggplot2)
if(!dir.exists(opt$out)){dir.create(opt$out,recursive = T)}
opt$out<-paste(opt$out,"/",opt$prefix,sep="")

map<-read.table(opt$map,quote="",row.names = 1,na.strings = "",comment.char="",check.names=F,stringsAsFactors=F, header = TRUE, sep = "\t")
group<-na.omit(map[opt$group])
label_order<-rownames(group)[order(group)]


otu<-import_qiime(otufilename = opt$otu)

tb_otu<-data.frame(otu@otu_table)
tax_otu<-data.frame(otu@tax_table)
if(opt$by){
  group<-group[match(colnames(tb_otu),rownames(group)),]
  tb_otu<-apply(tb_otu,1,function(x){tapply(x,INDEX = group,sum)})
}

otu<-data.frame(tax_otu,OTU=rownames(tax_otu),tb_otu,check.rows = T,check.names = F)

spe_count<-function(data,class_level_column=c(1:8),cut=0){
  require(dplyr)
  sam<-data[,-class_level_column]
  class<-data[,class_level_column]
  func2<-function(z){return(sum(z>cut))}
  sub_ct<-function(x){
    func<-function(y){return(tapply(y,x,sum))}
    ss<-apply(sam,2,func)
    ifelse(length(unique(x))==1,return(as.numeric(ss>cut)),
           return(apply(ss,2,func2)))
  }
  out<-as.data.frame(apply(class,2,sub_ct))
  out<-data.frame(Sample=colnames(sam),out,check.names = F)
  return(out)
}

tb1<-spe_count(data = otu)
bardata1<-melt(tb1,id="Sample")

p<-ggplot(bardata1,aes(x=Sample,y=value,fill=variable))+
  geom_bar(stat = "identity")+
  guides(fill=guide_legend(reverse = F))+
  ylab("The number of taxas")+labs(fill="")+
  theme_bw()+
  scale_x_discrete(limits=label_order)+
  theme(text = element_text(size = 10),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        axis.line = element_line(),panel.border =  element_blank(),
        axis.text.x = element_text(angle = 45,size = 10,hjust = 1))+
  scale_y_continuous(expand = c(0, 0))

p1=length(unique(bardata1$Sample))

ggsave(plot = p,paste(opt$out,"number_of_taxas.png",sep = ""),width = p1*0.3+3,height = 7,dpi = 300)
write.csv(tb1,paste(opt$out,"number_of_taxas.csv",sep = ""),row.names = F)

otu_asign_count<-function(data,class_level_column=c(1:8),cut=0){
  require(dplyr)
  require(stringr)
  sam<-data[,-class_level_column]
  class<-data[,class_level_column]
  sub_ct<-function(x){
    x<-as.character(x)
    v1<-is.na(str_extract(x,"unclass"))+is.na(str_extract(x,"norank"))+is.na(str_extract(x,"uncultured"))+!is.na(x) 
    v1<-v1==4
    func2<-function(z){
      v2<-z>cut
      return(sum(v1+v2==2))
    }
    return(apply(sam,2,func2))
  }
  out<-as.data.frame(apply(class,2,sub_ct))
  out<-data.frame(Sample=rownames(out),out)
  return(out)
}

tb2<-otu_asign_count(data = otu)
bardata2<-melt(tb2,id="Sample")

p<-ggplot(bardata2,aes(x=Sample,y=value,fill=variable))+
  geom_bar(stat = "identity")+
  guides(fill=guide_legend(reverse = F))+
  ylab("The number of OTUs")+labs(fill="")+
  theme_bw()+
  scale_x_discrete(limits=label_order)+
  theme(text = element_text(size = 10),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        axis.line = element_line(),panel.border =  element_blank(),
        axis.text.x = element_text(angle = 45,size = 10,hjust = 1))+
  scale_y_continuous(expand = c(0, 0))

ggsave(plot = p,paste(opt$out,"number_of_assigned_otus.png",sep = ""),width = p1*0.3+3,height = 7,dpi = 300)
write.csv(tb2,paste(opt$out,"number_of_assigned_otus.csv",sep = ""),row.names = F)


