library(optparse)

#######arguments
option_list <- list( 
    make_option(c("-i", "--input"),metavar="path", dest="otu",help="Specify the path of collapsed bacteria table",default=NULL),
    make_option(c("-m", "--map"),metavar="path",dest="map", help="Specify the path of mapping file",default=NULL),
    make_option(c("-e", "--exclude"),metavar="string",dest="ex", help="Specify the numeric variables excluded from plot and seprated by commas in mapping file",default="none"),
    make_option(c("-n", "--number"),metavar="int", dest="num",help="The number of most related species you want to plot, default is 20",default=30),
    make_option(c("-r", "--min-cor"),metavar="float", dest="minr",help="Min correlation coefficent to label significance, default is 0",default=0),
    make_option(c("-p", "--prefix"),metavar="str", dest="prefix",help="The prefix of output files, default if null",default=""),
    make_option(c("-o", "--output"),metavar="directory",dest="out", help="Specify the directory of output files",default="./")
    )

opt <- parse_args(OptionParser(option_list=option_list,description = "This script is used to relate the bacteria species and enviromental factors(numeric), and a heatmap is used to visualize the rank correlation coefficent"))

library(pheatmap)
library(psych)
library(stringr)


if(!dir.exists(opt$out)){dir.create(opt$out,recursive = T)}

opt$out<-paste(opt$out,"/",opt$prefix,sep="")
########prepare the cor data
ex<-str_split(opt$ex,",")[[1]]

dat <- read.table(as.character(opt$otu),quote="",comment.char="",check.names=F,stringsAsFactors=F, header = TRUE, sep = "\t")
dat<-dat[dat[,1]!="Others"&dat[,1]!="unclassified",]
dat<-dat[!duplicated(dat[,1]),]
rownames(dat)=dat[,1]
annotation_row<-str_extract(dat[,length(dat)],"p__[^;]+")
dat=dat[,-c(1,length(dat))]


map<-read.table(as.character(opt$map),sep="\t",na.strings="",header = T,row.names=1,comment.char = "",check.names = F,stringsAsFactors = F)
notstr=c()
for(i in 1:length(map)){
	notstr[i]=is.numeric(map[,i])
}

if(ex[1]!="none"){
    notstr<-notstr&!colnames(map)%in%ex
}

envdata<-map[colnames(map)[notstr]]

#envdata<-na.omit(envdata)
dat=t(dat)[match(rownames(envdata),colnames(dat)),]


######correlation statics

LEN<-ncol(envdata)
dat_for_cor<-data.frame(envdata,dat, check.rows = T, check.names = F)


cor_allft<-corr.test(dat_for_cor,method ="spearman",adjust="fdr")
cor_allft_r<-cor_allft$r
cor_allft_r<-cor_allft_r[-c(1:LEN),-c(LEN+1:dim(cor_allft_r)[2])]

cor_allft_p<-cor_allft$p
cor_allft_p<-cor_allft_p[-c(1:LEN),-c(LEN+1:dim(cor_allft_p)[2])]


write.table(cor_allft_r,paste(opt$out,"spearman_rank_correlation_matrix.txt",sep = ""),sep="\t")
write.table(cor_allft_p,paste(opt$out,"fdr_adjusted_p_value_matrix.txt",sep = ""),sep="\t")

cor_allft_p[abs(cor_allft_r)<opt$minr]<-1

ifcor<-colSums(t(cor_allft_p<0.05))
ifna<-colSums(t(is.na(cor_allft_p)))
order_pos<-order(ifcor,decreasing=TRUE)
order_pos<-order_pos[ifcor[order_pos]>=1&ifna[order_pos]==0]
selected_pos<-head(order_pos,opt$num)
cor_allft_r<-cor_allft_r[selected_pos,]
cor_allft_p<-cor_allft_p[selected_pos,]
annotation_row<-annotation_row[selected_pos]
#if(opt$num<dim(cor_allft_r)[1]){
#  ifcor<-colSums(t(cor_allft_p<0.05))
#  sel<-(ifcor>=sort(ifcor,T)[opt$num]&ifcor>0)
#  cor_allft_r<-cor_allft_r[sel,]
#  cor_allft_p<-cor_allft_p[sel,]
#  annotation_row<-annotation_row[sel]
#}else{
#  ifcor<-colSums(t(cor_allft_p<0.05))
#  sel<-(ifcor>0)
#  cor_allft_r<-cor_allft_r[sel,]
#  cor_allft_p<-cor_allft_p[sel,]
#  annotation_row<-annotation_row[sel]
#}


sig_label<-function(x){ifelse(x<0.001,"***",ifelse(x<0.01,"**",ifelse(x<0.05,"*","")))}
heat_s<-sig_label(cor_allft_p)
#heat_s<-cor_allft_p
#cor_allft_p[abs(cor_allft_r)<opt$minr]<-1
#for(i in 1:dim(cor_allft_p)[1]){
#  for(j in 1:dim(cor_allft_p)[2]){
#    ifelse(cor_allft_p[i,j]<0.05,
#           ifelse(cor_allft_p[i,j]>=0.01,heat_s[i,j]<-"*",
#                  ifelse(cor_allft_p[i,j]>=0.001,heat_s[i,j]<-"**",heat_s[i,j]<-"***")),
#           heat_s[i,j]<-"")
#  }
#}

rownames(cor_allft_r)<-str_extract(rownames(cor_allft_r),"[^;]{1,100}")
annotation_row =data.frame(Phylum=annotation_row)
rownames(annotation_row) = rownames(cor_allft_r) 
####corrlation plot
pdf(paste(opt$out,"Correlation_heatmap.pdf",sep = ""), height=5+0.5*dim(cor_allft_r)[1],width=6+1*dim(cor_allft_r)[2])
pheatmap(cor_allft_r,fontsize=12,annotation_row = annotation_row,border_color = "black",
         display_numbers = heat_s,fontsize_row =17,fontsize_col = 17,
         fontsize_number = 22,
         cluster_rows=T,clustering_distance_rows="correlation",
         cluster_cols=T,clustering_distance_cols="euclidean",
         clustering_method="centroid")
dev.off()
