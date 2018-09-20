library(optparse)

#######arguments
option_list <- list( 
    make_option(c("-i", "--input"),metavar="path", dest="otu",help="Specify the path of collapsed bacteria table",default=NULL),
    make_option(c("-m", "--map"),metavar="path",dest="map", help="Specify the path of mapping file",default=NULL),
    make_option(c("-e", "--exclude"),metavar="string",dest="ex", help="Specify the numeric variables excluded from plot and seprated by commas in mapping file",default="none"),
    make_option(c("-n", "--number"),metavar="int", dest="num",help="The number of most related species you want to plot, default is 20",default=20),
    make_option(c("-r", "--min-cor"),metavar="float", dest="minr",help="Min correlation coefficent to label significance, default is 0.4",default=0.4),
    make_option(c("-o", "--output"),metavar="directory",dest="out", help="Specify the directory of output files",default="./")
    )

opt <- parse_args(OptionParser(option_list=option_list,description = "This script is used to relate the bacteria species and enviromental factors(numeric), and a heatmap is used to visualize the rank correlation coefficent"))

library(pheatmap)
library(psych)
library(stringr)

########prepare the cor data
ex<-str_split(opt$ex,",")[[1]]

dat <- read.table(as.character(opt$otu),comment.char="",check.names=F,stringsAsFactors=F, header = TRUE, sep = "\t")
dat<-dat[dat[,1]!="Others"&dat[,1]!="unclassified",]
dat<-dat[!duplicated(dat[,1]),]
rownames(dat)=dat[,1]
annotation_row<-str_extract(dat[,length(dat)],"p__[^;]+")
dat=dat[,-c(1,length(dat))]


map<-read.table(as.character(opt$map),sep="\t",na.strings="",header = T,row.names=1,comment.char = "",check.names = F,stringsAsFactors = F)
notstr=c()
for(i in 1:length(map)){
	notstr[i]=!is.character(map[,i])
}
envdata<-map[,notstr]
if(ex[1]!="none"){envdata<-envdata[,!colnames(envdata)%in%ex]}
dat=t(dat)[match(rownames(envdata),rownames(t(dat))),]




######correlation statics
if(!dir.exists(opt$out)){dir.create(opt$out,recursive = T)}

LEN<-dim(envdata)[2]
dat_for_cor<-data.frame(envdata,dat, row.names = NULL, check.rows = T, check.names = F, stringsAsFactors = default.stringsAsFactors())


cor_allft<-corr.test(dat_for_cor,method ="spearman",adjust="fdr")
cor_allft_r<-data.frame(cor_allft$r,check.names = F)
cor_allft_r<-cor_allft_r[-c(1:LEN),]
cor_allft_r<-cor_allft_r[,-c(LEN+1:dim(cor_allft_r)[2])]

cor_allft_p<-data.frame(cor_allft$p,check.names = F)
cor_allft_p<-cor_allft_p[-c(1:LEN),]
cor_allft_p<-cor_allft_p[,-c(LEN+1:dim(cor_allft_p)[2])]


if(opt$num<dim(cor_allft_r)[1]){
  ifcor<-colSums(t(cor_allft_p<0.05))
  sel<-(ifcor>=sort(ifcor,T)[opt$num]&ifcor>0)
  cor_allft_r<-cor_allft_r[sel,]
  cor_allft_p<-cor_allft_p[sel,]
  annotation_row<-annotation_row[sel]
}else{
  ifcor<-colSums(t(cor_allft_p<0.05))
  sel<-(ifcor>0)
  cor_allft_r<-cor_allft_r[sel,]
  cor_allft_p<-cor_allft_p[sel,]
  annotation_row<-annotation_row[sel]
}

write.table(cor_allft_r,paste(opt$out,"/","spearman_rank_correlation_matrix.txt",sep = ""),sep="\t")
write.table(cor_allft_p,paste(opt$out,"/","fdr_adjusted_p_value_matrix.txt",sep = ""),sep="\t")


heat_s<-cor_allft_p
cor_allft_p[abs(cor_allft_r)<opt$minr]<-1
for(i in 1:dim(cor_allft_p)[1]){
  for(j in 1:dim(cor_allft_p)[2]){
    ifelse(cor_allft_p[i,j]<0.05,
           ifelse(cor_allft_p[i,j]>=0.01,heat_s[i,j]<-"*",
                  ifelse(cor_allft_p[i,j]>=0.001,heat_s[i,j]<-"**",heat_s[i,j]<-"***")),
           heat_s[i,j]<-"")
  }
}

rownames(cor_allft_r)<-str_extract(rownames(cor_allft_r),"[^;]{1,100}")
annotation_row =data.frame(Phylum=annotation_row)
rownames(annotation_row) = rownames(cor_allft_r) 
####corrlation plot
pdf(paste(opt$out,"/","Correlation_heatmap.pdf",sep = ""), height=5+0.5*dim(cor_allft_r)[1],width=7+0.3*dim(cor_allft_r)[2])
pheatmap(cor_allft_r,fontsize=12,annotation_row = annotation_row,border_color = "black",
         display_numbers = heat_s,fontsize_row =17,fontsize_col = 17,
         fontsize_number = 22,
         cluster_rows=T,clustering_distance_rows="correlation",
         cluster_cols=T,clustering_distance_cols="euclidean",
         clustering_method="centroid")
dev.off()
