library(optparse)

#######arguments
option_list <- list( 
    make_option(c("-i", "--input"),metavar="path", dest="otu",help="Specify the path of collapsed bacteria table",default=NULL),
    make_option(c("-m", "--map"),metavar="path",dest="map", help="Specify the path of mapping file",default=NULL),
    make_option(c("-e", "--exclude"),metavar="string",dest="ex", help="Specify the numeric variables excluded from plot and seprated by commas in mapping file",default="none"),
    make_option(c("-n", "--number"),metavar="int", dest="num",help="The number of most related species you want to plot, default is 20",default=30),
    make_option(c("-r", "--min-cor"),metavar="float", dest="minr",help="Min correlation coefficent to label significance, default is 0",default=0),
    make_option(c("-p", "--prefix"),metavar="str", dest="prefix",help="The prefix of output files, default if null",default=""),
    make_option(c("-W", "--cellwidth"),metavar="int", dest="cellwidth",help="cell width in points",default=32),
    make_option(c("-H", "--cellheight"),metavar="int", dest="cellheight",help="cell height in points",default=20),
    make_option(c("-C", "--cellcolors"),metavar="string",dest="cellcolors", help="Comma seprated 3 colors",default="default"),
    make_option(c("-S", "--fontsize"),metavar="int", dest="fontsize",help="font size in pixels",default=10),
    make_option(c("-A", "--adjustp"),metavar="logical", dest="adjustp",help="FDR adjust the p value or not", default=FALSE),
    make_option(c("-o", "--output"),metavar="directory",dest="out", help="Specify the directory of output files",default="./")
    )

opt <- parse_args(OptionParser(option_list=option_list,description = "This script is used to relate the bacteria species and enviromental factors(numeric), and a heatmap is used to visualize the rank correlation coefficent"))

library(pheatmap)
library(psych)
library(stringr)
library("RColorBrewer")


if(!dir.exists(opt$out)){dir.create(opt$out,recursive = T)}

opt$out<-paste(opt$out,"/",opt$prefix,sep="")
########prepare the cor data


dat <- read.table(as.character(opt$otu),quote="",comment.char="",check.names=F,stringsAsFactors=F, header = TRUE, sep = "\t")
dat<-dat[dat[,1]!="Others"&dat[,1]!="unclassified",]



dat<-dat[!duplicated(dat[,1]),]
rownames(dat)=dat[,1]
annotation_row<-str_extract(dat[,length(dat)],"p__[^;]+")

is_num=c()
for(i in 1:ncol(dat)){
  is_num[i]=is.numeric(dat[,i])
}
dat=dat[, is_num]


if(nrow(dat)>200){
    dat = dat[head(order(colSums(t(dat)),decreasing=TRUE), 200),]
}


map<-read.table(as.character(opt$map),sep="\t",na.strings="",header = T,row.names=1,comment.char = "",check.names = F,stringsAsFactors = F)
colnames(map)[is.na(colnames(map))]<-"NA"


notstr=c()
for(i in 1:length(map)){
	notstr[i]=is.numeric(map[,i])
}

#ex<-str_split(opt$ex,",")[[1]]
#if(ex[1]!="none"){
#    notstr<-notstr&!colnames(map)%in%ex
#}
if(opt$ex!="none"){
  pre_ex <- str_split(opt$ex,">")[[1]]
  if(length(pre_ex)>1){
    flag <- pre_ex[1]
    ex <- pre_ex[2]
  }else{
    flag<-"exclude"
    ex <- pre_ex[1]
  }
  ex<-str_split(ex,",")[[1]]
  if(flag=="exclude"){
    notstr=notstr&(!colnames(map)%in%ex)
  }else{
    notstr=notstr&(colnames(map)%in%ex)
  }
}


envdata<-map[colnames(map)[notstr]]

#envdata<-na.omit(envdata)
dat=t(dat)[match(rownames(envdata),colnames(dat)),]


######correlation statics

LEN<-ncol(envdata)
dat_for_cor<-data.frame(envdata,dat, check.rows = T, check.names = F)


cor_allft<-corr.test(dat_for_cor,method ="spearman", adjust=ifelse(as.logical(opt$adjustp), "fdr", "none"))


cor_allft_r<-cor_allft$r
cor_allft_r<-cor_allft_r[-c(1:LEN),-c(LEN+1:dim(cor_allft_r)[2])]

cor_allft_p<-cor_allft$p

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
cor_allft_p<-up_to_down(cor_allft_p)

cor_allft_p<-cor_allft_p[-c(1:LEN),-c(LEN+1:dim(cor_allft_p)[2])]



write.table(cor_allft_r,paste(opt$out,"spearman_rank_correlation_matrix.xls",sep = ""),sep="\t",col.names=NA)
write.table(cor_allft_p,paste(opt$out,"p_value_matrix.xls",sep = ""),sep="\t",col.names=NA)

cor_allft_p[abs(cor_allft_r)<opt$minr]<-1

ifcor<-colSums(t(cor_allft_p<0.05))
signif_p <- t(cor_allft_p)
signif_p[signif_p>=0.05]<-NA
ifminp <- colSums(signif_p, na.rm=T)
ifminp <- ifminp/ifcor
minp_factor <- 1/ifminp

ifna<-colSums(t(is.na(cor_allft_p)))
#print(ifna)
order_pos<-order(ifcor, minp_factor, decreasing=TRUE)
order_pos<-order_pos[ifcor[order_pos]>=1&ifna[order_pos]==0]
selected_pos<-head(order_pos, as.numeric(opt$num))
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

# p1<-max(nchar(rownames(cor_allft_r)))
# wd<-4+0.8*dim(cor_allft_r)[2]+p1*0.08
if(NA%in%annotation_row){
    # wd<-wd-2
    annotation_row=NA
}else{
    annotation_row =data.frame(Phylum=annotation_row)
    rownames(annotation_row) = rownames(cor_allft_r)     
}


cellwidth<-as.numeric(opt$cellwidth) #points
cellheight<-as.numeric(opt$cellheight) #points
fontsize<-as.numeric(opt$fontsize) #pixels?

if(opt$cellcolors=="default"){
  cellcolors<-colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100)
}else{
  cellcolors<-colorRampPalette(colors = str_split(opt$cellcolors, ",")[[1]])(100)
}
####corrlation plot

# pdf(paste(opt$out,"Correlation_heatmap.pdf",sep = ""), height=4+0.35*dim(cor_allft_r)[1],width=wd)
pheatmap(cor_allft_r, file=paste(opt$out,"Correlation_heatmap.pdf",sep = ""),
         annotation_row = annotation_row, border_color = "black",
         display_numbers = heat_s, fontsize=fontsize,
         fontsize_number = 22,
         color = cellcolors, cellwidth=cellwidth, cellheight=cellheight,
         cluster_rows=T, clustering_distance_rows="euclidean",
         cluster_cols=T, clustering_distance_cols="euclidean")
# dev.off()
