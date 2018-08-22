library(optparse)
library(pheatmap)
library(psych)
library(stringr)

#######arguments
option_list <- list( 
    make_option(c("-i", "--input"), dest="otu",help="Specify the path of collapsed bacteria table",default=NULL),
    make_option(c("-m", "--map"),dest="map", help="Specify the path of mapping file",default=NULL),
    make_option(c("-e", "--exclude"),dest="ex", help="Specify the numeric variables excluded from plot and seprated by commas in mapping file",default="none"),
    make_option(c("-n", "--number"), dest="num",help="How many species do you want to plot",default=20),
    make_option(c("-o", "--output"),dest="out", help="Specify the path of output files",default="./")
    )

opt <- parse_args(OptionParser(option_list=option_list))

########prepare the cor data
ex<-str_split(opt$ex,",")[[1]]
dat <- read.table(opt$input,comment.char="",check.names=F,stringsAsFactors=F, header = TRUE, sep = "\t")
dat<-dat[!duplicated(dat[,1]),]
rownames(dat)=dat[,1]
annotation_row<-str_extract(dat[,length(dat)],"p__[^;]+")
dat=dat[,-c(1,length(dat))]


map<-read.table(opt$map,header = T,row.names=1,comment.char = "",check.names = F,stringsAsFactors = F)
notstr=c()
for(i in 1:length(map)){
	notstr[i]=!is.character(map[,i])
}
envdata<-map[,notstr]
if(ex[1]!="none"){envdata<-envdata[,!colnames(envdata)%%in%%ex]}
dat=t(dat)[match(rownames(envdata),rownames(t(dat))),]


LEN<-dim(envdata)[2]

dat_for_cor<-data.frame(envdata,dat, row.names = NULL, check.rows = T, check.names = F, stringsAsFactors = default.stringsAsFactors())


cor_allft<-corr.test(dat_for_cor,method ="spearman",adjust="none")
cor_allft_r<-data.frame(cor_allft$r,check.names = F)
cor_allft_r<-cor_allft_r[-c(1:LEN),]
cor_allft_r<-cor_allft_r[,-c(LEN+1:Dim(cor_allft_r)[2])]

cor_allft_p<-data.frame(cor_allft$p,check.names = F)
cor_allft_p<-cor_allft_p[-c(1:LEN),]
cor_allft_p<-cor_allft_p[,-c(LEN+1:Dim(cor_allft_p)[2])]


heat_s<-cor_allft_p

for(i in 1:dim(cor_allft_p)[1]){
  for(j in 1:dim(cor_allft_p)[2]){
    ifelse(cor_allft_p[i,j]<0.05,
           ifelse(cor_allft_p[i,j]>=0.01,heat_s[i,j]<-"*",
                  ifelse(cor_allft_p[i,jj]>=0.001,heat_s[i,j]<-"**",heat_s[i,j]<-"***")),
           heat_s[i,j]<-"")
  }
}

rownames(cor_allft_r)<-str_extract(rownames(cor_allft_r),"[^;]{1,100}")
annotation_row =data.frame(Phylum=annotation_row)
rownames(annotation_row) = rownames(cor_allft_r) 

pdf(paste(out,"/","Correlation_heatmap.pdf", sep = ""))
pheatmap(cor_allft_r,fontsize=15,annotation_row = annotation_row,border_color = "purple",
         display_numbers = heat_s,
         cluster_rows=T,clustering_distance_rows="correlation",
         cluster_cols=T,clustering_distance_cols="euclidean",
         clustering_method="centroid")
dev.off()
