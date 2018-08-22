library("optparse")
library("pheatmap")
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
dat=dat[,-c(1,length(dat))]


map<-read.table("%s",header = T,row.names=1,comment.char = "",check.names = F,stringsAsFactors = F)
notstr=c()
for(i in 1:length(map)){
	notstr[i]=!is.character(map[,i])
}
envdata<-map[,notstr]
if(ex[1]!="none"){envdata<-envdata[,!colnames(envdata)%%in%%ex]}


dat=t(dat)[match(rownames(envdata),rownames(t(dat))),]

