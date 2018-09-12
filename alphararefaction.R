library(optparse)
library(reshape)
library(ggplot2)
library(stringr)

option_list <- list( 
    make_option(c("-i", "--input"), dest="i",help="Specify the path of collapsed bacteria table",default=NULL,metavar="path"),
    make_option(c("-o", "--output"), dest="o",help="The directory of output files",default=NULL,metavar="path")
  )


opt <- parse_args(OptionParser(option_list=option_list,description = "This script is used to plot rarefacation curve of alpha diversity and use the exported qiime2 rarefacation results path as input"))

if(!dir.exists(opt$o)){dir.create(opt$o,recursive = T)}

for (a in c("faith_pd","observed_otus","shannon")){
	dat <- read.table(paste(opt$i,"/",a,".csv",sep = ""),row.names=1,comment.char="",check.names=F,stringsAsFactors=F, header = TRUE, sep = ",")
	notstr=c()
	for(i in 1:length(dat)){
		notstr[i]=!is.character(dat[,i])
	}
	dat<-dat[,notstr]
	dat$SampleID<-rownames(dat)
	dat<-melt(dat,varnames="SampleID")
	dat$variable<-as.numeric(str_extract(dat$variable,"[0-9]+"))
	p <- ggplot(dat, aes(x = variable, y = value, color = SampleID)) +
	  geom_smooth(se=F, method = "lm",formula = y ~ log(x)) +  ### #使用log拟合
	  theme_bw() +xlab("Number of sequences")+ylab(a)+
	  theme(panel.grid=element_blank())
	ggsave(plot = p,filename = paste(opt$o,'/',a,"_rarefaction_curve.png",sep = ""),width = 8,height = 6,dpi = 300)
}
