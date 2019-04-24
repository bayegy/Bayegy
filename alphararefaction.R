library(optparse)


option_list <- list( 
    make_option(c("-i", "--input"), dest="i",help="Specify the directory of exported qiime2 alpha-rarefacation.qzv",default=NULL,metavar="directory"),
    make_option(c("-o", "--output"), dest="o",help="The directory of output files",default=NULL,metavar="directory"),
    make_option(c("-m", "--map"),metavar="path",dest="map", help="Specify the path of mapping file. Required",default=NULL),
    make_option(c("-p", "--prefix"),metavar="str", dest="prefix",help="The prefix of output files, default if null",default=""),
    make_option(c("-c", "--category"),metavar="string",dest="group", help="Category to compare. Required",default=NULL)

  )


opt <- parse_args(OptionParser(option_list=option_list,description = "This script is used to plot rarefacation curve of alpha diversity and use the exported qiime2 rarefacation results path as input"))

library(getopt)
library(reshape)
library(ggplot2)
library(stringr)
if(!dir.exists(opt$o)){dir.create(opt$o,recursive = T)}
opt$o<-paste(opt$o,"/",opt$prefix,sep="")

groupmean<-ifelse(is.null(opt$map)|is.null(opt$group),FALSE,TRUE)


for (a in c("faith_pd","observed_otus","shannon")){
	dat <- read.table(paste(opt$i,"/",a,".csv",sep = ""),row.names=1,comment.char="",check.names=F,stringsAsFactors=F, header = TRUE, sep = ",")
	notstr=c()
	for(i in 1:length(dat)){
		notstr[i]=!is.character(dat[,i])
	}
	dat<-dat[,notstr]
	if(groupmean){
	    map<-read.table(opt$map,quote="",row.names = 1,na.strings = "",comment.char="",check.names=F,stringsAsFactors=F, header = TRUE, sep = "\t")
	    group<-na.omit(map[c(opt$group,"Description")])
	    group<-group[order(rownames(group)),]
	    group<-group[order(group[,1]),]
	    dat<-dat[match(rownames(group),rownames(dat)),]
	    dat<-data.frame(apply(dat,2,function(x){tapply(x,INDEX = group[,1],mean)}),check.names=F)
	}
	p_width=ceiling(nrow(dat)/16)*2+6
	dat$SampleID<-rownames(dat)
	dat<-melt(dat,varnames="SampleID")
	dat$variable<-as.numeric(str_extract(dat$variable,"[0-9]+"))
	p <- ggplot(dat, aes(x = variable, y = value, color = SampleID)) +
	  geom_smooth(se=F, method = "lm",formula = y ~ log(x)) +  ### #使用log拟合
	  theme_bw() +xlab("Number of sequences")+ylab(a)+
	  theme(panel.grid=element_blank())+
	  guides(color=guide_legend(title = NULL))
	if(groupmean){
		base_dir<-normalizePath(dirname(get_Rscript_filename()))
		source(paste(base_dir,"/piputils/get_colors.R", sep = ""))
		groups_color<-get_colors(opt$group, opt$map)
		p<-p+scale_colour_manual(values = groups_color)
	}
	ggsave(plot = p,filename = paste(opt$o,a,"_rarefaction_curve.pdf",sep = ""),width = p_width,height = 6,dpi = 300)
}
