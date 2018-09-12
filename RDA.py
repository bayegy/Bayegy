#!/usr/bin/python
# -*- coding: utf-8 -*-
import argparse
import re,sys,os

#*********************************************************************** *********************************************************************************
#argument:
p = argparse.ArgumentParser(description="This script is used to plot RDA of species. The numeric enviroment factors must be encluded in maping file. The categories will be filterd before RDA")
p.add_argument('-i', '--input', dest = 'input', metavar = '<file>',
			help = 'taxonomic count data file')
p.add_argument('-o', '--output', dest = 'output', metavar = '<Directory>', default = './',
			help = 'given an output directory')
p.add_argument('-m', '--metadata', dest = 'meta', metavar = '<file>',
			help = 'sample metadata file')
p.add_argument('-g', '--group', dest = 'group', metavar = '<str>',
			help = 'column name in sample-metadata file')
p.add_argument('-n', '--number', dest = 'number', metavar = '<int>', default = '15',
			help = 'specify how many species to be display, defaulf is 15')
p.add_argument('-e', '--exclude', dest = 'exclude', metavar = '<str>', default = 'none',
			help = 'specify numeric variables excluded from rda seprated by commas,use "none" if all numeric variables is expected')
options = p.parse_args()

if not options.input and options.group:
	p.error("must have argument -i")
	sys.exit()
else:
	pass


rscript = open('rda.R', 'w')

print ('''
library(vegan)
library(stringr)
library(ggrepel)
library(ggplot2)
library(RColorBrewer)
ex<-str_split("%s",",")[[1]]
dat <- read.table("%s", header = TRUE, sep = "\\t",comment.char = "",check.names = F)
dat<-dat[!duplicated(dat[,1]),]

rownames(dat)=dat[,1]
map<-read.table("%s",header = T,na.strings="",row.names=1,comment.char = "",check.names = F,stringsAsFactors = F)
groups<-map["%s"]

notstr=c()
for(i in 1:length(map)){
	notstr[i]=!is.character(map[,i])
}
envdata<-map[,notstr]
if(ex[1]!="none"){envdata<-envdata[,!colnames(envdata)%%in%%ex]}

dat=dat[,-c(1,length(dat))]
dat=t(dat)[match(rownames(groups),rownames(t(dat))),]

#clean na
groups<-groups[!is.na(groups),]
dat<-dat[!is.na(groups),]
envdata<-envdata[!is.na(groups),]


dca <- decorana(veg = dat)
dcam <- max(dca$rproj)
if (dcam > 4){
	cca <- cca(formula=dat~.,data=envdata, scale = TRUE, na.action = na.exclude)
	pre <- "CCA"
}else{
	cca <- rda(formula=dat~.,data=envdata,scale = TRUE, na.action = na.exclude)
	pre <- "RDA"
}

path="%s"
ccascore <- scores(cca)
write.table(ccascore$sites, file = paste(path,"/", pre, ".sample.txt", sep = ""), sep = "\\t")
write.table(ccascore$species, file = paste(path,"/", pre, ".bacteria.txt", sep = ""), sep = "\\t")
envfit <- envfit(cca, envdata, permu = 2000, na.rm = TRUE)
rp <- cbind(as.matrix(envfit$vectors$r), as.matrix(envfit$vectors$pvals))
colnames(rp) <- c("r2", "Pr(>r)")
env <- cbind(envfit$vectors$arrows, rp)
write.table(as.data.frame(env), file = paste(path,"/", pre, ".envfit.txt", sep = ""),sep = "\\t")


new<-cca$CCA

#samples=data.frame(ccascore$sites)
samples=data.frame(new$u)
samples$id=rownames(samples)
samples=data.frame(samples,%s=groups)


#species=data.frame(ccascore$species)
species=data.frame(new$v)
species$id=rownames(species)
species<-species[match(colnames(dat),rownames(species)),]
sum<-scale(colSums(dat))[,1]
species$abundance<-sum-min(sum)+0.8
number<-%s+1
if (dim(dat)[2]<number){
  cut=-1
}else{
  cut=sort(colSums(dat),T)[number]
}
show_species<-species[colSums(dat)>cut,]
hide_species<-species[colSums(dat)<=cut,]
print(paste("The threshold of abundance is:",cut))



#envis=data.frame(envfit$vectors$arrows)#Envis seem to be the same length after permutation
envis=data.frame(new$biplot)

envis$id=rownames(envis)



pc1 = cca$CCA$eig[1]/sum(cca$CCA$eig) * 100
pc2 = cca$CCA$eig[2]/sum(cca$CCA$eig) * 100
xlab <- paste(pre, "1: ", round(pc1, digits = 2), "%%", sep = "")
ylab <- paste(pre, "2: ", round(pc2, digits = 2), "%%", sep = "")

if(pre=="RDA"){
envis[,1]<-envis[,1]*0.5
envis[,2]<-envis[,2]*0.5
p1<-ggplot(data=samples,aes(x=RDA1,y=RDA2)) + 
  geom_point(aes(x=RDA1,y=RDA2,color=%s,pch=%s),size=3) +
  geom_text_repel(aes(x=RDA1,y=RDA2,label=id),color="black",size=3)+  
  geom_text_repel(data=envis,aes(x=RDA1,y=RDA2,label=id),color="black",size=5) +  
  geom_hline(yintercept=0,linetype="dotted") + geom_vline(xintercept=0,linetype="dotted")+  
  theme_bw() + theme(panel.grid=element_blank())+xlab(xlab)+ylab(ylab)+
  labs(title=paste(pre," sample location plot",sep=""))+
  geom_segment(aes(x=0,y=0,xend = RDA1, yend = RDA2),data = envis,
                color=brewer.pal(8,"Accent")[6],size=1,
                arrow = arrow(length = unit(0.2,"cm")))+
  theme(title = element_text(size = 15,vjust = 0.5,hjust = 0.5),
        axis.text = element_text(size = 15),
        text = element_text(size = 15),
        axis.title = element_text(size = 15),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 15))

ggsave(paste(path,"/","%s_",pre,"_sample_location_plot.png",sep=""),plot=p1,width = 9,height = 7,dpi = 300)


p2<-ggplot(data=show_species,aes(x=RDA1,y=RDA2)) + 
  geom_point(aes(x=RDA1,y=RDA2,size=abundance),color=rainbow(3)[1]) + 
  geom_point(data=hide_species,aes(x=RDA1,y=RDA2,size=abundance),color="grey") +
  guides(size=F)+
  geom_text_repel(data=envis,aes(x=RDA1,y=RDA2,label=id),color="black",size=5) +
  geom_text_repel(aes(RDA1,RDA2,label=id)) +  
  geom_hline(yintercept=0,linetype="dotted") + geom_vline(xintercept=0,linetype="dotted")+  
  theme_bw() + theme(panel.grid=element_blank())+xlab(xlab)+ylab(ylab)+
  ggtitle(paste(pre," bacteria location plot",sep=""))+
  geom_segment(aes(x=0,y=0,xend = RDA1, yend = RDA2),data = envis,
                color=brewer.pal(8,"Accent")[6],size=1,
                arrow = arrow(length = unit(0.2,"cm")))+
  theme(title = element_text(size = 15,vjust = 0.5,hjust = 0.5),
        axis.text = element_text(size = 15),
        text = element_text(size = 15),
        axis.title = element_text(size = 15),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 15))

ggsave(paste(path,"/",pre,"_bacteria_location_plot.png",sep=""),plot=p2,width = 7,height = 7,dpi = 300)
}else{
envis[,1]<-envis[,1]*3.5
envis[,2]<-envis[,2]*3.5
p1<-ggplot(data=samples,aes(x=CCA1,y=CCA2)) + 
  geom_point(aes(x=CCA1,y=CCA2,color=%s,pch=%s),size=3) +
  geom_text_repel(aes(x=CCA1,y=CCA2,label=id),color="black",size=3)+  
  geom_text_repel(data=envis,aes(x=CCA1,y=CCA2,label=id),color="black",size=5) +  
  geom_hline(yintercept=0,linetype="dotted") + geom_vline(xintercept=0,linetype="dotted")+  
  theme_bw() + theme(panel.grid=element_blank())+xlab(xlab)+ylab(ylab)+
  labs(title=paste(pre," sample location plot",sep=""))+
  geom_segment(aes(x=0,y=0,xend = CCA1, yend = CCA2),data = envis,
                color=brewer.pal(8,"Accent")[6],size=1,
                arrow = arrow(length = unit(0.2,"cm")))+
  theme(title = element_text(size = 15,vjust = 0.5,hjust = 0.5),
        axis.text = element_text(size = 15),
        text = element_text(size = 15),
        axis.title = element_text(size = 15),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 15))

ggsave(paste(path,"/","%s_",pre,"_sample_location_plot.png",sep=""),plot=p1,width = 9,height = 7,dpi = 300)


p2<-ggplot(data=show_species,aes(x=CCA1,y=CCA2)) + 
  geom_point(aes(x=CCA1,y=CCA2,size=abundance),color=rainbow(3)[1]) + 
  geom_point(data=hide_species,aes(x=CCA1,y=CCA2,size=abundance),color="grey") +
  guides(size=F)+
  geom_text_repel(data=envis,aes(x=CCA1,y=CCA2,label=id),color="black",size=5) +
  geom_text_repel(aes(CCA1,CCA2,label=id)) +  
  geom_hline(yintercept=0,linetype="dotted") + geom_vline(xintercept=0,linetype="dotted")+  
  theme_bw() + theme(panel.grid=element_blank())+xlab(xlab)+ylab(ylab)+
  ggtitle(paste(pre," bacteria location plot",sep=""))+
  geom_segment(aes(x=0,y=0,xend = CCA1, yend = CCA2),data = envis,
                color=brewer.pal(8,"Accent")[6],size=1,
                arrow = arrow(length = unit(0.2,"cm")))+
  theme(title = element_text(size = 15,vjust = 0.5,hjust = 0.5),
        axis.text = element_text(size = 15),
        text = element_text(size = 15),
        axis.title = element_text(size = 15),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 15))

ggsave(paste(path,"/",pre,"_bacteria_location_plot.png",sep=""),plot=p2,width = 7,height = 7,dpi = 300)
}
'''
% (options.exclude,options.input,options.meta, options.group,options.output,options.group,options.number,options.group,options.group,options.group,options.group,options.group,options.group),
file = rscript)

rscript.close()

#os.system('/System/Pipline/DNA/DNA_Micro/16S_pipeline/16S_pipeline_V1.10/software/R-3.1.0/bin/Rscript rda.R')
os.system('Rscript rda.R')
os.remove('rda.R')
