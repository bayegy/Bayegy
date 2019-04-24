library(optparse)
#######arguments
option_list <- list( 
  make_option(c("-i", "--input"),metavar="path", dest="otu",help="Specify the path of collapsed bacteria table. Required",default=NULL),
  make_option(c("-m", "--map"),metavar="path",dest="map", help="Specify the path of mapping file. Required",default=NULL),
  make_option(c("-c", "--category"),metavar="string",dest="group", help="Category to compare. Required",default=NULL),
  make_option(c("-p", "--prefix"),metavar="str", dest="prefix",help="The prefix of output files, default if null",default=""),
  make_option(c("-n", "--number"),metavar="int", dest="num",help="The number of species needed to be plotted, default is 20",default=20),
  make_option(c("-e", "--exclude"),metavar="string",dest="ex", help="Specify the numeric variables excluded from plot and seprated by commas in mapping file",default="none"),
  make_option(c("-o", "--output"),metavar="directory",dest="out", help="Specify the directory of output files. default=./",default="./")
)

opt <- parse_args(OptionParser(option_list=option_list,description = "This script is used to compare the predicted pathway of level 3"))

if(!dir.exists(opt$out)){dir.create(opt$out,recursive = T)}
opt$out<-paste(opt$out,"/",opt$prefix,sep="")


library(vegan)
library(stringr)
library(ggrepel)
library(ggplot2)
library(RColorBrewer)
library(getopt)

base_dir<-normalizePath(dirname(get_Rscript_filename()))
source(paste(base_dir,"/piputils/get_colors.R", sep = ""))
groups_color<-get_colors(opt$group, opt$map)


ex<-str_split(opt$ex,",")[[1]]
dat <- read.table(opt$otu, header = TRUE, sep = "\\t",comment.char = "",check.names = F)
# dat[,2:(ncol(dat)-1)]=apply(dat[,2:(ncol(dat)-1)],2,function(x){x/sum(x)})

dat<-dat[!duplicated(dat[,1]),]

rownames(dat)=dat[,1]
map<-read.table(opt$map,header = T,na.strings="",sep = "\\t",row.names=1,comment.char = "",check.names = F,stringsAsFactors = F)

colnames(map)[is.na(colnames(map))]<-"NA"


notstr=c()
for(i in 1:length(map)){
  notstr[i]=is.numeric(map[,i])
}

if(ex[1]!="none"){
  notstr=notstr&(!colnames(map)%in%ex)
}
group_name=opt$group
sel=c(group_name,colnames(map)[notstr])

data=map[sel]
data=na.omit(data)

dat=dat[,-c(1,ncol(dat))]
dat=t(dat)[match(rownames(data),colnames(dat)),]


groups<-data[,1]
envdata<-data[sel[-1]]

dca <- decorana(veg = dat)

dcam <- max(dca$rproj)

if (dcam > 4){
  cca <- cca(formula=dat~.,data=envdata, scale = TRUE, na.action = na.exclude)
  pre <- "CCA"
}else{

  cca <- rda(formula=dat~.,data=envdata,scale = TRUE, na.action = na.exclude)
  pre <- "RDA"
}


rda_for_mcpp_group<-rda(dat,factor(groups),scale = TRUE, na.action = na.exclude)
mcpp_group<-anova(rda_for_mcpp_group,step=1000,perm.max=1000)
title_group<-paste("Permutation test",": P=",mcpp_group[[4]][1],sep="")
print(mcpp_group)
print(title_group)

mcpp_numeric<-anova(cca,step=1000,perm.max=1000)
title_all<-paste("Overall permutation test: P=",mcpp_numeric[[4]][1],sep="")
print(title_all)


path=opt$out
ccascore <- scores(cca)
write.table(ccascore$sites, file = paste(path,group_name,"_", pre, ".sample.txt", sep = ""), sep = "\\t")
write.table(ccascore$species, file = paste(path,group_name,"_", pre, ".bacteria.txt", sep = ""), sep = "\\t")
envfit <- envfit(cca, envdata, permu = 2000, na.rm = TRUE)
rp <- cbind(as.matrix(envfit$vectors$r), as.matrix(envfit$vectors$pvals))
colnames(rp) <- c("r2", "Pr(>r)")
env <- cbind(envfit$vectors$arrows, rp)
write.table(as.data.frame(env), file = paste(path, group_name, "_",pre, ".envfit.txt", sep = ""),sep = "\\t")


new<-cca$CCA

# samples=data.frame(new$u)
samples=data.frame(ccascore$sites)
range_sam<-max(samples)-min(samples)
# sum(rownames(samples)!=rownames(data))==0
samples$id=rownames(samples)
samples=data.frame(samples,Group=groups)



# species=data.frame(new$v)
species=data.frame(ccascore$species)
range_spe<-max(species)-min(species)
species$id=rownames(species)
species<-species[match(colnames(dat),rownames(species)),]
sum<-scale(colSums(dat))[,1]
species$abundance<-sum-min(sum)+0.8
number<-as.integer(opt$num)+1
if (dim(dat)[2]<number){
  cut=-1
}else{
  cut=sort(colSums(dat),T)[number]
}
show_species<-species[colSums(dat)>cut,]
hide_species<-species[colSums(dat)<=cut,]
print(paste("The threshold of abundance is:",cut))


# envis=data.frame(envfit$vectors$arrows)#Envis seem to be the same length after permutation
# envis=data.frame(new$biplot)
envis=data.frame(env)
envis[,1]<-envis[,1]*envis[,3]
envis[,2]<-envis[,2]*envis[,3]
range_env<-max(envis[,c(1,2)])-min(envis[,c(1,2)])
envis$id=rownames(envis)

pc1 = cca$CCA$eig[1]/sum(cca$CCA$eig) * 100
pc2 = cca$CCA$eig[2]/sum(cca$CCA$eig) * 100
xlab <- paste(pre, "1: ", round(pc1, digits = 2), "%", sep = "")
ylab <- paste(pre, "2: ", round(pc2, digits = 2), "%", sep = "")



ratio_spe_env<-range_spe/range_env
ratio_sam_env<-range_sam/range_env



if(pre=="RDA"){
envis[,1]<-envis[,1]*ratio_sam_env*0.85
envis[,2]<-envis[,2]*ratio_sam_env*0.85

p1<-ggplot(data=samples,aes(x=RDA1,y=RDA2)) +
  geom_point(aes(x=RDA1,y=RDA2,color=Group),size=3) +
  # geom_text_repel(aes(x=RDA1,y=RDA2,label=id),color="black",size=3)+
  geom_text_repel(data=envis,aes(x=RDA1,y=RDA2,label=id),color="black",size=5) +
  geom_hline(yintercept=0,linetype="dotted") + geom_vline(xintercept=0,linetype="dotted")+
  theme_bw() + theme(panel.grid=element_blank())+xlab(xlab)+ylab(ylab)+
  labs(title=title_group)+
  geom_segment(aes(x=0,y=0,xend = RDA1, yend = RDA2),data = envis,
                color=brewer.pal(8,"Accent")[6],size=0.4,
                arrow = arrow(length = unit(0.2,"cm")))+
  theme(title = element_text(size = 15,vjust = 0.5,hjust = 0.5),
        axis.text = element_text(size = 15),
        text = element_text(size = 15),
        axis.title = element_text(size = 15),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 15))+
  scale_colour_manual(values = groups_color)

ggsave(paste(path,group_name,"_",pre,"_sample_location_plot.pdf",sep=""),plot=p1,width = 8,height = 7,dpi = 300)


p2<-p1+geom_text_repel(aes(x=RDA1,y=RDA2,label=id),color="black",size=3)

ggsave(paste(path,group_name,"_",pre,"_sample_location_plot_with_labels.pdf",sep=""),plot=p2,width = 8,height = 7,dpi = 300)


envis[,1]<-envis[,1]*(1/ratio_sam_env)*ratio_spe_env
envis[,2]<-envis[,2]*(1/ratio_sam_env)*ratio_spe_env

p2<-ggplot(data=show_species,aes(x=RDA1,y=RDA2)) +
  geom_point(aes(x=RDA1,y=RDA2,size=abundance),color=rainbow(3)[1]) +
  geom_point(data=hide_species,aes(x=RDA1,y=RDA2,size=abundance),color="grey") +
  guides(size=F)+
  geom_text_repel(data=envis,aes(x=RDA1,y=RDA2,label=id),color="black",size=5) +
  geom_text_repel(aes(RDA1,RDA2,label=id)) +
  geom_hline(yintercept=0,linetype="dotted") + geom_vline(xintercept=0,linetype="dotted")+
  theme_bw() + theme(panel.grid=element_blank())+xlab(xlab)+ylab(ylab)+
  ggtitle(title_all)+
  geom_segment(aes(x=0,y=0,xend = RDA1, yend = RDA2),data = envis,
                color=brewer.pal(8,"Accent")[6],size=0.4,
                arrow = arrow(length = unit(0.2,"cm")))+
  theme(title = element_text(size = 15,vjust = 0.5,hjust = 0.5),
        axis.text = element_text(size = 15),
        text = element_text(size = 15),
        axis.title = element_text(size = 15),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 15))

# ggsave(paste(path,"%s_",pre,"_bacteria_location_plot.png",sep=""),plot=p2,width = 7,height = 7,dpi = 300)
ggsave(paste(path,group_name,"_",pre,"_bacteria_location_plot.pdf",sep=""),plot=p2,width = 7,height = 7,dpi = 300)
}else{
envis[,1]<-envis[,1]*ratio_sam_env*0.85
envis[,2]<-envis[,2]*ratio_sam_env*0.85

# envis[,1]<-envis[,1]*3.5
# envis[,2]<-envis[,2]*3.5
p1<-ggplot(data=samples,aes(x=CCA1,y=CCA2)) +
  geom_point(aes(x=CCA1,y=CCA2,color=Group),size=3) +
  # geom_text_repel(aes(x=CCA1,y=CCA2,label=id),color="black",size=3)+
  geom_text_repel(data=envis,aes(x=CCA1,y=CCA2,label=id),color="black",size=5) +
  geom_hline(yintercept=0,linetype="dotted") + geom_vline(xintercept=0,linetype="dotted")+
  theme_bw() + theme(panel.grid=element_blank())+xlab(xlab)+ylab(ylab)+
  labs(title=title_group)+
  geom_segment(aes(x=0,y=0,xend = CCA1, yend = CCA2),data = envis,
                color=brewer.pal(8,"Accent")[6],size=0.4,
                arrow = arrow(length = unit(0.2,"cm")))+
  theme(title = element_text(size = 15,vjust = 0.5,hjust = 0.5),
        axis.text = element_text(size = 15),
        text = element_text(size = 15),
        axis.title = element_text(size = 15),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 15))+
  scale_colour_manual(values = groups_color)

ggsave(paste(path,group_name,"_",pre,"_sample_location_plot.pdf",sep=""),plot=p1,width = 8,height = 7,dpi = 300)

p2<-p1+geom_text_repel(aes(x=CCA1,y=CCA2,label=id),color="black",size=3)

ggsave(paste(path,group_name,"_",pre,"_sample_location_plot_with_labels.pdf",sep=""),plot=p2,width = 8,height = 7,dpi = 300)



envis[,1]<-envis[,1]*(1/ratio_sam_env)*ratio_spe_env
envis[,2]<-envis[,2]*(1/ratio_sam_env)*ratio_spe_env

p2<-ggplot(data=show_species,aes(x=CCA1,y=CCA2)) +
  geom_point(aes(x=CCA1,y=CCA2,size=abundance),color=rainbow(3)[1]) +
  geom_point(data=hide_species,aes(x=CCA1,y=CCA2,size=abundance),color="grey") +
  guides(size=F)+
  geom_text_repel(data=envis,aes(x=CCA1,y=CCA2,label=id),color="black",size=5) +
  geom_text_repel(aes(CCA1,CCA2,label=id)) +
  geom_hline(yintercept=0,linetype="dotted") + geom_vline(xintercept=0,linetype="dotted")+
  theme_bw() + theme(panel.grid=element_blank())+xlab(xlab)+ylab(ylab)+
  ggtitle(title_all)+
  geom_segment(aes(x=0,y=0,xend = CCA1, yend = CCA2),data = envis,
                color=brewer.pal(8,"Accent")[6],size=0.4,
                arrow = arrow(length = unit(0.2,"cm")))+
  theme(title = element_text(size = 15,vjust = 0.5,hjust = 0.5),
        axis.text = element_text(size = 15),
        text = element_text(size = 15),
        axis.title = element_text(size = 15),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 15))

# ggsave(paste(path,"%s_",pre,"_bacteria_location_plot.png",sep=""),plot=p2,width = 7,height = 7,dpi = 300)
ggsave(paste(path,group_name, "_",pre,"_bacteria_location_plot.pdf",sep=""),plot=p2,width = 7,height = 7,dpi = 300)
}
