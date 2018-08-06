#R script for generating alpha dviersit comparison plots
ag<-commandArgs(T)


if(length(ag)<4){
print("Please input: (1) Full path of mapfile; 
      (2) Column name in mapfile you want to analyse; 
      (3) Full path of merged α diversity file; 
      (4) Path of the output files")
}else{

require(reshape)
require(ggplot2)
require(ggpubr)
library(dplyr)
library(ggsignif)

#map = "~/Desktop/WST/16S-pipeline/sample-metadata.tsv"
#group = "Group1"
#alpha= "~/Desktop/WST/16S-pipeline/alpha/alpha-summary.tsv"
#output="~/Desktop/plot.png"

  
map<-read.table(ag[1],header = T,row.names = 1,check.names = F,stringsAsFactors = F,sep = "\t",comment.char = "")
adiv<-read.table(ag[3],header = T,row.names = 1,check.names = F,stringsAsFactors = F,sep = "\t",comment.char = "")
  

a<-colnames(adiv)

map<-map[ag[2]]
colnames(map)<-"Group"
rowname_join<-function(x,y)
{
  ya<-data.frame(y[match(rownames(x),rownames(y)),])
  colnames(ya)<-colnames(y)
  colnames(ya)[colnames(ya)%in%colnames(x)]<-paste(colnames(ya)[colnames(ya)%in%colnames(x)],"_y",sep = "")
  out<-data.frame(x,ya,check.rows = T,check.names = F)
  return(out)
}

joinedtab<-rowname_join(map,adiv)
meltab<-melt(joinedtab,id.vars = "Group")
anno_dfa = compare_means(value ~ Group, group.by = "variable", data = meltab) %>% mutate(p.adj = format.pval(p.adj, digits = 2))
write.table(as.matrix(anno_dfa), paste(ag[4],"/","alpha_",ag[2],"_wilcox_compare_results.txt",sep=""), quote=FALSE, col.names=NA, sep="\t")


for(i in a){

unig<-unique(joinedtab["Group"])
p1<-dim(unig)[1]
unitab<-joinedtab[,colnames(joinedtab)%in%c("Group",i)]


colnames(unitab)[2]<-"value"

p2<-max(unitab[,2]) 
p3<-sd(unitab[,2])


anno_df<-compare_means(value ~ Group, data = unitab)
p4<-dim(anno_df)[1]
p.value.y.coord <- rep(p2,p4)

step.increase <- (1:p4)*0.5*p3
p.value.y.coord <- p.value.y.coord + step.increase
p5=max(p.value.y.coord)
anno_df<-mutate(anno_df,y_pos=p.value.y.coord)
# P-value y coordinates

alphaplotwithsig <- ggplot(unitab, aes(x=Group, y=value)) + geom_boxplot(fill=rainbow(7)[5]) + geom_point(aes(color="black"),position=position_jitterdodge())+
 ggsignif::geom_signif(data=anno_df,manual=T,aes(xmin=group1,xmax=group2,annotations=p.signif,y_position=y_pos))+guides(color=F)+
 stat_compare_means(label.y = p5+0.5*p3,label.x.npc="center")+ylab(i)+xlab("")+theme(text=element_text(size=15,face="bold"))

ggsave(paste(ag[4],"/",i,"_",ag[2],"_wilcox_compare_boxplot.png",sep=""), plot=alphaplotwithsig, height=6, width=p1*1.5+1, dpi = 300)
}
}



