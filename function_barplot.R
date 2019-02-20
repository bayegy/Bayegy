#utf-8
library(optparse)
#######arguments
option_list <- list( 
  make_option(c("-i", "--input"),metavar="path", dest="func",help="Specify the path of predicted KEGG pathway file of level 3. Required",default=NULL),
  make_option(c("-m", "--map"),metavar="path",dest="map", help="Specify the path of mapping file. Required",default=NULL),
  make_option(c("-c", "--category"),metavar="string",dest="group", help="Category to compare. Required",default=NULL),
  make_option(c("-l", "--log"),metavar = "logical",dest="log", help="If TRUE, log the data before comparing. Options of logical type are TRUE, T, FALSE, F. default = FALSE",default = "FALSE"),
  make_option(c("-j", "--adjust-p"),metavar = "logical",dest="ap", help="If TRUE, adjust the p value before barplot. default = TRUE",default = "TRUE"),
  make_option(c("-e", "--add-se"),metavar = "logical",dest="se", help="If TRUE, add SE error bar, otherwise add SD error bar. default = TRUE",default = "TRUE"),
  make_option(c("-a", "--alpha"),metavar = "float",dest="alpha", help="Alpha for significance. default=0.05",default=0.05),
  make_option(c("-b", "--output-by-L1"),metavar = "logical",dest="bl", help="IF TRUE, output barplot for each pathway of L1 level. If FALSE, output only one plot anyway. Or use auto to determine by the number of siginficant function (10)",default="auto"),
  make_option(c("-s", "--scale"),metavar = "logical",dest="scale", help="IF TRUE, scale the data before plot.",default="F"),
  make_option(c("-o", "--output"),metavar="directory",dest="out", help="Specify the directory of output files. default=./",default="./")
)

opt <- parse_args(OptionParser(option_list=option_list,description = "This script is used to compare the predicted pathway of level 3"))

library(stringr)
library(ggpubr)
library(reshape2)
library(agricolae)
library(ggplot2)

if(!dir.exists(opt$out)){dir.create(opt$out,recursive = T)}

func<-read.table(opt$func,comment.char="",quote = "",skip = 1,check.names=F,stringsAsFactors=F, header = TRUE, sep = "\t")
map<-read.table(opt$map,sep="\t",na.strings="",header = T,row.names=1,comment.char = "",check.names = F,stringsAsFactors = F)
group<-map[opt$group]
group<-na.omit(group)
#clean na of group

rownames(func)<-func[,ncol(func)]
func<-func[,-c(1,ncol(func))]
#清理数据
func<-func[,match(rownames(group),colnames(func))]
#clean na of func
func<-apply(func,2,function(x){x/sum(x)})
#求相对丰度
func<-t(func)
func<-func[,colSums(func>0)>=nrow(func)*0.5]
#剔除观测样本数小于50%总样本数的功能
if(as.logical(opt$log)){
  func<-func-min(func)
  min2<-min(func[func!=0])
  func<-log(func+min2,base = 10)-log(min2,base = 10)
}
#log转化，使得尽量符合正太分布
if(as.logical(opt$scale)){
  func<-scale(func,center=T,scale=T)
}
#normalize the data, make the plot easier to compare

func<-data.frame(t(func),check.names = F)
group<-group[match(colnames(func),rownames(group)),]
N_sample<-ncol(func)
sprintf("Proccesing anova of %d samples",N_sample)
data_for_anova<-melt(data.frame(group=group,t(func),check.names = FALSE),id.vars = "group")

anova_results<-compare_means(value~group,data = data_for_anova,method = "anova",group.by = "variable")
write.table(anova_results,paste(opt$out,'/',opt$group,"_all_pathway_anova_results.xls",sep = ""),sep = "\t",row.names = F)

if(as.logical(opt$ap)){
  pvalue<-anova_results$p.adj
}else{
  pvalue<-anova_results$p
}

siginficant_function<-as.character(anova_results$variable[pvalue<as.numeric(opt$alpha)])

uni_group<-sort(unique(group))
N_group<-length(uni_group)
my_duncan<-function(x){
  data_for_duncan<-data.frame(group=group,value=x)
  anova <- aov(value~group,data = data_for_duncan)
  plotdata<-duncan.test(anova,"group", alpha = 0.05)
  out<-plotdata$groups
  out<-out[match(uni_group,rownames(out)),]
  out1<-as.character(out[,2])
  names(out1)<-rownames(out)
  return(out1)
}



if(opt$bl=="auto"){
  lsig<-length(siginficant_function)
  if(lsig>10){
    how<-1
  }else{
    how<-2
  }
}else{
  if(as.logical(opt$bl)){
    how<-1
  }else{
    how<-2
  }
}


if(how==1){
  LF<-str_extract(rownames(func),'^[^;]+')
  unilf<-unique(LF)
  for (l1 in unilf){
    if(sum(LF==l1&rownames(func)%in%siginficant_function)>=1){
      sub_func<-func[LF==l1&rownames(func)%in%siginficant_function,]
      res<-apply(sub_func,1,my_duncan)
      LR<-str_extract(colnames(res),"[^;]+$")
      N_res<-ncol(res)
      sum_mean<-apply(sub_func,1,function(x){tapply(x,INDEX = group,mean)})
      sum_se<-apply(sub_func,1,function(x){tapply(x,INDEX = group,ifelse(as.logical(opt$se),function(x){sd(x)/sqrt(N_sample)},sd))})
      sum_mean<-sum_mean[match(uni_group,rownames(sum_mean)),]
      res<-res[match(uni_group,rownames(res)),]
      sum_se<-sum_se[match(uni_group,rownames(sum_se)),]
      bar_data<-data.frame(matrix(nrow = N_group*N_res,ncol = 5))
      colnames(bar_data)<-c("pathway","group","mean","se","sig")
      bar_data$pathway<-rep(LR,each=N_group)
      bar_data$group<-rep(uni_group,times=N_res)
      bar_data$mean<-as.vector(sum_mean)
      bar_data$se<-as.vector(sum_se)
      bar_data$sig<-as.vector(res)
      bar_data<-bar_data[order(bar_data$pathway),]
#      print(bar_data$group)
      p1<-ifelse(0.2*N_group*N_res+2<50,0.2*N_group*N_res+2,49.9)
      p2<-max(bar_data$mean)/18
      cud<-0.8/(2*N_group)
      pos<-rep(seq(from=0,to=N_res-1),each=N_group)
      cood_x<-rep(seq(0.6+cud,1.4-cud,length.out = N_group),times=N_res)+pos
      p<-ggplot(bar_data, aes(x=pathway, y=mean, fill=group)) + 
        geom_bar(width = 0.8,position=position_dodge(0.8), stat="identity") +
        geom_errorbar(aes(ymin=mean-se, ymax=mean+se),width=0.2,position=position_dodge(0.8))+
        geom_text(aes(x=cood_x,y=mean+se+p2,label=sig))+
        guides(fill=guide_legend(title = NULL))+
        xlab("")+ylab("Abundance")+theme_bw()+
        theme(text = element_text(size = 12),
              panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
              axis.line = element_line(),panel.border =  element_blank(),
              axis.text.x = element_text(angle = 90,size = 9,vjust = 1,hjust = 1))+scale_y_continuous(expand = c(0, 0))
      ggsave(plot = p,paste(opt$out,'/',opt$group,"_",l1,"_barplot_of_duncan.pdf",sep = ""),dpi=300,height = 7,width = p1)  
      if(as.logical(opt$se)){
        colnames(bar_data)<-c("KEGG pathway","Group","Mean","SE","Duncan significance")
      }else{
        colnames(bar_data)<-c("KEGG pathway","Group","Mean","SD","Duncan significance")
      }
      write.table(bar_data,paste(opt$out,'/',opt$group,"_",l1,"_duncan_results.xls",sep = ""),sep = "\t",row.names = F)
    }
  }
}else{
  if(sum(rownames(func)%in%siginficant_function)>=1){
    sub_func<-func[rownames(func)%in%siginficant_function,]
    res<-apply(sub_func,1,my_duncan)
    LR<-str_extract(colnames(res),"[^;]+$")
    N_res<-ncol(res)
    sum_mean<-apply(sub_func,1,function(x){tapply(x,INDEX = group,mean)})
    sum_se<-apply(sub_func,1,function(x){tapply(x,INDEX = group,ifelse(as.logical(opt$se),function(x){sd(x)/sqrt(N_sample)},sd))})
    sum_mean<-sum_mean[match(uni_group,rownames(sum_mean)),]
    res<-res[match(uni_group,rownames(res)),]
    sum_se<-sum_se[match(uni_group,rownames(sum_se)),]
    bar_data<-data.frame(matrix(nrow = N_group*N_res,ncol = 5))
    colnames(bar_data)<-c("pathway","group","mean","se","sig")
    bar_data$pathway<-rep(LR,each=N_group)
    bar_data$group<-rep(uni_group,times=N_res)
    bar_data$mean<-as.vector(sum_mean)
    bar_data$se<-as.vector(sum_se)
    bar_data$sig<-as.vector(res)
    bar_data<-bar_data[order(bar_data$pathway),]
    p1<-ifelse(0.2*N_group*N_res+2<50,0.2*N_group*N_res+2,49.9)
    p2<-max(bar_data$mean)/18
    cud<-0.8/(2*N_group)
    pos<-rep(seq(from=0,to=N_res-1),each=N_group)
    cood_x<-rep(seq(0.6+cud,1.4-cud,length.out = N_group),times=N_res)+pos
    p<-ggplot(bar_data, aes(x=pathway, y=mean, fill=group)) + 
      geom_bar(width = 0.8,position=position_dodge(0.8), stat="identity") +
      geom_errorbar(aes(ymin=mean-se, ymax=mean+se),width=0.2,position=position_dodge(0.8))+
      geom_text(aes(x=cood_x,y=mean+se+p2,label=sig))+
      guides(fill=guide_legend(title = NULL))+
      xlab("")+ylab("Abundance")+theme_bw()+
      theme(text = element_text(size = 12),
            panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
            axis.line = element_line(),panel.border =  element_blank(),
            axis.text.x = element_text(angle = 90,size = 9,vjust = 1,hjust = 1))+scale_y_continuous(expand = c(0, 0))
    ggsave(plot = p,paste(opt$out,'/',opt$group,"_all_significant_pathway_barplot_of_duncan.pdf",sep = ""),dpi=300,height = 7,width = p1)  
    if(as.logical(opt$se)){
      colnames(bar_data)<-c("KEGG pathway","Group","Mean","SE","Duncan significance")
    }else{
      colnames(bar_data)<-c("KEGG pathway","Group","Mean","SD","Duncan significance")
    }
    write.table(bar_data,paste(opt$out,'/',opt$group,"_all_significant_pathway_duncan_results.xls",sep = ""),sep = "\t",row.names = F)
  }
}


