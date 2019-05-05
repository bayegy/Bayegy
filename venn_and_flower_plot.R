
library(optparse)

option_list <- list( 
    make_option(c("-i", "--input"),metavar="path", dest="otu",help="Specify the path of otu table with taxonomy at last column",default=NULL),
    make_option(c("-m", "--map"),metavar="path",dest="map", help="Specify the path of mapping file",default=NULL),
    make_option(c("-c", "--category"),metavar="string",dest="group", help="Specify category name in mapping file",default="none"),
    make_option(c("-t", "--threshold"),metavar="int or float", dest="thresh",help="The threshold of abundance for the judgement of existence, default is 0",default=0),
    make_option(c("-s", "--skip"),metavar="logical",dest="skip", help="T(Skip the first line(e.g. comment line) while reading abundance table) or F(not skip first line)",default=TRUE),
    make_option(c("-o", "--output"),metavar="directory",dest="out", help="Specify the directory of output files",default="./")
    )
opt <- parse_args(OptionParser(option_list=option_list,description = "This script is used to plot venndiagram and flower diagram, and to display the special and common otus among groups"))

library(VennDiagram)
library(plotrix)
library(getopt)
if(!dir.exists(opt$out)){dir.create(opt$out,recursive = T)}
#ag<-commandArgs(T)
#if (length(ag)<3){
#	print("Please supply the following arguments:
#	1.the toxonomic otu table with the taxonomy at the last column;
#	2.the mapping file with the SampleID at first column
#	3.the column name in mapping file of group
#	4.the path of output files
#  5.the threshold of abundance for the judgement of existence")
#}else{

base_dir<-normalizePath(dirname(get_Rscript_filename()))
source(paste(base_dir,"/piputils/get_colors.R", sep = ""))
groups_color<-get_colors(opt$group, opt$map)


otu<-opt$otu
map<-opt$map
group<-opt$group
out<-opt$out
thresh<-opt$thresh

meta<-read.table(map,na.strings="",row.names=1,header = T,sep = "\t",comment.char = "",check.names = F,stringsAsFactors = F)

meta<-meta[group]
sk=ifelse(opt$skip,1,0)
data<-read.table(otu,header = T,skip = sk,sep = "\t",
                 comment.char = "",stringsAsFactors = F,check.names = F,row.names = 1)



taxonomy<-data[,dim(data)[2]]
data<-data[,-dim(data)[2]]
meta<-meta[match(colnames(data),rownames(meta)),]
data<-data[,!is.na(meta)]
meta<-meta[!is.na(meta)]

#print(rownames(meta))

abdc<-colSums(t(data))
func1<-function(x){
  return(tapply(x,INDEX = meta,sum)>thresh)
}

tb1<-apply(data,1,func1)
unig<-rownames(tb1)
d1<-colSums(tb1)
ll<-length(unig)

com<-sum(d1==ll)
com_spe<-data.frame(OTUID=colnames(tb1)[d1==ll],Taxonomy=taxonomy[d1==ll],SUM_Abundance=abdc[d1==ll])
write.table(com_spe,paste(out,"/",group,"_common_species.xls",sep = ""),row.names = F,sep = "\t")


flower_data<-c()
for(i in 1:ll){
  is_uni<-(d1==1&tb1[i,]==1)
  flower_data[i]<-sum(is_uni)
	uni_spe<-data.frame(OTUID=colnames(tb1)[is_uni],Taxonomy=taxonomy[is_uni],SUM_Abundance=abdc[is_uni])
  write.table(uni_spe,paste(out,"/",group,"_",unig[i],"_special_species.xls",sep = ""),row.names = F,sep = "\t")
}

if(length(unig)<=5){
ls1<-list()
for(i in 1:ll){ls1[[i]]<-colnames(tb1)[tb1[i,]]}
names(ls1)<-rownames(tb1)
venn.diagram(ls1,filename = paste(out,"/",group,"_Venn_plot.png",sep = ""),imagetype="png",alpha= 0.50,lwd =1.2,cat.cex=1.4,fill=groups_color,margin=0.15)

}


if(length(unig)>=3){
flower_plot <- function(sample, value, common, start, a, b,circ_r=1.5,ell_pos=2,tax_pos=4,num_pos=2, 
                        ellipse_col, 
                        circle_col = rgb(0, 162, 214, max = 255),
                        text_cex = 1) {
  par( bty = "n", ann = F, xaxt = "n", yaxt = "n", mar = c(1,1,1,1))
  plot(c(0,14),c(0,14),type="n")
  n   <- length(sample)
  deg <- 360 / n
  res <- lapply(1:n, function(t){
    draw.ellipse(x = 7 + ell_pos*cos((start + deg * (t - 1)) * pi / 180), 
                 y = 7 + ell_pos*sin((start + deg * (t - 1)) * pi / 180), 
                 col = ellipse_col[t],
                 border = ellipse_col[t],
                 a = a, b = b, angle = deg * (t - 1))
    text(x = 7+ num_pos * cos((start + deg * (t - 1)) * pi / 180),
         y = 7 + num_pos * sin((start + deg * (t - 1)) * pi / 180),
         value[t],cex = text_cex
    )
    
    if (deg * (t - 1) < 180 && deg * (t - 1) > 0 ) {
      text(x = 7 + tax_pos * cos((start + deg * (t - 1)) * pi / 180),
           y = 7 + tax_pos * sin((start + deg * (t - 1)) * pi / 180),
           sample[t],
           #srt = deg * (t - 1) - start,
           adj = 1,
           cex = text_cex
      )
      
    } else {
      text(x = 7 + tax_pos * cos((start + deg * (t - 1)) * pi / 180),
           y = 7 + tax_pos * sin((start + deg * (t - 1)) * pi / 180),
           sample[t],
           #srt = deg * (t - 1) + start,
           adj = 0,
           cex = text_cex
      )
    }     
  })
  draw.circle(x = 7, y = 7, r = circ_r, col = circle_col, border = circle_col)
  text(x=7,y=7,common,cex = text_cex)
}




pdf(paste(out,"/",group,"_Flower_plot.pdf",sep = ""),width = 6,height = 6)
flower_plot(sample=unig,value = flower_data,
            a=0.65,b=1.8,start =90,common = com,tax_pos = 4.3,num_pos = 2.5,
            ellipse_col = groups_color,
            ell_pos = 2.5,circ_r = 0.9)
dev.off()
}

#}
