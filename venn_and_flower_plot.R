library(VennDiagram)
library(plotrix)
ag<-commandArgs(T)
if (length(ag)<3){
	print("Please supply the following arguments:
	1.the toxonomic otu table with the taxonomy at the last column;
	2.the mapping file with the SampleID at first column
	3.the column name in mapping file of group
	4.the path of output files
  5.the threshold of abundance for criterial judgement standard of existence")
}else{
otu<-ag[1]
map<-ag[2]
group<-ag[3]
out<-ag[4]
thresh<-ag[5]
meta<-read.table(map,row.names=1,header = T,sep = "\t",comment.char = "",check.names = F,stringsAsFactors = F)

meta<-meta[group]

data<-read.table(otu,header = T,skip = 1,sep = "\t",
                 comment.char = "",stringsAsFactors = F,check.names = F,row.names = 1)

taxonomy<-data[,dim(data)[2]]
data<-data[,-dim(data)[2]]

meta<-t(meta)
meta<-meta[,match(colnames(data),colnames(meta))]

func1<-function(x){
  return(tapply(x,INDEX = meta,sum)>thresh)
}

tb1<-apply(data,1,func1)
unig<-rownames(tb1)
d1<-colSums(tb1)
ll<-length(unig)

com<-sum(d1==ll)
com_spe<-data.frame(OTUID=colnames(tb1)[d1==ll],Taxonomy=taxonomy[d1==ll])
write.table(com_spe,paste(out,"/",group,"_common_species.xls",sep = ""),row.names = F,sep = "\t")


flower_data<-c()
for(i in 1:ll){
  is_uni<-(d1==1&tb1[i,]==1)
  flower_data[i]<-sum(is_uni)
	uni_spe<-data.frame(OTUID=colnames(tb1)[is_uni],Taxonomy=taxonomy[is_uni])
  write.table(uni_spe,paste(out,"/",group,"_",unig[i],"_special_species.xls",sep = ""),row.names = F,sep = "\t")
}

if(length(unig)<=5){
ls1<-list()
for(i in 1:ll){ls1[[i]]<-colnames(tb1)[tb1[i,]]}
names(ls1)<-rownames(tb1)
venn.diagram(ls1,filename = paste(out,"/",group,"_Venn_plot.pdf",sep = ""),alpha= 0.50,lwd =1.2,cat.cex=1.4,fill=rainbow(length(ls1)),margin=0.15)

}


if(length(unig)>=3){
flower_plot <- function(sample, value, common, start, a, b,circ_r=1.5,ell_pos=2,tax_pos=4,num_pos=2, 
                        ellipse_col, 
                        circle_col = rgb(0, 162, 214, max = 255),
                        text_cex = 2) {
  par( bty = "n", ann = F, xaxt = "n", yaxt = "n", mar = c(1,1,1,1))
  plot(c(0,10),c(0,10),type="n")
  n   <- length(sample)
  deg <- 360 / n
  res <- lapply(1:n, function(t){
    draw.ellipse(x = 5 + ell_pos*cos((start + deg * (t - 1)) * pi / 180), 
                 y = 5 + ell_pos*sin((start + deg * (t - 1)) * pi / 180), 
                 col = ellipse_col[t],
                 border = ellipse_col[t],
                 a = a, b = b, angle = deg * (t - 1))
    text(x = 5 + num_pos * cos((start + deg * (t - 1)) * pi / 180),
         y = 5 + num_pos * sin((start + deg * (t - 1)) * pi / 180),
         value[t],cex = text_cex
    )
    
    if (deg * (t - 1) < 180 && deg * (t - 1) > 0 ) {
      text(x = 5 + tax_pos * cos((start + deg * (t - 1)) * pi / 180),
           y = 5 + tax_pos * sin((start + deg * (t - 1)) * pi / 180),
           sample[t],
           srt = deg * (t - 1) - start,
           adj = 1,
           cex = text_cex
      )
      
    } else {
      text(x = 5 + tax_pos * cos((start + deg * (t - 1)) * pi / 180),
           y = 5 + tax_pos * sin((start + deg * (t - 1)) * pi / 180),
           sample[t],
           srt = deg * (t - 1) + start,
           adj = 0,
           cex = text_cex
      )
    }     
  })
  draw.circle(x = 5, y = 5, r = circ_r, col = circle_col, border = circle_col)
  text(x=5,y=5,common,cex = text_cex)
}




pdf(paste(out,"/",group,"_Flower_plot.pdf",sep = ""),width = 6,height = 6)
flower_plot(sample=unig,value = flower_data,
            a=0.8,b=2,start =90,common = com,tax_pos = 4.93,num_pos = 2.7,
            ellipse_col = rainbow(ll),
            ell_pos = 2.7,circ_r = 1)
dev.off()
}

}
