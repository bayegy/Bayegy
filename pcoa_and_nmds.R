##########Library import

library(optparse)

option_list <- list( 
    make_option(c("-i", "--input"),metavar="path", dest="otu",help="Specify the path of otu table with taxonomy at last column",default=NULL),
    make_option(c("-m", "--map"),metavar="path",dest="map", help="Specify the path of mapping file",default=NULL),
    make_option(c("-c", "--category"),metavar="string",dest="group", help="Specify category name in mapping file",default="none"),
    make_option(c("-o", "--output"),metavar="directory",dest="out", help="Specify the directory of output files",default="./")
    )
opt <- parse_args(OptionParser(option_list=option_list,description = "This script is used to plot venndiagram and flower diagram, and to display the special and common otus among groups"))


if(!dir.exists(opt$out)){dir.create(opt$out,recursive = T)}
library("phyloseq")
library("scatterplot3d")
library("ggplot2")
library("ggrepel")

args <- commandArgs(trailingOnly = TRUE)

map = opt$map
category1=opt$group
otu=opt$otu



########################Using Phyloseq
gpt = import_qiime(otufilename=otu,mapfilename=map)


otu<-gpt@otu_table@.Data
sum_of_otus<-colSums(t(otu))
selected_otu<-names(sum_of_otus)[sum_of_otus>10]
gpt <- prune_taxa(selected_otu, gpt)


print("#Generate the NMDS plot for betadiversity")
for (distance_matrix in list(c('bray','bray_curtis'))){
  GP.ord <- ordinate(gpt, "NMDS", distance_matrix[1])
  NMDS_outputpdfname <- paste(opt$out,"/",category1,"_",distance_matrix[2], "_NMDS.pdf", sep="")
  NMDS_ordtxtname <- paste(opt$out,"/",category1,"_",distance_matrix[2], "_NMDS.ord.txt", sep="")
  pdf(NMDS_outputpdfname, width=7.6, height=6.6)
  p2 = plot_ordination(gpt, GP.ord, type="samples", color=category1) 
  p3 = p2  + geom_point(size=3) + geom_text_repel(aes(label=Description),hjust=0, vjust=2, size=4) + stat_ellipse()+theme(text = element_text(size = 15))
  print(p3 + ggtitle(distance_matrix[2]))
  dev.off()

  ####without names and ellipse
  pdf(paste(opt$out,"/",category1,"_",distance_matrix[2], "_NMDS_without_labels.pdf", sep=""), width=7.6, height=6.6)
  p2 = plot_ordination(gpt, GP.ord, type="samples", color=category1) 
  p3 = p2  + geom_point(size=3) +theme(text = element_text(size = 15))
  print(p3 + ggtitle(distance_matrix[2]))
  dev.off()


  write.table(as.matrix(GP.ord), NMDS_ordtxtname, quote=FALSE, col.names=NA, sep="\t")
}

print("#Generate the PCoA 2D plot for betadiversity")
for (distance_matrix in list(c('bray','bray_curtis'))){
  GP.ord <- ordinate(gpt, "PCoA", distance_matrix[1])
  PCoA_outputpdfname <- paste(opt$out,"/",category1,"_",distance_matrix[2], "_PCoA_2D.pdf", sep="")
  PCoA_ordtxtname <- paste(opt$out,"/",category1,"_",distance_matrix[2], "_PCoA_2D.ord.txt", sep="")


  pdf(PCoA_outputpdfname, width=7.6, height=6.6)
  p2 = plot_ordination(gpt, GP.ord, type="samples", color=category1) 
  p3 = p2  + geom_point(size=3) + geom_text_repel(aes(label=Description),hjust=0, vjust=2, size=4) + stat_ellipse()+theme(text = element_text(size = 15))
  print(p3 + ggtitle(distance_matrix[2]))
  dev.off()

  ######without names and ellipse
  pdf(paste(opt$out,"/",category1,"_",distance_matrix[2], "_PCoA_2D_without_labels.pdf", sep=""), width=7.6, height=6.6)
  p2 = plot_ordination(gpt, GP.ord, type="samples", color=category1) 
  p3 = p2  + geom_point(size=3)+theme(text = element_text(size = 15))
  print(p3 + ggtitle(distance_matrix[2]))
  dev.off()
  print("#Generate the PCoA 3D plot for betadiversity")
  ######pcoa 3d plot
  asign_rainbow_cor<-function(x){
    unig<-unique(x)
    color_out<-x
    for (i in 1:length(unig)) {
      color_out[color_out==unig[i]]<-rainbow(length(unig))[i]
    }
    return(color_out)
  }

  gp<-as.character(data.frame(sample_data(gpt))[category1][,1])
  tdata<-GP.ord$vectors[,1:3]
  eig<-data.frame(GP.ord$values)["Eigenvalues"][,1]
  lab<-paste("PC",c(1:3)," ",round((eig[1:3]/sum(eig))*100,digits=2),"%",sep="")
  pdf(paste(opt$out,"/",category1,"_",distance_matrix[2], "_PCoA_3D.pdf", sep=""), width=8, height=6.4)
  opar<-par(no.readonly=TRUE)
  par(fig=c(0,0.75,0,1))
  scatterplot3d(tdata,mar=c(2.2,2.2,0,0)+1,xlab=lab[1],ylab=lab[2],zlab=lab[3],color=asign_rainbow_cor(gp), grid=TRUE, box=F, type="h", lty.hplot=2, pch=19)
  par(fig=c(0.75,1,0,1),xpd=TRUE)
  legend("center", legend = unique(gp),bty = 'n',xpd = TRUE,horiz = FALSE,col = rainbow(length(unique(gp))), pch = 19, inset = -0.1)
  par(opar)
  dev.off()
  write.table(as.matrix(GP.ord), PCoA_ordtxtname, quote=FALSE, col.names=NA, sep="\t")

}

