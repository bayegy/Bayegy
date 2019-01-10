library(optparse)

#######arguments
option_list <- list( 
    make_option(c("-i", "--input"),metavar="path", dest="otu",help="OTU table with Consensus Lineage at last column",default=NULL),
    make_option(c("-m", "--map"),metavar="path",dest="map", help="Specify the path of mapping file. Required",default=NULL),
    make_option(c("-t", "--tree"),metavar="path",dest="tree", help="rooted tree. Required",default=NULL),
    make_option(c("-r", "--rep-seqs"),metavar="path",dest="rep", help="Representative sequences. Required",default=NULL),
    make_option(c("-c", "--category"),metavar="string",dest="group", help="Category to compare. Required",default=NULL),
    make_option(c("-p", "--prefix"),metavar="str", dest="prefix",help="The prefix of output files, default if null",default=""),
    make_option(c("-o", "--output"),metavar="directory",dest="out", help="Specify the directory of output files",default="./")
    )

opt <- parse_args(OptionParser(option_list=option_list,description = "This script is used to draw the relative stack barplot of species"))



if(!dir.exists(opt$out)){dir.create(opt$out,recursive = T)}
opt$out<-paste(opt$out,"/",opt$prefix,sep="")

library("phyloseq")
library("pheatmap")

cluster<-ifelse(is.null(opt$map)|is.null(opt$group),TRUE,FALSE)

print(cluster)
if(!cluster){
    map<-read.table(opt$map,quote="",row.names = 1,na.strings = "",comment.char="",check.names=F,stringsAsFactors=F, header = TRUE, sep = "\t")
    group<-na.omit(map[c(opt$group,"Description")])
    group<-group[order(rownames(group)),]
    group<-group[order(group[,1]),]
}


qiimedata = import_qiime(opt$otu, opt$map, opt$tree, opt$rep)

print("#calculate distance")
for (distance_matrix in list(c('bray','bray_curtis'), c('unifrac','unweighted_unifrac'), c('wunifrac','weighted_unifrac'))){
  beta_heatmap_outputpdfname <- paste(opt$out,distance_matrix[2], "_betadiversity_summary.pdf", sep="")
  #pdf(beta_heatmap_outputpdfname)
  Dist <- distance(qiimedata, method=distance_matrix[1])

  beta_outputtxtname <- paste(opt$out,distance_matrix[2], "_matrix.txt", sep="")
  write.table(as.matrix(Dist), beta_outputtxtname , quote=FALSE, col.names=NA, sep="\t")
  Dist_read<-read.table(beta_outputtxtname, head=T)
  if(cluster){
    pdf(beta_heatmap_outputpdfname,height = nrow(Dist_read)*0.3+3,width = ncol(Dist_read)*0.3+4)
    pheatmap(Dist_read,fontsize=10,border_color = "black",fontsize_row =10,
             cluster_rows=T,clustering_distance_rows="euclidean",
             cluster_cols=T,clustering_distance_cols="euclidean")
    dev.off()
  }else{
    print("no cluster heatmap")
    Dist_read<-data.frame(Dist_read)
    Dist_read<-Dist_read[match(rownames(group),rownames(Dist_read)),]
    Dist_read<-Dist_read[,match(rownames(group),colnames(Dist_read))]
    pdf(beta_heatmap_outputpdfname,height = nrow(Dist_read)*0.3+3,width = ncol(Dist_read)*0.3+4)
    pheatmap(Dist_read,fontsize=10,border_color = "black",fontsize_row =10,cluster_rows=F,cluster_cols=F)
    dev.off()
  }
}
