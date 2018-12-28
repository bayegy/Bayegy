
#utf-8
library("phyloseq")
library("pheatmap")
args <- commandArgs(trailingOnly = TRUE)
map = args[1]

this.dir <- getwd()
new.dir <- paste(this.dir, "/R_output", sep="")
#this.dir <- paste("~/Desktop/WST/16S-pipeline/", "", sep="")
setwd(new.dir)
txt = paste(this.dir, "/exported/feature-table.ConsensusLineage.txt", sep="")
#txt = args[3]
tre = paste(this.dir, "/exported/tree.rooted.nwk", sep="")
rs = paste(this.dir, "/exported/dna-sequences.fasta", sep="")
qiimedata = import_qiime(txt, map, tre, rs)

print("#calculate distance")
for (distance_matrix in list(c('bray','bray_curtis'), c('unifrac','unweighted_unifrac'), c('wunifrac','weighted_unifrac'))){
  beta_heatmap_outputpdfname <- paste(distance_matrix[2], "_betadiversity_summary.pdf", sep="")
  #pdf(beta_heatmap_outputpdfname)
  Dist <- distance(qiimedata, method=distance_matrix[1])

  beta_outputtxtname <- paste(distance_matrix[2], "_matrix.txt", sep="")
  write.table(as.matrix(Dist), beta_outputtxtname , quote=FALSE, col.names=NA, sep="\t")
  
  Dist_read<-read.table(beta_outputtxtname, head=T)
  pdf(beta_heatmap_outputpdfname,height = nrow(Dist_read)*0.3+3,width = ncol(Dist_read)*0.3+4)
  pheatmap(Dist_read,fontsize=10,border_color = "black",fontsize_row =10,
           cluster_rows=T,clustering_distance_rows="euclidean",
           cluster_cols=T,clustering_distance_cols="euclidean")
  dev.off()
  
}
