##########Library import
library("ape")
library("phyloseq")
library("ggplot2")
library("gplots")
library("vegan")
library("ggrepel")
library("pheatmap")
library("scatterplot3d")
library(getopt)


args <- commandArgs(trailingOnly = TRUE)

##########Data import
#txt = "/Users/chengguo/Desktop/Hengchuang/M231/exported/feature-table.ConsensusLineage.txt"
#tre = "/Users/chengguo/Desktop/Hengchuang/M231/exported/tree.rooted.nwk"
#rs = "/Users/chengguo/Desktop/Hengchuang/M231/exported/dna-sequences.fasta"
map = args[1]
category1=args[2]
ifline=as.logical(args[3])




base_dir<-normalizePath(dirname(get_Rscript_filename()))
source(paste(base_dir,"/piputils/get_colors.R", sep = ""))
groups_color<-get_colors(category1, map)

#category1 = "SampleType"
#map = "~/Desktop/Hengchuang/M122/M122_Mapping.tsv"

#category1 = "Group1"
#map = "/Users/chengguo/Desktop/Hengchuang/M231/M231_Mapping_2.tsv"

#category1 = "Group1"
#map = "~/Desktop/WST/16S-pipeline/sample-metadata.tsv"

this.dir <- getwd()
new.dir <- paste(this.dir, "/R_output", sep="")
#this.dir <- paste("~/Desktop/WST/16S-pipeline/", "", sep="")


setwd(new.dir)

txt = paste(this.dir, "/exported/feature-table.ConsensusLineage.txt", sep="")
#txt = args[3]
tre = paste(this.dir, "/exported/tree.rooted.nwk", sep="")
rs = paste(this.dir, "/exported/dna-sequences.fasta", sep="")

##mao
#txt = paste(this.dir, "/feature-table.ConsensusLineage.txt", sep="")
#tre = paste(this.dir, "/tree.rooted.nwk", sep="")
#rs = paste(this.dir, "/dna-sequences.fasta", sep="")


print("#Start reading realted files with Phyloseq")
print(map)
print(txt)
print(tre)
print(rs)


########################Using Phyloseq
gpt = import_qiime(txt, map, tre, rs)



otu<-gpt@otu_table@.Data
sum_of_otus<-colSums(t(otu))
selected_otu<-names(sum_of_otus)[sum_of_otus>10]
gpt <- prune_taxa(selected_otu, gpt)

#gpt <- subset_taxa(qiimedata, Kingdom=="Bacteria")
#gpt <- prune_taxa(names(sort(taxa_sums(gpt),TRUE)[1:30]), gpt)
#gpt <- prune_taxa(names(sort(taxa_sums(gpt),TRUE)), gpt)

#head(tax_table(gpt)[,2])
#head(tax_table(gpt)[1,])


# print("#Generate phylogenetic trees for common phylums")
# for (selected_phylum in c('Bacteroidetes','Firmicutes','Proteobacteria')){
#   print(paste("Making tree plots for", selected_phylum, sep=" "))
#   GP.chl <- subset_taxa(gpt, Phylum==selected_phylum)
#   phylogeny_outputpdfname <- paste(selected_phylum, ".phylogeny.pdf", sep="")
#   pdf( phylogeny_outputpdfname, width=12, height=14)
#   plot<-plot_tree(GP.chl, color=category1, shape="Family", label.tips="Genus", size="abundance", plot.margin=0.1, base.spacing=0.04, ladderize=TRUE, nodelabf=nodeplotblank)
#   print(plot+ ggtitle(selected_phylum))
#   dev.off()
# }

# pdf("Bacteria.phylogeny.pdf", width=12, height=14)
# plot<-plot_tree(gpt, color=category1,  label.tips="Family",  shape="Phylum", size="abundance", text.size=4, plot.margin=0.1, base.spacing=0.06, ladderize=TRUE, nodelabf=nodeplotblank)
# print(plot + ggtitle("Phylogenetics trees of Family of Bacteria"))
# dev.off()


print("#Generate the NMDS plot for betadiversity")
for (distance_matrix in list(c('bray','bray_curtis'), c('unifrac','unweighted_unifrac'), c('wunifrac','weighted_unifrac'))){
  GP.ord <- ordinate(gpt, "NMDS", distance_matrix[1])
  NMDS_outputpdfname <- paste(category1,"_",distance_matrix[2], "_NMDS.pdf", sep="")
  NMDS_ordtxtname <- paste(category1,"_",distance_matrix[2], "_NMDS.ord.txt", sep="")
  pdf(NMDS_outputpdfname, width=7.6, height=6.6)
  p2 = plot_ordination(gpt, GP.ord, type="samples", color=category1) 
  p3 = p2  + geom_point(size=3) + geom_text_repel(aes(label=Description),hjust=0, vjust=2, size=4) + stat_ellipse()+theme(text = element_text(size = 15))+theme_bw()+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.line = element_line(),panel.border =  element_blank())
  print(p3 + ggtitle(distance_matrix[2])+scale_colour_manual(values = groups_color))
  dev.off()

  ####without names and ellipse
  pdf(paste(category1,"_",distance_matrix[2], "_NMDS_without_labels.pdf", sep=""), width=7.6, height=6.6)
  p2 = plot_ordination(gpt, GP.ord, type="samples", color=category1)
  p3 = p2  + geom_point(size=3) +theme(text = element_text(size = 15))+theme_bw()+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.line = element_line(),panel.border =  element_blank())
  print(p3 + ggtitle(distance_matrix[2])+scale_colour_manual(values = groups_color))
  dev.off()


  write.table(as.matrix(GP.ord$points), NMDS_ordtxtname, quote=FALSE, col.names=NA, sep="\t")
}

print("#Generate the PCoA 2D plot for betadiversity")
for (distance_matrix in list(c('bray','bray_curtis'), c('unifrac','unweighted_unifrac'), c('wunifrac','weighted_unifrac'))){
  GP.ord <- ordinate(gpt, "PCoA", distance_matrix[1])
  PCoA_outputpdfname <- paste(category1,"_",distance_matrix[2], "_PCoA_2D.pdf", sep="")
  PCoA_ordtxtname <- paste(category1,"_",distance_matrix[2], "_PCoA.ord.txt", sep="")


  pdf(PCoA_outputpdfname, width=7.6, height=6.6)
  p2 = plot_ordination(gpt, GP.ord, type="samples", color=category1) 
  p3 = p2  + geom_point(size=3) + geom_text_repel(aes(label=Description),hjust=0, vjust=2, size=4) + stat_ellipse()+theme(text = element_text(size = 15))+theme_bw()+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.line = element_line(),panel.border =  element_blank())
  print(p3 + ggtitle(distance_matrix[2])+scale_colour_manual(values = groups_color))
  dev.off()

  ######without names and ellipse
  pdf(paste(category1,"_",distance_matrix[2], "_PCoA_2D_without_labels.pdf", sep=""), width=7.6, height=6.6)
  p2 = plot_ordination(gpt, GP.ord, type="samples", color=category1) 
  p3 = p2  + geom_point(size=3)+theme(text = element_text(size = 15))+theme_bw()+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.line = element_line(),panel.border =  element_blank())
  print(p3 + ggtitle(distance_matrix[2])+scale_colour_manual(values = groups_color))
  dev.off()
  print("#Generate the PCoA 3D plot for betadiversity")
  ######pcoa 3d plot
  asign_rainbow_cor<-function(x){
    unig<-unique(x)
    color_out<-x
    for (i in 1:length(unig)) {
      color_out[color_out==unig[i]]<-groups_color[i]
    }
    return(color_out)
  }

  gp<-as.character(data.frame(sample_data(gpt))[category1][,1])
  gp_ord<-order(gp)
  gp<-gp[gp_ord]
  tdata<-GP.ord$vectors[,1:3][gp_ord,]
  eig<-data.frame(GP.ord$values)["Eigenvalues"][,1]
  lab<-paste("Axis.",c(1:3)," [",round((eig[1:3]/sum(eig))*100,digits=2),"%]",sep="")
  pdf(paste(category1,"_",distance_matrix[2], "_PCoA_3D.pdf", sep=""), width=8, height=6.4)
  opar<-par(no.readonly=TRUE)
  par(fig=c(0,0.75,0,1))
  if(ifline){
    scatterplot3d(tdata,mar=c(2.2,2.2,0,0)+1,xlab=lab[1],ylab=lab[2],zlab=lab[3],color=asign_rainbow_cor(gp), grid=TRUE, box=F, type="h", lty.hplot=2, pch=19)
  }else{
    scatterplot3d(tdata,mar=c(2.2,2.2,0,0)+1,xlab=lab[1],ylab=lab[2],zlab=lab[3],color=asign_rainbow_cor(gp), grid=TRUE, box=F, pch=19)
  }
  par(fig=c(0.75,1,0,1),xpd=TRUE)
  legend("center", legend = unique(gp), bty = 'n',xpd = TRUE,horiz = FALSE,col = groups_color, pch = 19, inset = -0.1)
  par(opar)
  dev.off()
  write.table(as.matrix(tdata), PCoA_ordtxtname, quote=FALSE, col.names=NA, sep="\t")

}




####################Using mixOmics for PLS-DA plot
#plsdatxt = paste(this.dir, "/R_output/feature-table.PLSDA.txt", sep="")


library(knitr)
knitr::opts_chunk$set(dpi = 100, echo= TRUE, warning=FALSE, message=FALSE, fig.align = 'center', 
                      fig.show=TRUE, fig.keep = 'all', out.width = '50%') 

## ----message = FALSE-----------------------------------------------------
library(mixOmics)

## ------------------------------------------------------------------------
#srbct <- load("/Users/chengguo/Downloads/PLSDA_SRBCT/result-SRBCT-sPLSDA.RData")
X = read.table(txt,skip=1, head=TRUE,comment.char = "",row.names = 1,sep = "\t",check.names=F)
taxonomy<-X[,length(X)]
X<-X[,-length(X)]
#head(X)
tX<-t(X)


#rownames(tX)<-colnames(X)
#head(tX)
A = read.table(map, header = T,row.names = 1,comment.char = "",sep = "\t",check.names = F,na.strings = "")
Y = A[category1][,1]

tX<-tX[match(rownames(A),rownames(tX)),]
taxonomy<-taxonomy[colSums(tX)>10]
tX<-tX[,colSums(tX)>10]

#print(dim(tX)[2])
#print(length(taxonomy))


pca.srbct = pca(tX, ncomp = 3, center = TRUE, scale = TRUE)
#pca.srbct #outputs the explained variance per component
plot(pca.srbct)  # screeplot of the eingenvalues (explained variance per component)
pdf(paste(category1,"_","PCA_plot.pdf",sep=""), width = 7.6, height = 6.6)
plotIndiv(pca.srbct, group = Y,col.per.group = groups_color, ind.names = FALSE, legend = TRUE, ellipse = TRUE, title = 'PCA plot')
dev.off()
#write.table(as.matrix(pca.srbct), “PCA_ord.txt”, quote=FALSE, col.names=NA, sep="\t")


srbct.plsda <- plsda(tX, Y)  # set ncomp to 10 for performance assessment later
plsda.vip <- vip(srbct.plsda)
write.table(data.frame(OTUID=rownames(plsda.vip),plsda.vip,Taxonomy=taxonomy),paste(category1,"_","PLSDA_Variable_importance_in_projection.txt"),row.names = F,sep="\t")

pdf(paste(category1,"_","PLSDA_AUC_plot.pdf",sep=""), width = 6, height = 6)
auroc(srbct.plsda, roc.comp = 2)
dev.off()

pdf(paste(category1,"_","PLSDA_comp_plot.pdf",sep=""), width = 6, height = 4)
plotIndiv(srbct.plsda ,col.per.group = groups_color, comp = 1:2, group = Y, ellipse.level = 0.75,size.xlabel = 15, size.ylabel = 15,size.axis = 15,size.legend = 15,size.legend.title = 15,ind.names = FALSE, title = "Supervised PLS-DA on OTUs",abline = T,legend = TRUE,ellipse = T)
dev.off()
#write.table(as.matrix(srbct.plsda), “PLSDA_ord.txt”, quote=FALSE, col.names=NA, sep="\t")



##################Use ggplot to draw alpha diversity boxplot. 
library("ggplot2") # load related packages
#read files.
alphadatxt = paste(this.dir, "/alpha/alpha-summary.tsv", sep="")
#alphameta = paste(this.dir, "/alpha/sample-metadata_alphadiversity.txt", sep="")
design = read.table(map, header=T,row.names= 1,comment.char="", check.names=F,sep="\t",stringsAsFactors = F,na.strings = "") 
alpha = read.table(alphadatxt, header=T, row.names= 1, sep="\t")

# merge information for script
index = cbind(design, alpha[match(rownames(design), rownames(alpha)), ])
# run shannon, observed_otus, faith_pd separately as the aes function is not accepting variables!!! Hard coded for Group1 as well. Really bad script.

for(alpha_index in colnames(alpha)){
  p = ggplot(index, aes_string(x=category1, y=alpha_index, color=category1)) + geom_boxplot(alpha=1, outlier.size=0, size=0.7, width=0.5, fill="transparent") +  geom_jitter( position=position_jitter(0.17), size=1, alpha=0.7) + labs(x="Groups", y=paste(alpha_index," index",sep = ""))+theme_bw()+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.line = element_line(),panel.border =  element_blank())+theme(axis.text.x = element_text(angle = 45,size = 10,hjust = 1))+scale_colour_manual(values = groups_color)
  ggsave(paste(category1,"_alpha_diversity_",alpha_index,".boxplot.pdf", sep=""), p, width = 6, height = 3)
}


#p = ggplot(index, aes_string(x=category1, y="observed_otus", color=category1)) + geom_boxplot(alpha=1, outlier.size=0, size=0.7, width=0.5, fill="transparent") +  geom_jitter( position=position_jitter(0.17), size=1, alpha=0.7) + labs(x="Groups", y="observed_otus index")+theme_bw()+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.line = element_line(),panel.border =  element_blank())+theme(axis.text.x = element_text(angle = 45,size = 10,hjust = 1))
#ggsave(paste(category1,"_","alpha_diversity_observed_otus.boxplot.pdf", sep=""), p, width = 6, height = 3)

#p = ggplot(index, aes_string(x=category1, y="shannon", color=category1)) + geom_boxplot(alpha=1, outlier.size=0, size=0.7, width=0.5, fill="transparent") +  geom_jitter( position=position_jitter(0.17), size=1, alpha=0.7) + labs(x="Groups", y="shannon index")+theme_bw()+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.line = element_line(),panel.border =  element_blank())+theme(axis.text.x = element_text(angle = 45,size = 10,hjust = 1))
#ggsave(paste(category1,"_","alpha_diversity_shannon.boxplot.pdf", sep=""), p, width = 6, height = 3)

#p = ggplot(index, aes_string(x=category1, y="faith_pd", color=category1)) + geom_boxplot(alpha=1, outlier.size=0, size=0.7, width=0.5, fill="transparent") +  geom_jitter( position=position_jitter(0.17), size=1, alpha=0.7) + labs(x="Groups", y="faith_pd index")+theme_bw()+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.line = element_line(),panel.border =  element_blank())+theme(axis.text.x = element_text(angle = 45,size = 10,hjust = 1))
#ggsave(paste(category1,"_","alpha_diversity_faith_pd.boxplot.pdf", sep=""), p, width = 6, height = 3)
