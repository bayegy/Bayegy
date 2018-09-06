##########Library import
library("ape")
library("phyloseq")
library("ggplot2")
library("gplots")
library("vegan")
library("ggrepel")
library("pheatmap")

args <- commandArgs(trailingOnly = TRUE)

##########Data import
#txt = "/Users/chengguo/Desktop/Hengchuang/M231/exported/feature-table.ConsensusLineage.txt"
#tre = "/Users/chengguo/Desktop/Hengchuang/M231/exported/tree.rooted.nwk"
#rs = "/Users/chengguo/Desktop/Hengchuang/M231/exported/dna-sequences.fasta"
map = args[1]
category1=args[2]

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
tre = paste(this.dir, "/exported/tree.rooted.nwk", sep="")
rs = paste(this.dir, "/exported/dna-sequences.fasta", sep="")

print("#Start reading realted files with Phyloseq")
print(map)
print(txt)
print(tre)
print(rs)


########################Using Phyloseq
qiimedata = import_qiime(txt, map, tre, rs)

gpt <- subset_taxa(qiimedata, Kingdom=="Bacteria")
#gpt <- prune_taxa(names(sort(taxa_sums(gpt),TRUE)[1:30]), gpt)
gpt <- prune_taxa(names(sort(taxa_sums(gpt),TRUE)), gpt)

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
for (distance_matrix in c('bray', 'unifrac', 'jaccard', 'wunifrac')){
  GP.ord <- ordinate(gpt, "NMDS", distance_matrix)
  NMDS_outputpdfname <- paste(category1,"_",distance_matrix, "_NMDS.pdf", sep="")
  NMDS_ordtxtname <- paste(distance_matrix, "_NMDS.ord.txt", sep="")
  pdf(NMDS_outputpdfname, width=6.6, height=6.6)
  p2 = plot_ordination(gpt, GP.ord, type="samples", color=category1) 
  p3 = p2  + geom_point(size=3) + geom_text_repel(aes(label=Description),hjust=0, vjust=2, size=4) + stat_ellipse()+theme(text = element_text(size = 15))
  print(p3 + ggtitle(distance_matrix))
  dev.off()
  write.table(as.matrix(GP.ord), NMDS_ordtxtname, quote=FALSE, col.names=NA, sep="\t")
}

print("#Generate the PCoA 2D plot for betadiversity")
for (distance_matrix in c('bray', 'unifrac', 'jaccard', 'wunifrac')){
  GP.ord <- ordinate(gpt, "PCoA", distance_matrix)
  PCoA_outputpdfname <- paste(category1,"_",distance_matrix, "_PCoA_2D.pdf", sep="")
  PCoA_ordtxtname <- paste(distance_matrix, "_PCoA_2D.ord.txt", sep="")
  pdf(PCoA_outputpdfname, width=6.6, height=6.6)
  p2 = plot_ordination(gpt, GP.ord, type="samples", color=category1) 
  p3 = p2  + geom_point(size=3) + geom_text_repel(aes(label=Description),hjust=0, vjust=2, size=4) + stat_ellipse()+theme(text = element_text(size = 15))
  print(p3 + ggtitle(distance_matrix))
  dev.off()
  write.table(as.matrix(GP.ord), PCoA_ordtxtname, quote=FALSE, col.names=NA, sep="\t")
}


print("#calculate distance")
for (distance_matrix in c('bray', 'unifrac', 'jaccard', 'wunifrac')){
  beta_heatmap_outputpdfname <- paste(distance_matrix, "_betadiversity_summary.pdf", sep="")
  #pdf(beta_heatmap_outputpdfname)
  Dist <- distance(qiimedata, method=distance_matrix)

  beta_outputtxtname <- paste(distance_matrix, "_matrix.txt", sep="")
  write.table(as.matrix(Dist), beta_outputtxtname , quote=FALSE, col.names=NA, sep="\t")
  
  Dist_read<-read.table(beta_outputtxtname, head=T)
  pdf(beta_heatmap_outputpdfname)
  pheatmap(Dist_read,fontsize=10,border_color = "black",fontsize_row =10,
           cluster_rows=T,clustering_distance_rows="euclidean",
           cluster_cols=T,clustering_distance_cols="euclidean",
           clustering_method="centroid")
  dev.off()
  
}


####################Using mixOmics for PLS-DA plot
plsdatxt = paste(this.dir, "/R_output/feature-table.PLSDA.txt", sep="")


library(knitr)
knitr::opts_chunk$set(dpi = 100, echo= TRUE, warning=FALSE, message=FALSE, fig.align = 'center', 
                      fig.show=TRUE, fig.keep = 'all', out.width = '50%') 

## ----message = FALSE-----------------------------------------------------
library(mixOmics)

## ------------------------------------------------------------------------
#srbct <- load("/Users/chengguo/Downloads/PLSDA_SRBCT/result-SRBCT-sPLSDA.RData")
X = read.table(plsdatxt, head=TRUE,comment.char = "",row.names = 1,sep = "\t")
#head(X)
tX<-t(X)
#head(tX)
A = read.table(map, header = T,row.names = 1,comment.char = "",sep = "\t")
A=A[match(rownames(tX),rownames(A)),]
Y = A[category1][,1]

pca.srbct = pca(tX, ncomp = 3, center = TRUE, scale = TRUE)
#pca.srbct #outputs the explained variance per component
plot(pca.srbct)  # screeplot of the eingenvalues (explained variance per component)
pdf(paste(category1,"_","PCA_plot.pdf",sep=""), width = 6.6, height = 6.6)
plotIndiv(pca.srbct, group = Y, ind.names = FALSE, legend = TRUE, ellipse = TRUE, title = 'PCA plot')
dev.off()
#write.table(as.matrix(pca.srbct), “PCA_ord.txt”, quote=FALSE, col.names=NA, sep="\t")


srbct.plsda <- plsda(tX, Y)  # set ncomp to 10 for performance assessment later
plsda.vip <- vip(srbct.plsda)
write.csv(plsda.vip,paste(category1,"_","PLSDA_Variable_importance_in_projection.txt"))

pdf(paste(category1,"_","PLSDA_AUC_plot.pdf",sep=""), width = 5, height = 4)
auroc(srbct.plsda, roc.comp = 2)
dev.off()

pdf(paste(category1,"_","PLSDA_comp_plot.pdf",sep=""), width = 6, height = 6)
plotIndiv(srbct.plsda , comp = 1:2, group = Y, ellipse.level = 0.75,size.xlabel = 15, size.ylabel = 15,size.axis = 15,size.legend = 15,size.legend.title = 15,ind.names = FALSE, title = "Supervised PLS-DA on OTUs",abline = T,legend = TRUE,ellipse = T)
dev.off()
#write.table(as.matrix(srbct.plsda), “PLSDA_ord.txt”, quote=FALSE, col.names=NA, sep="\t")



##################Use ggplot to draw alpha diversity boxplot. 
library("ggplot2") # load related packages
#read files.
alphadatxt = paste(this.dir, "/alpha/alpha-summary.tsv", sep="")
alphameta = paste(this.dir, "/alpha/sample-metadata_alphadiversity.txt", sep="")
design = read.table(alphameta, header=T, row.names= 1, sep="\t") 
alpha = read.table(alphadatxt, header=T, row.names= 1, sep="\t")

# merge information for script
index = cbind(alpha, design[match(rownames(alpha), rownames(design)), ]) 
# run shannon, observed_otus, faith_pd separately as the aes function is not accepting variables!!! Hard coded for Group1 as well. Really bad script.
p = ggplot(index, aes_string(x=category1, y="observed_otus", color=category1)) + geom_boxplot(alpha=1, outlier.size=0, size=0.7, width=0.5, fill="transparent") +  geom_jitter( position=position_jitter(0.17), size=1, alpha=0.7) + labs(x="Groups", y="observed_otus index")
ggsave(paste(category1,"_","alpha_diversity_observed_otus.boxplot.pdf", sep=""), p, width = 5, height = 3)

p = ggplot(index, aes_string(x=category1, y="shannon", color=category1)) + geom_boxplot(alpha=1, outlier.size=0, size=0.7, width=0.5, fill="transparent") +  geom_jitter( position=position_jitter(0.17), size=1, alpha=0.7) + labs(x="Groups", y="shannon index")
ggsave(paste(category1,"_","alpha_diversity_shannon.boxplot.pdf", sep=""), p, width = 5, height = 3)

p = ggplot(index, aes_string(x=category1, y="faith_pd", color=category1)) + geom_boxplot(alpha=1, outlier.size=0, size=0.7, width=0.5, fill="transparent") +  geom_jitter( position=position_jitter(0.17), size=1, alpha=0.7) + labs(x="Groups", y="faith_pd index")
ggsave(paste(category1,"_","alpha_diversity_faith_pd.boxplot.pdf", sep=""), p, width = 5, height = 3)
