library(optparse)
#######arguments
option_list <- list(
  make_option(c("-i", "--input"),metavar="path", dest="input", help="Abundance table. Required",default=NULL),
  make_option(c("-a", "--alpha"),metavar="path", dest="alpha", help="alpha index table.",default=NULL),
  make_option(c("-t", "--tree"),metavar="tree.nwk", dest="tree", help="tree with nwk format",default=NULL),
  make_option(c("-m", "--map"),metavar="path",dest="map", help="Specify the path of mapping file. Required",default=NULL),
  make_option(c("-c", "--category"),metavar="string",dest="category", help="Category to compare. Required",default=NULL),
  make_option(c("-C", "--colors"),metavar="string",dest="colors", help="Comma seprated group colors.",default=NULL),
  make_option(c("-e", "--ellipse"),metavar="logical",dest="ellipse", help="draw ellipse or not", default=NULL),
  make_option(c("-l", "--line"),metavar = "logical",dest="line", help="If TRUE, plot dotted line in 3D pcoa plot",default = FALSE),
  make_option(c("-k", "--skip"),metavar = "logical",dest="skip", help="IF TRUE, skip the first line.",default=FALSE),
  make_option(c("--output-pcoa"),metavar="directory",dest="pcoa", help="Specify the directory of output pcoa files. default is not output",default=NULL),
  make_option(c("--output-nmds"),metavar="directory",dest="nmds", help="Specify the directory of output nmds files. default is not output",default=NULL),
  make_option(c("--output-matrix"),metavar="directory",dest="matrix", help="Specify the directory of output distance matrix files. default is not output",default=NULL),
  make_option(c("--output-plsda"),metavar="directory",dest="plsda", help="Specify the directory of output plsda files. default is not output",default=NULL),
  make_option(c("--output-pca"),metavar="directory",dest="pca", help="Specify the directory of output pca files. default is not output",default=NULL),
  make_option(c("--output-alpha-heatmap"),metavar="directory",dest="heatmap", help="Specify the directory of output heatmap files. default is not output",default=NULL)
)

opt <- parse_args(OptionParser(option_list=option_list, description = "Base visulization of alpha and beta index"))


optl <- list()

for(idx in c('pcoa', 'nmds', 'matrix', 'plsda', 'pca', 'heatmap')){
  drt <- opt[[idx]]
  not_null <- !is.null(drt)
  optl[[idx]] <- not_null
  if(not_null){
    if(!dir.exists(drt)){dir.create(drt,recursive = T)}
    if(!endsWith(drt, '/')){drt<-opt[[idx]]<-paste(drt, '/', sep = "")}
    TMP_DIR <- drt
  }
}

choose<-function(condition,choice1,choice2){
  if(condition){
    return(choice1)
  }else{
    return(choice2)
  }
}


##########Library import
library("ape")
library("phyloseq")
library("ggplot2")
# library("gplots")
library("vegan")
library("ggrepel")
library("pheatmap")
library("scatterplot3d")
library("stringr")
library(getopt)

library('cowplot')

input = opt$input
map = opt$map
category = opt$category
skip = as.logical(opt$skip)

tmp_map_exists <- FALSE
df_map <- read.table(map, sep="\t", na.strings="", header = TRUE, comment.char = "", check.names = F, stringsAsFactors = F)
if(!"Description"%in%colnames(df_map)){
  df_map$Description <- df_map[, 1]
  map <- paste(TMP_DIR, "mappgin_file.txt", sep = "")
  write.table(df_map, map, quote = FALSE, sep = "\t", na="", row.names = FALSE)
  tmp_map_exists <- TRUE
}



if(is.null(opt$colors)){
  base_dir<-normalizePath(dirname(get_Rscript_filename()))
  source(paste(base_dir,"/piputils/get_colors.R", sep = ""))
  groups_color<-get_colors(category, map)
}else{
  groups_color<-str_split(opt$colors, ",")[[1]]
}



if(optl$nmds||optl$matrix||optl$pcoa){
  ########################Using Phyloseq
  gpt = import_qiime(otufilename = input, mapfilename=map, treefilename=opt$tree)
  otu<-gpt@otu_table@.Data
  sum_of_otus<-colSums(t(otu))
  selected_otu<-names(sum_of_otus)[sum_of_otus>0]
  gpt <- prune_taxa(selected_otu, gpt)
}


if(optl$nmds||optl$matrix){
  print("#Generate the NMDS plot for betadiversity")
  for (distance_matrix in list(c('bray','bray_curtis'), c('unifrac','unweighted_unifrac'), c('wunifrac','weighted_unifrac'))){
    GP.distance <- distance(physeq = gpt,
                            method = distance_matrix[1],
                            type = "samples")
    if(optl$matrix){
      Dist <- as.matrix(GP.distance)
      write.table(Dist, paste(opt$matrix, distance_matrix[2], "_matrix.xls", sep="") , quote=FALSE, col.names=NA, sep="\t")
    }
    if(optl$nmds){
      GP.ord <- ordinate(gpt, "NMDS", distance_matrix[1])
      nmds <- metaMDS( comm = as.dist(GP.distance) )
      GP.stress <- round(nmds$stress,4)
      title <- paste(distance_matrix[2], '   NMDS stress:', as.character(GP.stress), sep = ' ')
      pdf(paste(opt$nmds, category,"_",distance_matrix[2], "_NMDS.pdf", sep=""), width=7.6, height=6.6)
      p2 = plot_ordination(gpt, GP.ord, type="samples", color=category)
      p3 = p2  + geom_point(size=3) + geom_text_repel(aes(label=Description),hjust=0, vjust=2, size=4) + choose(is.null(opt$ellipse)||as.logical(opt$ellipse), stat_ellipse(), theme())+theme(text = element_text(size = 15))+theme_bw()+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.line = element_line(),panel.border =  element_blank())
      print(p3 + ggtitle(title) + scale_colour_manual(values = groups_color))
      dev.off()

      ####without names and ellipse
      pdf(paste(opt$nmds, category,"_",distance_matrix[2], "_NMDS_without_labels.pdf", sep=""), width=7.6, height=6.6)
      p2 = plot_ordination(gpt, GP.ord, type="samples", color=category)
      p3 = p2  + geom_point(size=3) +theme(text = element_text(size = 15))+theme_bw()+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.line = element_line(),panel.border =  element_blank())
      if(!is.null(opt$ellipse)&&as.logical(opt$ellipse)){
        p3 = p3 + stat_ellipse()
      }
      print(p3 + ggtitle(title)+scale_colour_manual(values = groups_color))
      dev.off()
      write.table(as.matrix(GP.ord$points), paste(opt$nmds,category,"_",distance_matrix[2], "_NMDS.ord.xls", sep=""), quote=FALSE, col.names=NA, sep="\t")
    }
  }
}

if(optl$pcoa){
  print("#Generate the PCoA 2D plot for betadiversity")
  for (distance_matrix in list(c('bray','bray_curtis'), c('unifrac','unweighted_unifrac'), c('wunifrac','weighted_unifrac'))){
    GP.ord <- ordinate(gpt, "PCoA", distance_matrix[1])
    pdf(paste(opt$pcoa, category, "_", distance_matrix[2], "_PCoA_2D.pdf", sep=""), width=7.6, height=6.6)
    p2 = plot_ordination(gpt, GP.ord, type="samples", color=category)
    p3 = p2  + geom_point(size=3) + geom_text_repel(aes(label=Description),hjust=0, vjust=2, size=4) + choose(is.null(opt$ellipse)||as.logical(opt$ellipse), stat_ellipse(), theme())+theme(text = element_text(size = 15))+theme_bw()+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.line = element_line(),panel.border =  element_blank())
    print(p3 + ggtitle(distance_matrix[2])+scale_colour_manual(values = groups_color))
    dev.off()

    ######without names and ellipse
    pdf(paste(opt$pcoa, category,"_",distance_matrix[2], "_PCoA_2D_without_labels.pdf", sep=""), width=7.6, height=6.6)
    p2 = plot_ordination(gpt, GP.ord, type="samples", color=category)
    p3 = p2  + geom_point(size=3)+theme(text = element_text(size = 15))+theme_bw()+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.line = element_line(),panel.border =  element_blank())
    if(!is.null(opt$ellipse)&&as.logical(opt$ellipse)){
      p3 = p3 + stat_ellipse()
    }
    print(p3 + ggtitle(distance_matrix[2])+scale_colour_manual(values = groups_color))
    dev.off()
    print("#Generate the PCoA 3D plot for betadiversity")
    ######pcoa 3d plot
    asign_color<-function(x){
      unig<-unique(x)
      color_out<-x
      for (i in 1:length(unig)) {
        color_out[color_out==unig[i]]<-groups_color[i]
      }
      return(color_out)
    }
    gp<-as.character(data.frame(sample_data(gpt))[category][,1])
    gp_ord<-order(gp)
    gp<-gp[gp_ord]
    tdata<-GP.ord$vectors[,1:3][gp_ord,]
    eig<-data.frame(GP.ord$values)["Eigenvalues"][,1]
    lab<-paste("Axis.",c(1:3)," [",round((eig[1:3]/sum(eig))*100,digits=2),"%]",sep="")

    pdf(paste(opt$pcoa, category,"_",distance_matrix[2], "_PCoA_3D.pdf", sep=""), width=8, height=6.4)
    opar<-par(no.readonly=TRUE)
    par(fig=c(0,0.75,0,1))
    if(as.logical(opt$line)){
      scatterplot3d(tdata,mar=c(2.2,2.2,0,0)+1,xlab=lab[1],ylab=lab[2],zlab=lab[3],color=asign_color(gp), grid=TRUE, box=F, type="h", lty.hplot=2, pch=19)
    }else{
      scatterplot3d(tdata,mar=c(2.2,2.2,0,0)+1,xlab=lab[1],ylab=lab[2],zlab=lab[3],color=asign_color(gp), grid=TRUE, box=F, pch=19)
    }
    par(fig=c(0.75,1,0,1),xpd=TRUE)
    legend("center", legend = unique(gp), bty = 'n',xpd = TRUE,horiz = FALSE,col = groups_color, pch = 19, inset = -0.1)
    par(opar)
    dev.off()
    write.table(as.matrix(tdata), paste(opt$pcoa, category,"_",distance_matrix[2], "_PCoA.ord.xls", sep=""), quote=FALSE, col.names=NA, sep="\t")
  }
}


library(knitr)
knitr::opts_chunk$set(dpi = 100, echo= TRUE, warning=FALSE, message=FALSE, fig.align = 'center',
                      fig.show=TRUE, fig.keep = 'all', out.width = '50%')

if(optl$pca||optl$plsda){
  ## ----message = FALSE-----------------------------------------------------
  library(mixOmics)
  ## ------------------------------------------------------------------------
  #srbct <- load("/Users/chengguo/Downloads/PLSDA_SRBCT/result-SRBCT-sPLSDA.RData")
  X = read.table(input, skip=ifelse(skip, 1, 0), head=TRUE, comment.char = "", row.names = 1,sep = "\t",check.names=F)
  taxonomy<-X[,length(X)]
  X<-X[,-length(X)]
  tX<-t(X)
  A = read.table(map, header = T,row.names = 1,comment.char = "",sep = "\t",check.names = F,na.strings = "")
  Y = A[category][,1]
  tX<-tX[match(rownames(A),rownames(tX)),]
  taxonomy<-taxonomy[colSums(tX)>10]
  tX<-tX[,colSums(tX)>10]
  tX <- scale(tX)
}


if(optl$pca){
  df.pca <- prcomp(tX, center = TRUE, scale. = TRUE)
  label <- rownames(df.pca$x)
  len_g <- length(unique(Y))
  ss <- apply(df.pca$x, 2, var)
  pct <- ss/sum(ss)
  xylab <- paste(c("PC1 (", "PC2 ("), round(pct[1:2]*100, 2), "%)", sep = "")
  p<-ggplot(data = data.frame(scale(df.pca$x)), aes(x=PC1, y=PC2))+
    theme_classic()+
    labs(x=xylab[1], y=xylab[2])+
    # geom_abline(data = data.frame(slope=c(0,91), intercept=c(0,0)), mapping = aes(slope=slope, intercept=intercept), size=0.1)+
    geom_hline(yintercept = 0, size=0.1) + geom_vline(xintercept = 0, size=0.1) +
    # scale_x_continuous(limits = c(-3, 3)) + scale_y_continuous(limits = c(-3, 3))+
    geom_point(size=2, mapping = aes(color= Y))+
    geom_text_repel(label=label, size=2) +
    scale_color_manual(values = groups_color) +
    # stat_ellipse(level = 0.95, type = "norm", size=0.1, segments = 300)+
    theme(legend.title = element_text(size = 0), text = element_text(size = 9))
  # browser()
  ggsave(plot=p, file=paste(opt$pca, category,"_","PCA_plot.pdf",sep=""), width=7, height=6.6)


  p<-ggplot(data = data.frame(scale(df.pca$x)), aes(x=PC1, y=PC2))+
    theme_classic()+
    labs(x=xylab[1], y=xylab[2])+
    # geom_abline(data = data.frame(slope=c(0,91), intercept=c(0,0)), mapping = aes(slope=slope, intercept=intercept), size=0.1)+
    geom_hline(yintercept = 0, size=0.1) + geom_vline(xintercept = 0, size=0.1) +
    # scale_x_continuous(limits = c(-3, 3)) + scale_y_continuous(limits = c(-3, 3))+
    geom_point(size=2, mapping = aes(color= Y))+
    # geom_text_repel(label=label, size=2) +
    scale_color_manual(values = groups_color) +
    # stat_ellipse(level = 0.95, type = "norm", size=0.1, segments = 300)+
    theme(legend.title = element_text(size = 0), text = element_text(size = 9))
  # browser()
  ggsave(plot=p, file=paste(opt$pca, category,"_","PCA_plot_without_lables.pdf",sep=""), width=7.6, height=6.6)
}


if(optl$plsda){
  srbct.plsda <- plsda(tX, Y)  # set ncomp to 10 for performance assessment later
  plsda.vip <- vip(srbct.plsda)
  write.table(data.frame(OTUID=rownames(plsda.vip),plsda.vip,Taxonomy=taxonomy),paste(opt$plsda, category, "_", "PLSDA_Variable_importance_in_projection.xls"),row.names = F,sep="\t")

  pdf(paste(opt$plsda, category,"_","PLSDA_AUC_plot.pdf",sep=""), width = 9, height = 6)
  auroc(srbct.plsda, roc.comp = 2)
  dev.off()

  pdf(paste(opt$plsda, category,"_","PLSDA_comp_plot.pdf",sep=""), width = 9, height = 8)
  plotIndiv(srbct.plsda ,col.per.group = groups_color, comp = 1:2, group = Y, ellipse.level = 0.75,size.xlabel = 15, size.ylabel = 15,size.axis = 15,size.legend = 15,size.legend.title = 15,ind.names = FALSE, title = "Supervised PLS-DA on OTUs",abline = T,legend = TRUE,ellipse = T)
  dev.off()
}


if(optl$heatmap){
  design = read.table(map, header=T,row.names= 1,comment.char="", check.names=F,sep="\t",stringsAsFactors = F,na.strings = "") 
  alpha = read.table(opt$alpha, header=T, row.names= 1, sep="\t")
  # merge information for script
  index = cbind(design, alpha[match(rownames(design), rownames(alpha)), ])
  # run shannon, observed_otus, faith_pd separately as the aes function is not accepting variables!!! Hard coded for Group1 as well. Really bad script.

  for(alpha_index in colnames(alpha)){
    p = ggplot(index, aes_string(x=category, y=alpha_index, color=category)) + geom_boxplot(alpha=1, outlier.size=0, size=0.7, width=0.5, fill="transparent") +  geom_jitter( position=position_jitter(0.17), size=1, alpha=0.7) + labs(x="Groups", y=paste(alpha_index," index",sep = ""))+theme_bw()+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.line = element_line(),panel.border =  element_blank())+theme(axis.text.x = element_text(angle = 45,size = 10,hjust = 1))+scale_colour_manual(values = groups_color)
   # ggsave(paste(opt$heatmap, category,"_alpha_diversity_",alpha_index,".boxplot.pdf", sep=""), p, width = 6, height = 3)
    filename = paste(opt$heatmap, category,"_alpha_diversity_",alpha_index,".boxplot.pdf", sep="")
    print('save_plot(filename = filename, plot = p , base_height = 6, base_width = NULL)')
    print('# save_plot() function is from cowplot package.')
    save_plot(filename = filename, plot = p , base_height = 6, base_width = NULL)
  }
}


if(tmp_map_exists){
  system(sprintf("rm %s", map))
}
