#!/home/admin1/miniconda2/bin/Rscript
#utf-8
df_map<-function(df, func){
  out<-c()
  for(i in 1:ncol(df)){
    out[i]<-func(df[, i])
  }
  return(out)
}

choose<-function(condition,choice1,choice2){
  if(condition){
    return(choice1)
  }else{
    return(choice2)
  }
}

#'util
or<-function(choice1,choice2){
  if(length(choice1)==1){
      if(is.na(choice1)|choice1==FALSE){
        return(choice2)
      }
  }else{
      if(length(choice1)==0){
        return(choice2)
      }
  }
  return(choice1)
}


library(optparse)
#######arguments
option_list <- list(
  make_option(c("-i", "--input"),metavar="path", dest="input",help="Abundance table. Required",default=NULL),
  make_option(c("-m", "--map"),metavar="path",dest="map", help="Specify the path of mapping file. Required",default=NULL),
  make_option(c("-c", "--category"),metavar="string",dest="category", help="Category to compare. Required",default=NULL),
  make_option(c("-t", "--colors"),metavar="string",dest="colors", help="Comma seprated group colors.",default=NULL),
  make_option(c("-l", "--label"),metavar="logical",dest="label", help="Label sample name or not", default=FALSE),
  make_option(c("-f", "--fontsize"),metavar="int",dest="fontsize", help="font size", default=15),
  make_option(c("-e", "--ellipse"),metavar="logical",dest="ellipse", help="draw ellipse or not", default=TRUE),
  make_option(c("-v", "--line"),metavar="logical",dest="line", help="draw vertical line in 3D plot or not", default=TRUE),
  make_option(c("-p", "--prefix"),metavar="str", dest="prefix",help="The prefix of output files, default if ''",default=""),
  make_option(c("-o", "--output"),metavar="directory",dest="out", help="Specify the directory of output files. default=./",default="./")
)

opt <- parse_args(OptionParser(option_list=option_list,description = "This script is used to compare the predicted pathway of level 3"))

library(getopt)
library(ggplot2)
library(ggrepel)
library(stringr)
library("scatterplot3d")
if(!dir.exists(opt$out)){dir.create(opt$out,recursive = T)}

opt$out<-paste(opt$out,"/",opt$prefix, sep="")

if(is.null(opt$colors)){
  base_dir<-normalizePath(dirname(get_Rscript_filename()))
  source(paste(base_dir,"/piputils/get_colors.R", sep = ""))
  groups_color<-get_colors(opt$category, opt$map)
}else{
  groups_color<-str_split(opt$colors, ",")[[1]]
}


input<-read.table(opt$input, na.strings="", comment.char="", row.names=1, quote = "", check.names=F, stringsAsFactors=F, header = TRUE, sep = "\t")
map<-read.table(opt$map, na.strings="", comment.char="", row.names=1, quote = "", check.names=F, stringsAsFactors=F, header = TRUE, sep = "\t")


input <- input[, df_map(input, is.numeric)]
input <- t(input)
# print(opt$category)
groups <- map[opt$category]
groups <- na.omit(groups)

input <- input[match(rownames(groups), rownames(input)), ]
groups <- groups[, 1]

input <- input[, colSums(input!=0)>0]
input <- scale(input)
# print(input)
df.pca <- prcomp(input, center = TRUE, scale. = TRUE)


label <- choose(as.logical(opt$label), rownames(df.pca$x), NA)
# len_g <- length(unique(groups))

ss <- apply(df.pca$x, 2, var)
pct <- ss/sum(ss)

xylab <- paste(c("PC1 (", "PC2 (", "PC3 ("), round(pct[1:3]*100, 2), "%)", sep = "")

p<-ggplot(data = data.frame(df.pca$x), aes(x=PC1, y=PC2, color= groups))+
  theme_classic()+
  labs(x=xylab[1], y=xylab[2])+
  # geom_abline(data = data.frame(slope=c(0,91), intercept=c(0,0)), mapping = aes(slope=slope, intercept=intercept), size=0.1)+
  geom_hline(yintercept = 0, size=0.1, linetype="dashed") + geom_vline(xintercept = 0, size=0.1, linetype="dashed") +
  # scale_x_continuous(limits = c(-3, 3)) + scale_y_continuous(limits = c(-3, 3))+
  geom_point(size=4)+
  geom_text_repel(label=label, size=3) +
  scale_color_manual(values = groups_color) +
  choose(as.logical(opt$ellipse), stat_ellipse(level = 0.95, type = "norm", size=0.1, segments = 300), theme())+
  theme(legend.title = element_text(size = 0), text = element_text(size = as.numeric(opt$fontsize)))

ggsave(plot=p, file=paste(opt$out, "PCA.pdf", sep = ""), width=6, height=5.6)
write.table(as.matrix(df.pca$x), file=paste(opt$out, "pca_points_ordinates.xls", sep = ""), quote=FALSE, col.names=NA, sep="\t")
# groups <- groups[match(rownames(input), rownames(groups)), ]

# print("#Generate the PCoA 3D plot for betadiversity")
######pcoa 3d plot
asign_color<-function(x){
  unig<-sort(unique(x))
  color_out<-x
  for (i in 1:length(unig)) {
    color_out[color_out==unig[i]]<-groups_color[i]
  }
  return(color_out)
}

# gp<-groups
# gp_ord<-order(gp)
# gp<-gp[gp_ord]
tdata<-df.pca$x[, 1:3]
# eig<-data.frame(GP.ord$values)["Eigenvalues"][,1]
# xylab<-paste("Axis.",c(1:3)," [",round((eig[1:3]/sum(eig))*100,digits=2),"%]",sep="")
pdf(paste(opt$out, "PCA_3D.pdf", sep=""), width=8, height=6.4)
opar<-par(no.readonly=TRUE)
par(fig=c(0,0.75,0,1))
if(as.logical(opt$line)){
  scatterplot3d(tdata,mar=c(2.2,2.2,0,0)+1,xlab=xylab[1],ylab=xylab[2],zlab=xylab[3],color=asign_color(groups), grid=TRUE, box=F, type="h", lty.hplot=2, pch=19)
}else{
  scatterplot3d(tdata,mar=c(2.2,2.2,0,0)+1,xlab=xylab[1],ylab=xylab[2],zlab=xylab[3],color=asign_color(groups), grid=TRUE, box=F, pch=19)
}
par(fig=c(0.75,1,0,1),xpd=TRUE)
legend("center", legend = sort(unique(groups)), bty = 'n',xpd = TRUE,horiz = FALSE,col = groups_color, pch = 19, inset = -0.1)
par(opar)
dev.off()