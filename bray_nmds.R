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
  make_option(c("-o", "--output"),metavar="directory",dest="out", help="Specify the directory of output files. default=./",default="./")
)

opt <- parse_args(OptionParser(option_list=option_list,description = "This script is used to compare the predicted pathway of level 3"))

library(getopt)
library(ggplot2)
library(ggrepel)
library(stringr)
library(vegan)
options(scipen = 200)

if(!dir.exists(opt$out)){dir.create(opt$out,recursive = T)}

if(is.null(opt$colors)){
  base_dir<-normalizePath(dirname(get_Rscript_filename()))
  source(paste(base_dir,"/piputils/get_colors.R", sep = ""))
  groups_color<-get_colors(opt$group, opt$map)
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

# NMDS does not allow negative values
# input <- scale(input)
# df.pca <- prcomp(input, center = TRUE, scale. = TRUE)
nmds <- metaMDS(input, distance = 'bray', k = 2)

label <- choose(as.logical(opt$label), rownames(nmds$points), NA)
# len_g <- length(unique(groups))

# ss <- apply(df.pca$x, 2, var)
# pct <- ss/sum(ss)

# xylab <- paste(c("PC1 (", "PC2 ("), round(pct[1:2]*100, 2), "%)", sep = "")
xylab <- c("MDS1", "MDS2")

p<-ggplot(data = data.frame(nmds$points), aes(x=MDS1, y=MDS2, color= groups))+
  theme_classic()+
  labs(x=xylab[1], y=xylab[2], title = paste('Stress =', round(nmds$stress, 4)))+
  # geom_abline(data = data.frame(slope=c(0,91), intercept=c(0,0)), mapping = aes(slope=slope, intercept=intercept), size=0.1)+
  geom_hline(yintercept = 0, size=0.1, linetype="dashed") + geom_vline(xintercept = 0, size=0.1, linetype="dashed") +
  # scale_x_continuous(limits = c(-3, 3)) + scale_y_continuous(limits = c(-3, 3))+
  geom_point(size=4)+
  geom_text_repel(label=label, size=3) +
  scale_color_manual(values = groups_color) +
  choose(as.logical(opt$ellipse), stat_ellipse(level = 0.95, type = "norm", size=0.1, segments = 300), theme())+
  theme(legend.title = element_text(size = 0),
        text = element_text(size = as.numeric(opt$fontsize)),
        plot.title = element_text(hjust = 0.5))


ggsave(plot=p, file=paste(opt$out, "/NMDS.pdf", sep = ""), width=6, height=5.6)
# save the bray-curtis distance matrix
dis<-vegdist(input, method = "bray")
write.table(as.matrix(dis), file=paste(opt$out, "/bray_curtis_distance_matrix.xls", sep = ""), quote=FALSE, col.names=NA, sep="\t")
write.table(as.matrix(nmds$points), file=paste(opt$out, "/nmds_points_ordinates.xls", sep = ""), quote=FALSE, col.names=NA, sep="\t")
# save the bray-curtis distance matrix
# groups <- groups[match(rownames(input), rownames(groups)), ]