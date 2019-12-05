#!usr/bin/env Rscript
library(optparse)

#######arguments
option_list <- list( 
  make_option(c("-i", "--input"),metavar="path", dest="otu",help="Specify the path of otu table",default=NULL),
  make_option(c("-t", "--tree"),metavar="path", dest="tree",help="rooted tree.",default=NULL),
  make_option(c("-m", "--map"),metavar="path",dest="map", help="Specify the path of mapping file. Required",default=NULL),
  make_option(c("-c", "--category"),metavar="string",dest="group", help="Category to compare. Required",default=NULL),
  make_option(c("-n", "--number"),metavar="int", dest="num",help="The number of species needed to be plotted, default is 50",default=50),
  make_option(c("-o", "--output"),metavar="directory",dest="out", help="Specify the directory of output files",default="./")
)
opt <- parse_args(OptionParser(option_list=option_list,description = "This script is used to draw the relative stack barplot of species"))



library("phyloseq")
library("ggtree")
library("ggplot2")

tree <- opt$tree
otu <- opt$otu
map <- opt$map
category <- opt$group
num <- as.numeric(opt$num)
out <- opt$out

if(!dir.exists(out)){
  dir.create(out, recursive = T)
}

tre_object <- import_qiime(otufilename = otu, treefilename = tree, mapfilename = map)
df1 <- tre_object@otu_table@.Data
sum_of_otus<-colSums(t(df1))
selected_otu <- base::tapply(X=sum_of_otus, INDEX=tre_object@tax_table@.Data[, 'Genus'], 
                             FUN = function(v){return(names(v)[which.max(v)])})

selected_otu <- as.vector(selected_otu)
sum_of_otus <- sum_of_otus[selected_otu]
top <- head(order(sum_of_otus, decreasing = T), num)
selected_otu <- selected_otu[top]
tre_object <- prune_taxa(selected_otu, tre_object)

df2 <- t(tre_object@otu_table@.Data)
df2 <- scale(df2)
groups <- tre_object@sam_data[, category][[1]]
groups_sum <- apply(df2, 2, function(x){return(tapply(x, groups, mean))})
groups_sum <- t(groups_sum)

# p <- ggtree(tre_object, aes(color=Phylum))+
#   geom_tiplab(size=4,align=TRUE, linesize=.5,aes(label=Genus))
#   # scale_color_discrete(breaks = t,name="Phylum")
# p1<-gheatmap(p, groups_sum)

genus <- tre_object@tax_table@.Data[, 'Genus']
max_char_len <- max(nchar(genus))
par1 <- length(levels(groups))
uni_phylum <- unique(tre_object@tax_table@.Data[, 'Phylum'])


p = ggtree(tre_object, aes(color=Phylum))+
  geom_tiplab(size=4,align=TRUE, linesize=.5,aes(label=Genus))+
  scale_color_discrete(breaks = uni_phylum,name="Phylum")

wd=par1*0.2
ofs=max_char_len*0.02 + 0.02*wd
p1<-gheatmap(p, groups_sum, offset = ofs, width=wd, hjust=0.5,font.size=2,colnames_offset_y=-0.4,colnames_angle=75)+theme(legend.position = "right",text=element_text(size=15),axis.ticks=element_blank())
file_name=sprintf("%s/%s_phylogenetic_tree", out, category)

df3 <- data.frame(groups_sum, tre_object@tax_table@.Data)

write.table(df3, paste(file_name,"_table.xls",sep = ""), sep = "\t", quote=FALSE, col.names=NA)
ggsave(p1, file=paste(file_name,"_heatmap.pdf",sep=""), width=(wd+1)*4+3+max_char_len*0.03, height=10)


