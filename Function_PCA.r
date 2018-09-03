library("ggbiplot")
library("stringr")
# 菌群数据实战
# 读入实验设计
args <- commandArgs(trailingOnly = TRUE)
KEGG_function_txt = args[1]
meta_txt = args[2]
category_1=args[3]

this.dir <- getwd()
setwd(this.dir)

#setwd("~/Desktop")
#KEGG_function_txt = "feature-table.metagenome.L2.PCA.txt"
#meta_txt = "sample-metadata.tsv"
#category_1= "Group1"

design = read.table(meta_txt, header=T, row.names= 1, sep="\t",check.names = F) 

# 读取OTU表
otu_table = read.delim(KEGG_function_txt, row.names= 1,  header=T, sep="\t",check.names=F)
otu_table[otu_table == 0] = 0.0001
#otu_table

# 过滤数据并排序
idx = rownames(design) %in% colnames(otu_table) 
sub_design = design[idx,]
count = otu_table[, rownames(sub_design)]

# 基于OTU表PCA分析
#otu.pca <- prcomp(t(count), scale. = TRUE)

# 绘制PCA图，并按组添加椭圆
#p<-ggbiplot(otu.pca, obs.scale = 1, var.scale = 1, groups = sub_design$Group1, ellipse = TRUE,var.axes = F)


# 显著高丰度菌的影响

# 转换原始数据为百分比



norm = t(t(count)/colSums(count,na=T)) * 100 # normalization to total 100


#colnames(norm)<-colnames(count)
#rownames(norm)<-rownames(count)
# 筛选mad值大于0.5的OTU
mad.5 = norm[apply(norm,1,mad)>0.5,]
# 另一种方法：按mad值排序取前6波动最大的OTUs
mad.5 = head(norm[order(apply(norm,1,mad), decreasing=T),],n=6)
# 计算PCA和菌与菌轴的相关性
tmad<-t(mad.5)
#colnames(tmad)<-rownames(mad.5)
#rownames(tmad)<-colnames(mad.5)


otu.pca <- prcomp(tmad)

print(paste("Making PCA plots for", KEGG_function_txt, sep=" "))
KEGG_function_txt<-str_replace(KEGG_function_txt,"PCA.txt","")
###########The name is bad here use tools::file_path_sans_ext("ABCD.csv") to obtain the basename in the future.
PCA_plot_outputpdfname <- paste(KEGG_function_txt, category_1,".PCA.pdf", sep="")
pdf( PCA_plot_outputpdfname, width=12, height=12)
ggbiplot(otu.pca, obs.scale = 1, var.scale = 1, groups = sub_design[[category_1]], ellipse = TRUE,var.axes = T,varname.adjust=1,varname.size=3)
print(plot)
dev.off()


