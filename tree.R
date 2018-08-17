
library("ggtree")
library("stringr")
otu_table<-read.table("./exported/feature-table.taxonomy.txt",header = T,skip=1,row.names = 1,check.names = F,stringsAsFactors = F,sep = "\t",comment.char = "")
tax<-otu_table[,dim(otu_table)[2]]


groupInfo<-str_extract(tax,"p__[^;]{1,100}")
groupInfo[is.na(groupInfo)]<-"Unclassfied_phylum"
groupInfo <- split(rownames(otu_table), groupInfo)

groupInfo1<-str_extract(tax,"g__[^;]{1,100}")
groupInfo1[is.na(groupInfo1)]<-"Unclassfied_genus"
groupInfo1 <- split(rownames(otu_table), groupInfo1)


otu_table<-otu_table[,-dim(otu_table)[2]]
metadata<-read.table("/home/admin1/16S_pipline_testfile/database/sample-metadata.tsv",header = T,row.names=1,check.names = F,stringsAsFactors = F,sep = "\t",comment.char = "")
metadata<-metadata["Group2"]
metagroup<-metadata[,1][match(colnames(otu_table),rownames(metadata))]
otu_table=scale(t(otu_table))

mysum=function(x){
	return(tapply(x,INDEX = metagroup,FUN = sum))
}


data=t(apply(otu_table,2,mysum))
tree <- read.tree("AdditionalPhylogeny//tree.nwk")
tree<- groupOTU(tree, groupInfo,group_name = "Phylum")
tree<- groupOTU(tree, groupInfo1,group_name = "taxa")
pa1<-length(unique(metagroup))
pa<-c(0,-0.26,-0.026,-0.01,0,0.01,0.01,0.01,0.01)
p = ggtree(tree,aes(color=Phylum))+geom_tiplab(size=3, align=TRUE, linesize=.5,aes(label=taxa))
pdf(file="AdditionalPhylogeny//Group2_hylogenetic_tree_heatmap.pdf", width=10, height=10)
gheatmap(p, data, offset = 0.1, width=1, hjust=-3,colnames_offset_y=-0.3,colnames_offset_x=pa[pa1])+theme_tree2()+theme(legend.position = "right",text=element_text(size=15,face="bold"))
dev.off()

