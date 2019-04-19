# coding:utf-8
import argparse
import re
import sys
import os

# argument:
p = argparse.ArgumentParser(description="This script is used to plot phylogentic tree of top abundant otus")
p.add_argument('-i', '--input', dest='input', metavar='<path>', default=False,
               help='Specify feature_table.txt with taxonomy at last column and otu id at first column')
p.add_argument('-m', '--metadata', dest='metadata', metavar='<path>', default=False,
               help='Specify metadata with sample id at first column')
p.add_argument('-g', '--group', dest='group', metavar='<str>', default=False,
               help='Column name of group in metadata')
p.add_argument('-r', '--repseqs', dest='repseqs', metavar='<path>', default=False,
               help='Specify representive sequences file after masking and aligning')
p.add_argument('-o', '--outdir', dest='outdir', metavar='<directory>', default='./',
               help='Specify the output directory')
p.add_argument('-b', '--by-group-mean', dest='mean', action='store_true',
               help="Use group mean to plot heatmap")
p.add_argument('-n', '--number', dest='num', metavar='<int>', default=30,
               help='How many most abundant species do you want to analyze')

options = p.parse_args()


#os.system("if [ ! -d %s ];then mkdir -p %s;fi"%(options.outdir,options.outdir))
if not os.path.exists(options.outdir):
    os.makedirs(options.outdir)

########select otu####


with open('tree.R', 'w') as rscript:
    print('''
otu_table<-read.table("%s",header = T,skip=1,row.names = 1,check.names = F,stringsAsFactors = F,sep = "\\t",comment.char = "")
metadata<-read.table("%s",na.strings="",header = T,row.names=1,check.names = F,stringsAsFactors = F,sep = "\\t",comment.char = "")
metadata<-metadata["%s"]

####clean na
otu_table<-otu_table[,-dim(otu_table)[2]]
sel<-(!is.na(metadata))
sel1<-sel[match(colnames(otu_table),rownames(metadata))]
otu_table<-otu_table[,sel1]

otusum<-colSums(t(otu_table))
hold<-sort(otusum,T)[%s]
out<-data.frame(rownames(otu_table)[otusum>=hold])
write.table(out,"%s/selected_features.txt",sep = "",row.names = F,col.names = F,quote = F)
'''
          % (options.input, options.metadata, options.group, options.num, options.outdir),
          file=rscript)
os.system('Rscript tree.R')


#######align and mask#####
#os.system("qiime tools export  %s --output-dir %s/" % (options.repseqs, options.outdir))


#######select rep-seqs####
with open('%s/selected_features_reseqs.fasta' % (options.outdir), 'w') as fout:
    s_otuid = open('%s/selected_features.txt' % (options.outdir), 'r')
    s_otuid = s_otuid.read()

    s_otuid = re.split('\n', s_otuid)
    sn = []
    ln = 1
    for line in open(options.repseqs, 'r'):

        line = re.sub('\n$', '', line)
        line = re.sub('^>', '', line)
        if ln in sn:
            fout.write(line + "\n")

        if line in s_otuid:
            fout.write('>' + line + "\n")
            sn.append(ln + 1)
        ln += 1


######form tree#####
os.system("qiime tools import --input-path %s/selected_features_reseqs.fasta --output-path %s/selected_features_reseqs.qza --type 'FeatureData[AlignedSequence]'&&\
qiime phylogeny fasttree --i-alignment %s/selected_features_reseqs.qza --o-tree %s/selected_unrooted-tree.qza&&\
qiime phylogeny midpoint-root --i-tree %s/selected_unrooted-tree.qza --o-rooted-tree %s/selected_rooted-tree.qza&&\
qiime tools export %s/selected_rooted-tree.qza --output-dir %s/" % (options.outdir, options.outdir, options.outdir, options.outdir, options.outdir, options.outdir, options.outdir, options.outdir))

ifmean = 'T' if options.mean else 'F'
######visualize tree####
with open('tree.R', 'w') as rscript:
    print('''
library("ggtree")
library("stringr")
otu_table<-read.table("%s",header = T,skip=1,row.names = 1,check.names = F,stringsAsFactors = F,sep = "\\t",comment.char = "")
metadata<-read.table("%s",na.strings="",header = T,row.names=1,check.names = F,stringsAsFactors = F,sep = "\\t",comment.char = "")
gname="%s"

metadata=metadata[order(metadata[gname]),]


metadata<-metadata[gname]
metagroup<-metadata[,1]
#print(otu_table)
otu_table[,1:(ncol(otu_table)-1)]<-otu_table[,1:(ncol(otu_table)-1)][,match(rownames(metadata),colnames(otu_table)[-ncol(otu_table)])]

####clean na
otu_table<-otu_table[,c(!is.na(metagroup),T)]
metagroup<-metagroup[!is.na(metagroup)]


otusum<-colSums(t(otu_table[,-dim(otu_table)[2]]))
hold<-sort(otusum,T)[%s]

otu_table<-otu_table[otusum>=hold,]

tax<-otu_table[,dim(otu_table)[2]]

groupInfo<-str_extract(tax,"p__[^;]{1,100}")
groupInfo[is.na(groupInfo)]<-"Unclassfied_phylum"
leg<-as.factor(groupInfo)
groupInfo <- split(rownames(otu_table), groupInfo)



groupInfo1<-str_extract(tax,"g__[^;]{1,100}")
groupInfo1[is.na(groupInfo1)]<-"Unclassfied_genus"
groupInfo1 <- split(rownames(otu_table), groupInfo1)


otu_table<-otu_table[,-dim(otu_table)[2]]
#
otu_table=scale(t(otu_table))

mysum=function(x){
  return(tapply(x,INDEX = metagroup,FUN = sum))
}



if(%s){
  data=t(apply(otu_table,2,mysum))
}else{
  data=t(otu_table)
}

par1<-ncol(data)

tree <- read.tree("%s/tree.nwk")
tree<- groupOTU(tree, groupInfo,group_name = "Phylum")
tree<- groupOTU(tree, groupInfo1,group_name = "taxa")



t<-levels(leg)

par(mar=c(7, 4, 4, 2) + 0.1)
p = ggtree(tree,aes(color=Phylum))+
  geom_tiplab(size=4,align=TRUE, linesize=.5,aes(label=taxa))+
  scale_color_discrete(breaks = t,name="Phylum")
wd=par1*0.1
ofs=0.19
p1<-gheatmap(p, data, offset = ofs, width=wd, hjust=0.5,font.size=2,colnames_offset_y=-0.4,colnames_angle=75)+theme(legend.position = "right",text=element_text(size=15),axis.ticks=element_blank())
file_name="%s/%s_phylogenetic_tree"

write.table(data,paste(file_name,"_table.txt",sep = ""),sep = "\t",quote=FALSE,col.names=NA)
ggsave(p1,file=paste(file_name,"_heatmap.pdf",sep=""), width=(wd+1)*4+3, height=10)
'''
          % (options.input, options.metadata, options.group, options.num, ifmean, options.outdir, options.outdir, str(options.group)),
          file=rscript)

os.system('Rscript tree.R')
