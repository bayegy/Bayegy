
from __future__ import print_function,division
from optparse import OptionParser
from collections import defaultdict
import re,sys,os

#argument:
usage = 'python %prog [-i] [-m] [-g] [-r]'
p = OptionParser(usage = usage)
p.add_option('-i','--input',dest='input',metavar='*.txt',default=False,
			help='specify feature_table.txt with taxonomy at last column and otu id at first column')
p.add_option('-m','--metadata',dest='metadata',metavar='*.txt',default=False,
			help='specify metadata with sample id at first column')
p.add_option('-g','--group',dest='group',metavar='[group name]',default=False,
			help='column name of group in metadata')
p.add_option('-r','--repseqs',dest='repseqs',metavar='FeatureTable[Sequence]',default=False,
			help='specify representive sequences file before mask and align')
p.add_option('-o','--outdir',dest='outdir',metavar='[Directory]',default='./',
			help='specify the output directory')
p.add_option('-c','--cut',dest='cut',metavar='int',default=100,
			help='give the threshold of N of reads for otu you want to show')

(options,args) = p.parse_args()


os.system("if [ ! -d %s ];then mkdir %s;fi"%(options.outdir,options.outdir))


########select otu####


with open('tree.R', 'w') as rscript:
	print ('''
otu_table<-read.table("%s",header = T,skip=1,row.names = 1,check.names = F,stringsAsFactors = F,sep = "\\t",comment.char = "")
otu_table<-otu_table[,-dim(otu_table)[2]]
otusum<-colSums(t(otu_table))
out<-data.frame(rownames(otu_table)[otusum>as.numeric("%s")])
write.table(out,"%s/selected_features.txt",sep = "",row.names = F,col.names = F,quote = F)
'''
% (options.input,options.cut,options.outdir),
file = rscript)
os.system('Rscript tree.R')



#######align and mask#####
os.system("qiime alignment mafft   --i-sequences %s   --o-alignment %s/aligned-rep-seqs.qza&&\
	qiime alignment mask   --i-alignment %s/aligned-rep-seqs.qza   --o-masked-alignment %s/masked-aligned-rep-seqs.qza&&\
	qiime tools export  %s/masked-aligned-rep-seqs.qza --output-dir %s/"%(options.repseqs,options.outdir,options.outdir,options.outdir,options.outdir,options.outdir))


#######select rep-seqs####
with open('%s/selected_features_reseqs.fasta'%(options.outdir),'w') as fout:
	s_otuid=open('%s/selected_features.txt'%(options.outdir),'r')
	s_otuid=s_otuid.read()

	s_otuid=re.split('\n',s_otuid)
	sn=[]
	ln=1
	for line in open("%s/aligned-dna-sequences.fasta"%(options.outdir),'r'):
		
		line=re.sub('\n$','',line)
		line=re.sub('^>','',line)
		if ln in sn:
			fout.write(line+"\n")
		
		if line in s_otuid:
			fout.write('>'+line+"\n")
			sn.append(ln+1)
		ln+=1


######form tree#####
os.system("qiime tools import --input-path %s/selected_features_reseqs.fasta --output-path %s/selected_features_reseqs.qza --type 'FeatureData[AlignedSequence]'&&\
qiime phylogeny fasttree --i-alignment %s/selected_features_reseqs.qza --o-tree %s/selected_unrooted-tree.qza&&\
qiime phylogeny midpoint-root --i-tree %s/selected_unrooted-tree.qza --o-rooted-tree %s/selected_rooted-tree.qza&&\
qiime tools export %s/selected_rooted-tree.qza --output-dir %s/"%(options.outdir,options.outdir,options.outdir,options.outdir,options.outdir,options.outdir,options.outdir,options.outdir))



######visualize tree####
with open('tree.R', 'w') as rscript:
	print ('''
library("ggtree")
library("stringr")
otu_table<-read.table("%s",header = T,skip=1,row.names = 1,check.names = F,stringsAsFactors = F,sep = "\\t",comment.char = "")
tax<-otu_table[,dim(otu_table)[2]]
groupInfo <- split(rownames(otu_table), str_extract(tax,"p__[^;]{1,100}"))

otu_table<-otu_table[,-dim(otu_table)[2]]
metadata<-read.table("%s",header = T,row.names=1,check.names = F,stringsAsFactors = F,sep = "\\t",comment.char = "")
metadata<-metadata["%s"]
metagroup<-metadata[,1][match(colnames(otu_table),rownames(metadata))]
otu_table=scale(t(otu_table))

mysum=function(x){
	return(tapply(x,INDEX = metagroup,FUN = sum))
}


data=t(apply(otu_table,2,mysum))
tree <- read.tree("%s/tree.nwk")
tree<- groupOTU(tree, groupInfo)
p = ggtree(tree,aes(color=group))+geom_tiplab(size=0, align=TRUE, linesize=.1)
pdf(file="%s/%s", width=10, height=10)
gheatmap(p, data, offset = -0.2, width=1, hjust=-.3)+theme_tree2()+theme(legend.position = "right",text=element_text(size=15,face="bold"))
dev.off()
'''
% (options.input,options.metadata,options.group,options.outdir,options.outdir,options.group+'_PhylogeneticTreeHeatmap.pdf'),
file = rscript)

os.system('Rscript tree.R')


