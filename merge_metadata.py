import sys,re,os

with open('Rscript.R','w') as Rscript:
	print('''
map=read.table('%s',header = T,row.names = 1,stringsAsFactors = F,check.names = F,comment.char = "",sep = "\\t")
meta=read.table('%s',header = T,row.names = 1,stringsAsFactors = F,check.names = F,comment.char = "",sep = "\\t")
rowname_join<-function(x,y)
{
  ya<-data.frame(y[match(rownames(x),rownames(y)),])
  colnames(ya)<-colnames(y)
  colnames(ya)[colnames(ya)%%in%%colnames(x)]<-paste(colnames(ya)[colnames(ya)%%in%%colnames(x)],"_y",sep = "")
  out<-data.frame(x,ya,check.rows = T,check.names = F)
  return(out)
}
out<-rowname_join(map,meta)
out<-data.frame(SMAPLEID=rownames(out),out)
colnames(out)[1]<-"#SampleID"
write.table(out,'metadata.tsv',sep = '\\t',quote = F,row.names=F)
'''%(sys.argv[1],sys.argv[2]),file=Rscript)

os.system('Rscript Rscript.R')