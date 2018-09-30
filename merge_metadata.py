#coding='utf-8'
import sys
import re
import os

with open('Rscript.R', 'w') as Rscript:
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
out<-data.frame(SMAPLEID=rownames(out),out,check.names=F)
out$description<-out$SMAPLEID
colnames(out)[1]<-"#SampleID"
out<-out[out$Farm!='F1',]
write.table(out,'pre_map.txt',sep = '\\t',quote = F,row.names=F,na='')
''' % (sys.argv[1], sys.argv[2]), file=Rscript)

os.system('Rscript Rscript.R')
