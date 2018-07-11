
import sys,os,re

n=int(sys.argv[1])
info = os.getcwd() 
path=info+'/'+sys.argv[2]

fout=open(info+'/aligned_file_tree.txt','w')

for line in open(path):
	line=line.replace('\n','')
	if re.search('[a-zA-Z0-9]',line)!=None:
		line=line.ljust(n,'-')
	fout.write(line+'\n')

fout.close()