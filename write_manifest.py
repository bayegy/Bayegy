from __future__ import print_function
from optparse import OptionParser
import re,sys,os

#*********************************************************************** 青年才俊红凡凡 *********************************************************************************
#argument:
usage = '%prog -[i]'
p = OptionParser(usage = usage)
p.add_option('-i', '--input', dest = 'input', metavar = '<file>',
			help = 'Path of demuxed PE files')
p.add_option('-m', '--map', dest = 'meta', metavar = '<File>',
			help = 'raw sample metadata files')
p.add_option('-s', '--sample-id-patern', dest = 'sp', default=r'\.([^\.]+)\.[12]\.fq$',metavar = '<str>',
			help = 'column name in sample-metadata file')
p.add_option('-f', '--forward-file-pattern', dest = 'fp',default = r'\.1\.fq$', metavar = '<int>',
			help = 'specify how many species to be display, defaulf is 15')
p.add_option('-r', '--reverse-file-pattern', dest = 'rp',default = r'\.2\.fq$', metavar = '<int>',
			help = 'specify how many species to be display, defaulf is 15')
p.add_option('-o', '--output', dest = 'out', metavar = '<str>', default = 'none',
			help = 'specify numeric variables excluded from rda seprated by commas,use "none" if all numeric variables is expected')
(options,args) = p.parse_args()

if not os.path.exists(options.out):
	os.makedirs(options.out)

id_ds={}
with open(options.meta,'r') as meta,open(options.out+'/'+'sample-metadata.tsv','w') as outmeta:
		for line in enumerate(meta):
			li=re.split('\t',re.sub('\n$','',line[1]))
			il=len(li)-1
			if line[0]==0:
				li[0]="#SampleID"
				li[il]="Description"
				outmeta.write('\t'.join(li)+'\n')
			else:
				id_ds[li[0]]=li[il]
				li[0]=li[il]
				outmeta.write('\t'.join(li)+'\n')




sp=re.compile(options.sp)
fp=re.compile(options.fp)
rp=re.compile(options.rp)

fout = open(options.out+'/'+'manifest.txt', 'w') # 合并内容到该文件 
nfile=0
fout.write('sample-id,absolute-filepath,direction\n')
for root, dirs, files in os.walk(options.input): 
	if len(dirs)==0:
		root=os.path.abspath(root) 
		for fl in files: 
			if re.search(sp,fl):
				info = "%s/%s" % (root,fl) 
				sn = id_ds[re.search(sp,fl).group(1)]
				if re.search(rp,fl):
					fout.write("%s,%s,reverse\n"%(sn,info))
					nfile+=1
				elif re.search(fp,fl):
					fout.write("%s,%s,forward\n"%(sn,info))
					nfile+=1
					
fout.close() 
emm = "    %s fastq file were writed" % (nfile) 
print(emm)


