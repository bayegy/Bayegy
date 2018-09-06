# -*- coding: utf-8 -*-

from __future__ import print_function
from optparse import OptionParser
import re,sys,os

#argument:
usage = '%prog -[i]'
p = OptionParser(usage = usage)
p.add_option('-m', '--metadata', dest = 'metadata', metavar = '*.txt',
			help = 'taxonomic count biom file')
p.add_option('-g', '--group', dest = 'group', metavar = 'Group',
			help = 'Column name of group in metadata')
p.add_option('-l', '--level', dest = 'level', metavar = 'level', default = 'Genus',
			help = 'Specify the taxonomy level')

(options,args) = p.parse_args()

with open(options.metadata,'r') as infile:
	ln=1
	g=[]

	for line in infile:
		line=re.sub('\n$','',line)
		#print(line)
		line=re.split('\t',line)
		if ln==1:
			for i in range(len(line)):
				if line[i].lstrip().rstrip()==options.group:
					gn=i
					ln+=1
		else:
			g.append(line[gn])

unig=list(set(g))
unig.sort()
if '' in unig:
	unig.remove("")

for i in range(len(unig)):
	for j in range(len(unig)):
		if i<j:
			print(unig[i])
			print(unig[j])
			os.system('differential_abundance.py -i exported/DiffAbundance/tax/otu_table.even_%s.taxonomy.biom -o exported/DiffAbundance/DESeq2_%s_Between_%s_and_%s_DiffAbundance_%s.txt  -a DESeq2_nbinom -m %s -c %s -x %s -y %s -d'%(options.level,options.group,unig[i],unig[j],options.level,options.metadata,options.group,unig[i],unig[j]))
