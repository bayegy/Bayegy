# -*- coding:utf-8 -*-
import re,sys,os
import argparse
class classfiletree(object):
	def __init__(self,father_filetree,inpath,outpath):
		if not os.path.exists(outpath):
			os.makedirs(outpath)
		self.outpath=os.path.abspath(outpath)
		self.inpath=inpath
		self.rawtree=self.outpath+'/raw_filetree.txt'
		self.alignedfiletree=self.outpath+'/aligned_file_tree.txt'
		self.child_filetree=self.outpath+'/child_filtree.txt'
		self.cleaned_filetree=self.outpath+'/cleaned_filetree.txt'
		self.father=father_filetree

	def get_tree(self):
		os.system("tree /F %s > %s/raw_filetree.txt"%(self.inpath,self.outpath))
	
	def align_filetree(self,raw_tree,length=100):
		n=int(length)
		fout=open(self.alignedfiletree,'w')
		for line in open(raw_tree):
			line=line.replace('\n','')
			if re.search('[a-zA-Z0-9]',line):
				line=line.ljust(n,'·')
			fout.write(line+'\n')
		fout.close()

	def inherit_filetree(self,father,mask):
		print(mask)
		mask=re.split(',',mask)
		mask=[re.compile(i) for i in mask]
		print(mask)
		ann={}
		for line in open(father,'r',encoding='utf-8'):
			line=re.sub('^[^0-9a-zA-Z_]+','',line.strip())
			line=re.sub('Between_[^_]+_and_[^_]+','',line)
			line=re.split('·+',line)
			for i in mask:
				line[0]=re.sub(i,'',line[0])
			#print(line[0])
			if not line[0]=="":
				try:
					if not (line[0] in ann.keys() or line[1]==""):
						ann[line[0]]=line[1]
				except IndexError:
					pass

		with open(self.child_filetree,'w') as fout:
			for line in open(self.alignedfiletree,'r'):
				line=line.strip()
				li=re.split('·+',line)[0]
				li=re.sub('^[^0-9a-zA-Z_]+','',li)
				li=re.sub('Between_[^_]+_and_[^_]+','',li)
				for i in mask:
					li=re.sub(i,'',li)
				#print(li)
				try:
					fout.write(line+ann[li]+'\n')
				except KeyError:
					fout.write(line+'\n')

	def cleanfiletree(self,dirfiletree):
		fout=open(self.cleaned_filetree,'w')
		l=[]
		n=1
		for line in open(dirfiletree):
			if n in l:
				if re.search('·+[^·]+\n$',line)!=None:
					fout.write(line)
					n=n+1
				elif re.search('[ |│]+$',line)!=None:
					l.append(n+1)
					n=n+1
			else:
				if re.search('·+[^·]+\n$',line)!=None:
					fout.write(line)
					n=n+1
				elif re.search('[ |│]+$',line)!=None:
					fout.write(line)
					l.append(n+1)
					n=n+1
					
		fout.close()

if __name__=="__main__":
	p = argparse.ArgumentParser(description="This script is used to form readme of results")
	p.add_argument('-i', '--input', dest = 'input', metavar = '<directory>',
				help = 'Results file path')
	p.add_argument('-o', '--output', dest = 'output', metavar = '<directory>', default = './',
				help = 'given an output directory')
	p.add_argument('-f', '--father-filetree', dest = 'father', metavar = '<path>',
				help = 'the filetree to be inherited')
	p.add_argument('-m', '--mask', dest = 'mask', metavar = '<str>',
				help = 'Strings in filename to be masked when aligning. Use comma as delimeter',default='Group\d*')
	options = p.parse_args()
	ex=classfiletree(father_filetree=options.father,inpath=options.input,outpath=options.output)
	ex.get_tree()
	ex.align_filetree(raw_tree=ex.rawtree)
	ex.inherit_filetree(father=ex.father,mask=options.mask)
	ex.cleanfiletree(dirfiletree=ex.child_filetree)
