#!/usr/bin/python
# -*- coding: utf-8 -*-
from __future__ import print_function,division
from optparse import OptionParser
from collections import defaultdict
import re,sys,os

#*********************************************************************** 青年才俊红凡凡 *********************************************************************************
#argument:
usage = '%prog -[i]'
p = OptionParser(usage = usage)
p.add_option('-i','--input',dest='input',metavar='*.absoult.xls',
			help='taxonomic count data file')
p.add_option('-o','--output',dest='output',metavar='*.heatmap.pdf',default='permanova.pdf',
			help='given an output file name, default is permanova.pdf')
p.add_option('-g','--group',dest='group',metavar='group',default=False,
			help='annotate the columns with the same color from this group file,default is null')
p.add_option('-e', '--env', dest = 'env', metavar = 'env.list',
			help = 'environment parameter file')
p.add_option('-n', '--number', dest = 'number', metavar = 'top number or percent number', default = '0.01',
			help = 'show the number or percent you need, default is relative abundance more than 1%')
p.add_option('-m', '--model', dest = 'model', metavar = 't or p', action = 'store', default = 'p',
			help = "choose species according to top of rank abundance or percent of rank abundance,'t'=top and 'p'=percent default is percent")
p.add_option('-s', '--show', dest = 'show', metavar = 'True or False', action = 'store_true', default = False,
			help="push those numbers in all samples below this threshold into 'Others' default is False")
p.add_option('--cexRow',dest='cexRow',default='0.6',
			help='positive numbers, used for the row axis labeling. default is 0.6')
p.add_option('--cexCol',dest='cexCol',default='1.1',
			help='positive numbers, used for the column axis labeling. default is 0.8')
p.add_option('--srtCol',dest='srtCol',default='315',
			help='angle of column labels, in degrees from horizontal')
p.add_option('--ncol',dest='ncol',default='2',
			help='The number of columns in which to set the legend items, if group. default is 1')
p.add_option('--width',dest='width',default='10',
			help='the width of the graphics region in inches.  The default values are 10.')
p.add_option('--height',dest='height',default='10',
			help='the height of the graphics region in inches.  The default values are 10.')
(options,args) = p.parse_args()

if not (options.input and options.env):
	p.error("must have argument -i -e")
	sys.exit()
else:
	pass

#********************************************************************* I AM A LONELY LINE ********************************************************************************
#data file:
infile = open(options.input)
outfile = open('data.txt', 'w')
line_dict = {}
total_dict = defaultdict(int)
percent_dict = defaultdict(float)
def tree():
	return defaultdict(tree)
items = tree()
all = 0
n = 1
spe = []
for i in infile:
	isplit = re.split(r'\t',i.rstrip())
	if n == 1:
		rawsample = re.split(r'\t',i.rstrip())
		print ('\t%s' % '\t'.join(isplit[1:-1]), file = outfile)
	else:
		if isplit[0] in spe:
			asplit = re.split(r'\t', line_dict[isplit[0]])
			temp = isplit[0] + '\t'
			for a in rawsample[1:-1]:
				tp = int(asplit[rawsample.index(a)]) + int(isplit[rawsample.index(a)])
				temp = temp + str(tp) + '\t'
			line_dict[isplit[0]] = temp.rstrip()
		else:
			spe.append(isplit[0])
			line_dict[isplit[0]] = '\t'.join(isplit[:-1])
		for a in rawsample[1:-1]:
			total_dict[isplit[0]] += int(isplit[rawsample.index(a)])
		all += total_dict[isplit[0]]
	n += 1
infile.close()

if options.model == 'p':
	for i in line_dict.keys():
		percent_dict[i] = float(total_dict[i] / all)
	if options.show:
		items['Others'] = defaultdict(int)
		for i in percent_dict.keys():
			isplit = re.split(r'\t',line_dict[i])
			if percent_dict[i] >= float(options.number) and i != 'Others':
				print ('%s' % line_dict[i],file = outfile)
			else:
				for a in rawsample[1:-1]:
					items['Others'][a] += int(isplit[rawsample.index(a)])
		tmp = 'Others' + '\t'
		tmptotal = 0
		for b in rawsample[1:-1]:
			tmp += str(items['Others'][b]) + '\t'
			tmptotal += int(items['Others'][b])
		if tmptotal != 0:
			print ('%s' % tmp.rstrip(), file = outfile)
		else:
			print ('NOTICE: no species relative abundance less than %s' % options.number)
	else:
		for i in percent_dict.keys():
			if percent_dict[i] >= float(options.number):
				print ('%s' % line_dict[i],file = outfile)

else:
	glist = []
	x = 0
	xx = 0
	dict_sort = sorted(total_dict.items(),key = lambda d:d[1],reverse = True)
	for k,v in dict_sort:
		glist.append(k)
	leng = len(glist)
	if leng < int(options.number):
		print ('line numbers less than %s, used all' % options.number)
		while x <= leng - 1:
			print ('%s' % line_dict[glist[x]], file = outfile)
			x += 1
	else:
		while xx <= int(options.number) - 1:
			if glist[xx] != 'Others':
				print ('%s' % line_dict[glist[xx]],file = outfile)
			else:
				print ('There is Others')
			xx += 1
		l = int(options.number)
		items['Others'] = defaultdict(int)
		asplit = re.split(r'\t',line_dict['Others'])
		while l <= leng - 1:
			isplit = re.split(r'\t',line_dict[glist[l]])
			for a in rawsample[1:-1]:
				items['Others'][a] += int(isplit[rawsample.index(a)])
			l += 1
		try:
			for a in rawsample[1:-1]:
				items['Others'][a] += int(asplit[rawsample.index(a)])
		except:
			pass
		tmp = 'Others' + '\t'
		for b in rawsample[1:-1]:
			tmp += str(items['Others'][b]) + '\t'
		if options.show:
			print ('%s' % tmp.rstrip(),file = outfile)
		else:
			pass
outfile.close()

samgroup = []
groupdict = defaultdict(list)
if options.group:
	with open(options.group) as ingroup:
		for i in ingroup:
			isplit = i.rstrip().split('\t')
			if isplit[0] == 'sample' or re.match(r'#', isplit[0]):
				pass
			else:
				samgroup.append(isplit[1])
				groupdict[isplit[1]].append(isplit[0])
	
	group = set(samgroup)

	with open('spe.txt', 'w') as outfile:
		with open('data.txt') as infile:
			for i in infile:
				isplit = i.rstrip().split('\t')
				if isplit[0] == '':
					rawsample = i.rstrip().split('\t')
					print ('species\tgroup', file = outfile)
				else:
					gtotal = []
					gdict = {}
					for a in group:
						total = 0
						for b in groupdict[a]:
							total += int(isplit[rawsample.index(b)])
						gtotal.append(total)
						gdict[total] = a
					value = sorted(gtotal, reverse = True)[0]
					print ('%s\t%s' % (isplit[0], gdict[value]), file = outfile)

#********************************************************************* I AM A LONELY LINE ********************************************************************************
#Rscript:
path = os.getcwd()
with open('permanova.R', 'w') as rscript:
	print ('''
library(vegan)
library(gplots)
dat <- read.table("%s/data.txt", header = TRUE, row.names = 1, sep = "\\t")
env <- read.table("%s/%s", header = TRUE, row.names = 1, sep = "\\t")'''
% (path, path, options.env),
file = rscript)
	if options.group:
		print ('''
groups <- read.table("%s/spe.txt", header = TRUE, sep = "\\t")
rownames(groups) <- groups$species
group <- unique(groups$group)
cols <- c("red", "green", "blue", "purple", "cyan", "darkblue", "seagreen", "steelblue","DarkTurquoise", "Sienna", "Chocolate", "BlueViolet", "Magenta", "brown", "gray", "darkred", "pink", "orange")
col <- rep("black", dim(dat)[1])
for (i in 1:length(group)){
	col[groups[rownames(dat), 2] == group[i]] <- cols[i]
}''' % path,
file = rscript)
	else:
		print ('warning: no group file input')
		print ('col <- rep("black", dim(dat)[1])', file = rscript)
	print ('''
y <- t(env)
b <- matrix(0, nrow = nrow(dat), ncol = nrow(y))
cor <- matrix(0, nrow = nrow(dat), ncol = nrow(y))
colnames(cor) <- colnames(b) <- rownames(y)
rownames(cor) <- rownames(b) <- rownames(dat)
for (a in 1:nrow(dat)){
	d <- dat[a,] + 0.0001
	corx <- as.numeric(dat[a,])
	x <- as.matrix(d)
	x <- t(x)
	braydist <- vegdist(x, method = "bray")
	for (i in 1:nrow(y)){
		permanResult <- adonis(braydist ~ y[i, ], permutations = 1000)
		b[a, i] = permanResult$aov.tab$`Pr(>F)`[1]
		cory <- as.numeric(y[i,])
		correlation <- cor.test(corx, cory, method = "pearson")
		cor[a, i] <- correlation$estimate
	}
}
write.table(
	cor,
	file = "%s/correlation.txt",
	quote = FALSE,
	col.names = TRUE,
	row.names = TRUE,
	sep = "\\t")
write.table(
	b,
	file = "%s/p_value.txt",
	quote = FALSE,
	col.names = TRUE,
	row.names = TRUE,
	sep = "\\t")
cell <- b
for(i in 1:dim(cell)[1]){
	for(j in 1:dim(cell)[2]){
		if(cell[i,j] >= 0.01 & cell[i,j]< 0.05){
			cell[i,j] <- "+"
		}
		else if(cell[i,j] < 0.01){
			cell[i,j] <- "*"
		}
		else{
			cell[i,j] <- ""
		}
	}
}
pdf("%s/%s", height = %s, width = %s)
cexRow <- %s
cexCol <- %s
srtCol <- %s
ncol <- %s
heatmap.2(
	cor,
	col = bluered, 
	scale = "none", 
    trace = "none", 
	offsetRow = 0,
	key.title = "Spearman correlation",
	key.xlab = "",
	key.ylab = "",
	key.ytickfun = NULL,
	keysize = 0.9,
	cellnote = cell,
	notecol = "black",
	colRow = col,
	cexRow = cexRow,
	cexCol = cexCol,
	srtCol = srtCol,
	offsetCol = 0,
	adjCol = c(0,1),
	margins = c(5, 10)
	)
	'''
% (path, path, path, options.output, options.height, options.width, options.cexRow, options.cexCol, options.srtCol, options.ncol),
file = rscript)
	if options.group:
		print ('''
legend(
	x = "topright",
	legend = group,
	col = cols[1:length(group)],
	pch = 15,
	horiz = FALSE,
	bty = "n",
	xpd = TRUE,
	ncol = ncol,
	cex = 0.6
	)
dev.off()''',
file = rscript)
	else:
		print ('dev.off()', file = rscript)

os.system('Rscript permanova.R')
