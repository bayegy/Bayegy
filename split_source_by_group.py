# coding:utf-8
import argparse
import re
import sys
import os
from qiime2 import Artifact
from pandas import *
# argument:
p = argparse.ArgumentParser(description="This script is used to plot phylogentic tree of top abundant otus")
p.add_argument('-i', '--input', dest='otu', metavar='<path>', default=False,
               help='table.qza')
p.add_argument('-m', '--metadata', dest='map', metavar='<path>', default=False,
               help='Specify metadata with sample id at first column')
p.add_argument('-c', '--category', dest='group', metavar='<str>', default=False,
               help='Column name of group in mapping file, according to this group, the missing samples and zero features will be removed')
p.add_argument('-r', '--repseqs', dest='rep', metavar='<path>', default=False,
               help='representative seqs of features (.qza)')
p.add_argument('-t', '--taxonomy', dest='taxon', metavar='<path>', default=False,
               help='table.qza')
p.add_argument('--min-abundance', dest='mina', metavar='int', default=1,
               help='pass this to filter out the otu table and representative sequences less than the specified value')
p.add_argument('-o', '--out', dest='out', metavar='<directory>',
               default='./', help='Specify the output directory')

opt = p.parse_args()
if not os.path.exists(opt.out):
  os.makedirs(opt.out)
opt.out = opt.out + '/'
mapd = read_csv(opt.map, sep='\t')

mapd = mapd.ix[mapd[opt.group].notna(), :]


def savepd(tb, name):
  tb.to_csv(opt.out + name, sep='\t', index=False)


def iterin(a, b):
  return list(map(lambda x: x in list(b), a))


mapd = mapd.ix[:, mapd.notna().sum() > 2]
savepd(mapd, 'mapping_file.txt')

print("Min abundance is %d" % (int(opt.mina)))
otu = Artifact.load(opt.otu)
otu = otu.view(DataFrame)
otu = otu.ix[iterin(otu.index, mapd.ix[:, 0]), :]
otu = otu.ix[:, otu.sum() >= int(opt.mina)]

rep = Artifact.load(opt.rep)
rep = rep.view(Series)
rep = rep.ix[iterin(rep.index, otu.columns)]


taxon = Artifact.load(opt.taxon)
taxon = taxon.view(DataFrame)
taxon = taxon.ix[iterin(taxon.index, otu.columns), :]

otu = Artifact.import_data("FeatureTable[Frequency]", otu)
otu.save(opt.out + 'table.qza')


rep = Artifact.import_data("FeatureData[Sequence]", rep)
rep.save(opt.out + 'rep-seqs.qza')

taxon = Artifact.import_data("FeatureData[Taxonomy]", taxon)
taxon.save(opt.out + 'taxonomy.qza')
