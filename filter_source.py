#!/usr/bin/env python3
# coding:utf-8
import argparse
import re
import sys
import os
from qiime2 import Artifact
from pandas import *
# import pdb
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
p.add_argument('-f', '--filter', dest='filter', metavar='<flag:taxon1,taxon2>', default=False,
               help='Taxon to be filtered. use "keep:taxon1,taxon2" to keep noly taxon1, taxon2. or use "exclude:taxon1,taxon2" to filter taxon1, taxon2.')
p.add_argument('--min-abundance', dest='mina', metavar='int', default=1,
               help='pass this to filter out the otu table and representative sequences less than the specified value')
p.add_argument('-o', '--out', dest='out', metavar='<directory>',
               default='./', help='Specify the output directory')

opt = p.parse_args()
if not os.path.exists(opt.out):
    os.makedirs(opt.out)
opt.out = opt.out + '/'


def savepd(tb, name):
    tb.to_csv(opt.out + name, sep='\t', index=False)


def iterin(a, b):
    return list(map(lambda x: x in list(b), a))


def iterfind(a, b):
    out = []
    for i in a:
        found = False
        for j in b:
            if re.search(j, i, flags=re.IGNORECASE):
                found = True
                break
        out.append(found)
    return out


# load data
otu = Artifact.load(opt.otu)
otu = otu.view(DataFrame)
rep = Artifact.load(opt.rep)
rep = rep.view(Series)
taxon = Artifact.load(opt.taxon)
taxon = taxon.view(DataFrame)


if opt.map and opt.group:
    # load mapping file
    mapd = read_csv(opt.map, sep='\t')
    # filter mapping file
    mapd = mapd.ix[mapd[opt.group].notna(), :]
    mapd = mapd.ix[:, mapd.notna().sum() > 1]
    savepd(mapd, 'mapping_file.txt')
    # filter otu table
    otu = otu.ix[iterin(otu.index, mapd.ix[:, 0]), :]

otu = otu.ix[:, otu.sum() >= int(opt.mina)]
otu_set1 = otu.columns

# pdb.set_trace()

# filter taxa
if opt.filter:
    flag, ftaxa = opt.filter.split(':')
    ftaxa = ftaxa.split(',')
    if flag == "keep":
        taxon = taxon.ix[iterfind(taxon.ix[:, 0], ftaxa), :]
    else:
        taxon = taxon.ix[[not b for b in iterfind(taxon.ix[:, 0], ftaxa)], :]
otu_set2 = taxon.index

# calculate inner set
inner_otuset = [o for o in otu_set1 if o in otu_set2]


otu = otu.ix[:, iterin(otu.columns, inner_otuset)]
taxon = taxon.ix[iterin(taxon.index, inner_otuset), :]
rep = rep.ix[iterin(rep.index, inner_otuset)]
# save data
otu = Artifact.import_data("FeatureTable[Frequency]", otu)
otu.save(opt.out + 'table.qza')
rep = Artifact.import_data("FeatureData[Sequence]", rep)
rep.save(opt.out + 'rep-seqs.qza')
taxon = Artifact.import_data("FeatureData[Taxonomy]", taxon)
taxon.save(opt.out + 'taxonomy.qza')
