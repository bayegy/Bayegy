#!/usr/bin/env python3
# coding:utf-8
import argparse
import re
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
p.add_argument('--raw', metavar='<directory>',
               default=None, help='output csv tables before filter')
p.add_argument('--tsv', metavar='<path>',
               default=None, help='output otu table with taxonomy at last column')
p.add_argument('--consensus', metavar='<path>',
               default=None, help='output otu table with Consensus Lineage at last column')
p.add_argument('--even', metavar='<path>',
               default=None, help='output even otu table with taxonomy at last column')
p.add_argument('--evenconsensus', metavar='<path>',
               default=None, help='output even otu table with Consensus Lineage at last column')
p.add_argument('--fmap', metavar='<path>',
               default=None, help='output mapping file after filter')
p.add_argument('--odepth', metavar='<path>',
               default=None, help='Output depth used in even tabale calculation')
opt = p.parse_args()
if not os.path.exists(opt.out):
    os.makedirs(opt.out)
opt.out = opt.out + '/'


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


def clean_tax(tax):
    for ambi in ['uncultured', 'Ambiguous', 'unidentified', 'unassigned', 'unclassified']:
        tax = re.sub(';[^;]*' + ambi + '.*', '', tax, flags=re.I)
    tax = re.sub('; ', ';', tax)
    tax = re.sub(r'[\-| ]+', '_', tax)
    tax = re.sub(';$', '', tax)
    return(tax)


def ch_col(df, index, val):
    columns = list(df.columns)
    columns[index] = val
    df.columns = columns
    return df


# load otu table from qiime2 Artifact
if opt.otu.endswith('.qza'):
    otu = Artifact.load(opt.otu)
    otu = otu.view(DataFrame)
    taxon = Artifact.load(opt.taxon)
    taxon = taxon.view(DataFrame)
# otherwise load otu table from tsv table
else:
    skip_rows = []
    with open(opt.otu) as f:
        for num, line in enumerate(f):
            if len(line.split('\t')) < 2:
                skip_rows.append(num)
    otu = read_csv(opt.otu, sep='\t', skiprows=skip_rows, index_col=0)
    otu_t = otu.iloc[:, -1]
    taxon = DataFrame({
        'Taxon': otu_t.copy(),
        'Confidence': [0.9] * otu_t.shape[0]
    }, index=otu_t.index.copy())
    taxon.index.set_names("Feature ID", inplace=True)
    otu.drop(otu.columns[-1], axis=1, inplace=True)
    otu = otu.T
    otu.columns.set_names(None, inplace=True)

if opt.raw:
    otu.T.to_csv(os.path.join(opt.raw, "otu_tale_before_filter.csv"))
    taxon.to_csv(os.path.join(opt.raw, "taxonomy_before_filter.csv"))

# load rep seqs from qiime2 Artifact
if opt.rep.endswith('.qza'):
    rep = Artifact.load(opt.rep)
# otherwise import rep seqs from fasta
else:
    rep = Artifact.import_data('FeatureData[Sequence]', opt.rep)
rep = rep.view(Series)

if opt.map and opt.group:
    # load mapping file
    mapd = read_csv(opt.map, sep='\t')
    # filter mapping file
    mapd = mapd.loc[mapd[opt.group].notna(), :]
    mapd = mapd.loc[:, mapd.notna().sum() > 1]
    mapd.to_csv(
        opt.fmap if opt.fmap else os.path.join(opt.out, 'mapping_file.txt'),
        sep='\t', index=False
    )
    # filter otu table
    otu = otu.loc[iterin(otu.index, mapd.iloc[:, 0]), :]

otu = otu.loc[:, otu.sum() >= int(opt.mina)]
otu_set1 = otu.columns

# pdb.set_trace()

# filter taxa
if opt.filter:
    flag, ftaxa = opt.filter.split(':')
    ftaxa = ftaxa.split(',')
    if flag == "keep":
        taxon = taxon.loc[iterfind(taxon.iloc[:, 0], ftaxa), :]
    else:
        taxon = taxon.loc[[not b for b in iterfind(taxon.iloc[:, 0], ftaxa)], :]

# clean taxa
taxon.iloc[:, 0] = [clean_tax(t) for t in taxon.iloc[:, 0]]

otu_set2 = taxon.index

# calculate inner set
inner_otuset = [o for o in otu_set1 if o in otu_set2]

# loc is much more faster
# otu = otu.iloc[:, iterin(otu.columns, inner_otuset)]
otu = otu.loc[:, inner_otuset]
# taxon = taxon.iloc[iterin(taxon.index, inner_otuset), :]
taxon = taxon.loc[inner_otuset, :]
# rep = rep.iloc[iterin(rep.index, inner_otuset)]
rep = rep.loc[inner_otuset]
# save data
if opt.tsv or opt.consensus:
    df = otu.T
    df = df.join(taxon)
    df.drop(df.columns[-1], axis=1, inplace=True)
    df.index.set_names("#OTU ID", inplace=True)
    df = ch_col(df, -1, "taxonomy")
    # print(df)
    df.to_csv(opt.tsv, sep='\t')
    if opt.consensus:
        df = ch_col(df, -1, "Consensus Lineage")
        df.to_csv(opt.consensus, sep='\t')

otu_qza = Artifact.import_data("FeatureTable[Frequency]", otu)
otu_qza.save(opt.out + 'table.qza')

if opt.even or opt.evenconsensus or opt.odepth:
    from qiime2.plugins import feature_table
    depth = int(otu.T.sum().min())
    # if depth < 5000:
    #     print("Min depth among samples is too small({}); reset sample depth to: 5000".format(depth))
    #     depth = 5000
    if opt.odepth:
        with open(opt.odepth, 'w') as f:
            f.write(str(depth))
    print("sample depth is: {}".format(depth))
    if opt.even or opt.evenconsensus:
        rarefy_result = feature_table.methods.rarefy(table=otu_qza, sampling_depth=depth)
        df = rarefy_result.rarefied_table.view(DataFrame).T
        df = df.join(taxon)
        df.drop(df.columns[-1], axis=1, inplace=True)
        df.index.set_names("#OTU ID", inplace=True)
        if opt.even:
            df = ch_col(df, -1, "taxonomy")
            df.to_csv(opt.even, sep='\t')
        if opt.evenconsensus:
            df = ch_col(df, -1, "Consensus Lineage")
            df.to_csv(opt.evenconsensus, sep='\t')

rep = Artifact.import_data("FeatureData[Sequence]", rep)
rep.save(opt.out + 'rep-seqs.qza')
taxon = Artifact.import_data("FeatureData[Taxonomy]", taxon)
taxon.save(opt.out + 'taxonomy.qza')
