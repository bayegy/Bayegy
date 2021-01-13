# -*- coding: utf-8 -*-
import argparse
import re
import sys
import os
# argument:
p = argparse.ArgumentParser()
p.add_argument('-i', '--input', dest='input', metavar='<path>',
               help='Path of abundance biom file')
p.add_argument('-m', '--metadata', dest='metadata', metavar='*.txt',
               help='taxonomic count biom file')
p.add_argument('-g', '--group', dest='group', metavar='Group',
               help='Column name of group in metadata')
p.add_argument('-l', '--level', dest='level', metavar='level', default='Genus',
               help='Specify the taxonomy level')
p.add_argument('-o', '--output', dest='out', metavar='<directory>', default='./',
               help='The path of output files')


options = p.parse_args()
if not os.path.exists(options.out):
    os.makedirs(options.out)


with open(options.metadata, 'r') as infile:
    ln = 1
    g = []

    for line in infile:
        line = re.sub('\n$', '', line)
        # print(line)
        line = re.split('\t', line)
        if ln == 1:
            for i in range(len(line)):
                if line[i].lstrip().rstrip() == options.group:
                    gn = i
                    ln += 1
        else:
            g.append(line[gn])

unig = list(set(g))
unig.sort()
if '' in unig:
    unig.remove("")


for i in range(len(unig)):
    for j in range(len(unig)):
        if i < j:
            print("Comparing %s to %s" % (unig[i], unig[j]))
            os.system("differential_abundance.py -i %s -o %s/DESeq2_%s_Between_%s_and_%s_DiffAbundance_%s.txt  -a DESeq2_nbinom -m %s -c %s -x '%s' -y '%s' -d" %
                      (options.input, options.out, options.group, unig[i], unig[j], options.level, options.metadata, options.group, unig[i], unig[j]))
