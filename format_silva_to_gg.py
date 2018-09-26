# -*- coding: utf-8 -*-
import argparse
import re
import sys
import os

p = argparse.ArgumentParser(
    description="This script is used to transform the taxonomy derived from Silva to GreenGene")
p.add_argument('-i', '--input', dest='input', metavar='<path>', help='Taxonomic qza')
p.add_argument('-o', '--output', dest='output', metavar='<directory>', default='./', help='Given an output directory')
options = p.parse_args()

if not os.path.exists(options.output):
  os.makedirs(options.output)


def trans(tar):
  level = {"D_0": "k", "D_1": "p", "D_2": "c", "D_3": "o", "D_4": "f", "D_5": "g", "D_6": "s"}
  for t in level.keys():
    tar = re.sub(t, level[t], tar)
  tar = re.sub(' ', '-', tar)
  tar = re.sub(';D_7.+', '', tar)
  tar = re.sub(';[^;]+uncultured.+', '', tar)
  return(tar)


os.system('qiime tools export --input-path %s --output-path %s' % (options.input, options.output))

with open(options.output + '/' + 'taxonomy.tsv', 'r') as infile, open(options.output + '/' + 'transformed_taxonomy.tsv', 'w') as outfile:
  for line in enumerate(infile):
    if line[0] == 0:
      outfile.write(line[1])
    else:
      line = re.split('\t', line[1].strip())
      line[1] = trans(line[1])
      outfile.write('\t'.join(line) + '\n')


os.system('qiime tools import --input-format TSVTaxonomyFormat --output-path %s/%s --input-path %s/transformed_taxonomy.tsv --type FeatureData[Taxonomy]' % (
    options.output, options.input.split("/")[-1], options.output))
