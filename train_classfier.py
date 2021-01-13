#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import argparse
import re
import sys
import os

p = argparse.ArgumentParser()

p.add_argument('-i', '--input', dest='input', metavar='<path>', help='Reference set')
p.add_argument('-t', '--taxonomy', dest='taxonomy', metavar='<path>', help='Reference taxonomy')
p.add_argument('-f', '--forward-primer', dest='fp', metavar='<str>', default=False, help='Forward primer')
p.add_argument('-r', '--reverse-primer', dest='rp', metavar='<str>', default=False, help='Reverse primer')
p.add_argument('-o', '--output', dest='output', metavar='<directory>', default='./', help='Given an output directory')
p.add_argument('-p', '--prefix', dest='prefix', metavar='<str>', default='', help='Set outfiles prefix')


options = p.parse_args()
if not os.path.exists(options.output):
    os.makedirs(options.output)

options.output = options.output + '/' + options.prefix

if not options.fp or not options.rp:
    os.system('''qiime tools import --type 'FeatureData[Sequence]'   --input-path %s  --output-path %srep-set.qza&&\
qiime tools import  --type 'FeatureData[Taxonomy]'  --input-format HeaderlessTSVTaxonomyFormat  --input-path %s --output-path %sref-taxonomy.qza&&\
qiime feature-classifier fit-classifier-naive-bayes --i-reference-reads %srep-set.qza --i-reference-taxonomy %sref-taxonomy.qza --o-classifier %sclassifier.qza''' %
              (options.input, options.output, options.taxonomy, options.output,
               options.output, options.output, options.output))
elif options.taxonomy:
    os.system('''qiime tools import --type 'FeatureData[Sequence]'   --input-path %s  --output-path %srep-set.qza&&\
qiime tools import  --type 'FeatureData[Taxonomy]'  --input-format HeaderlessTSVTaxonomyFormat  --input-path %s --output-path %sref-taxonomy.qza&&\
qiime feature-classifier extract-reads  --i-sequences %srep-set.qza  --p-f-primer %s   --p-r-primer %s --o-reads %sref-seqs.qza&&\
qiime feature-classifier fit-classifier-naive-bayes --i-reference-reads %sref-seqs.qza --i-reference-taxonomy %sref-taxonomy.qza --o-classifier %sclassifier.qza''' %
              (options.input, options.output, options.taxonomy, options.output,
               options.output, options.fp, options.rp, options.output, options.output, options.output, options.output))
else:
    os.system('''qiime tools import --type 'FeatureData[Sequence]'   --input-path %s  --output-path %srep-set.qza&&\
qiime feature-classifier extract-reads  --i-sequences %srep-set.qza  --p-f-primer %s   --p-r-primer %s --o-reads %sref-seqs.qza''' %
              (options.input, options.output,
               options.output, options.fp, options.rp, options.output))
