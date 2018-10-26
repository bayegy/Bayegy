# coding:utf-8
import argparse
import re
import sys
import os
import pandas as pd

# argument:
p = argparse.ArgumentParser(
    description="This script is used to form manifest of 16S demuxed sequences for the use of qiime2")
p.add_argument('-i', '--input', dest='input', metavar='<path>',
               help='Path of demuxed PE files')
p.add_argument('-m', '--map', dest='meta', metavar='<path>',
               help='raw sample metadata files. The sample ID must in first and last column. Two columns do not have to be same, but the IDs in first column must be respond to --sample-id-patern')
p.add_argument('-s', '--sample-id-patern', dest='sp', default=r'raw\.split\.(.+)\.[12]\.fq$', metavar='<regular expression>',
               help='The regular expression of sample ID in file names. You must use \(\) to contain sample ID expression')
p.add_argument('-f', '--forward-file-pattern', dest='fp', default=r'\.1\.fq$', metavar='<regular expression>',
               help='The regular expression representing forward sequence in file names')
p.add_argument('-r', '--reverse-file-pattern', dest='rp', default=r'\.2\.fq$', metavar='<regular expression>',
               help='The regular expression representing reverse sequence in file names')
p.add_argument('-o', '--output', dest='out', metavar='<directory>', default='./',
               help='The path of output files')
options = p.parse_args()

if not os.path.exists(options.out):
  os.makedirs(options.out)


for root, dirs, files in os.walk(options.input):
  if len(dirs) == 0:
    root = os.path.abspath(root)
    for fl in files:
      if re.search(sp, fl):
        try:
          info = "%s/%s" % (root, fl)
          pre_id = re.search(sp, fl).group(1)


id_ds = {}
with open(options.meta, 'r') as meta:
  for line in enumerate(meta):
    li = re.split('\t', line[1].strip())
    if line[0] > 0:
      id_ds[li[0].strip()] = li[len(li) - 1].strip()


print("############################################################writting manifest\n\nThe regular expression for matching sample ID is %s, you should change -s if no fastq files were writed" % (options.sp))
print("\nThe regular expression for matching forward sequences is %s" % (options.fp))
print("\nThe regular expression for matching reverse sequence is %s \nyou need to change -r  and -f if the direction of sequence is wrong in manifest" % (options.rp))
sp = re.compile(options.sp)
fp = re.compile(options.fp)
rp = re.compile(options.rp)

fout = open(options.out + '/' + 'manifest.txt', 'w')
nfile = 0
ff = 0
id_sets = []
fout.write('sample-id,absolute-filepath,direction\n')
for root, dirs, files in os.walk(options.input):
  if len(dirs) == 0:
    root = os.path.abspath(root)
    for fl in files:
      if re.search(sp, fl):
        try:
          info = "%s/%s" % (root, fl)
          pre_id = re.search(sp, fl).group(1)
          sn = id_ds[pre_id]
          if re.search(rp, fl):
            fout.write("%s,%s,reverse\n" % (sn, info))
            # id_sets.append(pre_id)
            nfile += 1
          elif re.search(fp, fl):
            fout.write("%s,%s,forward\n" % (sn, info))
            id_sets.append(pre_id)
            nfile += 1
        except KeyError:
          ff += 1
fout.close()


with open(options.meta, 'r') as mp, open(options.out + '/' + 'sample-metadata.tsv', 'w') as outfile:
  miss = 0
  for line in enumerate(mp):
    li = re.split('\t', line[1].strip())
    li = [l.strip() for l in li]
    fc = len(li) - 1
    if line[0] == 0:
      li[0] = "#SampleID"
      li[fc] = "Description"
      outfile.write('\t'.join(li) + '\n')
    else:
      if li[0] in id_sets:
        outfile.write('\t'.join(li) + '\n')
      else:
        miss += 1


emm = "\n    %s fastq files were writed" % (nfile)
emm1 = "\n    %s fastq files were filterd as the ids of samples were not found in mapping file" % (ff)
print(emm)
print(emm1)
print('\n    %s samples were removed from sample-metadata, as the ids of samples were not found in forward fastq files\' names' % (miss))

# open(options.out + '/' + 'sample-metadata.tsv', 'w') as outmeta
#      li[0] = "#SampleID"
#      li[il] = "Description"
