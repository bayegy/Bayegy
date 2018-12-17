# coding:utf-8
import argparse
import re
import sys
import os
import pandas as pd
import math
from functools import reduce
# argument:
p = argparse.ArgumentParser(
    description="This script is used to form manifest of 16S demuxed sequences for the use of qiime2")
p.add_argument('-i', '--input', dest='input', metavar='<path>',
               help='Path of demuxed PE files')
p.add_argument('-m', '--map', dest='meta', metavar='<path>', default=False,
               help='raw sample metadata files. The sample ID must in first and last column. Two columns do not have to be same, but the IDs in first column must be respond to --sample-id-patern')
p.add_argument('-s', '--sample-id-patern', dest='sp', default=r'raw\.split\.(.+)\.[12]\.fq$', metavar='<regular expression>',
               help='The regular expression of sample ID in file names. You must use \(\) to contain sample ID expression')
p.add_argument('-f', '--forward-file-pattern', dest='fp', default=r'\.1\.fq$', metavar='<regular expression>',
               help='The regular expression representing forward sequence in file names')
p.add_argument('-r', '--reverse-file-pattern', dest='rp', default=r'\.2\.fq$', metavar='<regular expression>',
               help='The regular expression representing reverse sequence in file names')
p.add_argument('-n', '--groups', dest='number', default=3, metavar='{int,groups}',
               help='Number of groups, OR names of groups seprated by commas')
p.add_argument('-c', '--count-of-samples', dest='count', default=30, metavar='int',
               help='Count of samples')
p.add_argument('-d', '--delimmiter', dest='sep', default='', metavar='string',
               help='delimmiter between group and number in sample id')
p.add_argument('-o', '--output', dest='out', metavar='<directory>', default='./',
               help='The path of output files')
options = p.parse_args()

if not os.path.exists(options.out):
  os.makedirs(options.out)


sp = re.compile(options.sp)
fp = re.compile(options.fp)
rp = re.compile(options.rp)

pre_id = []
for root, dirs, files in os.walk(options.input):
  if len(dirs) == 0:
    root = os.path.abspath(root)
    for fl in files:
      if re.search(sp, fl):
        pre_id.append(re.search(sp, fl).group(1))


pre_id = list(set(pre_id))
pre_id.sort()

groups = ["C", "T1", "T2", "T3", "T4", "M1", "M2", "M3", "A", "B", "C", "D", "E", "F", "G", "H", "L"]


id_ds = {}
if options.meta:
  with open(options.meta, 'r') as meta:
    for line in enumerate(meta):
      li = re.split('\t', line[1].strip())
      if line[0] > 0:
        id_ds[pre_id[line[0] - 1]] = li[len(li) - 1].strip()
else:
  try:
    nog = int(options.number)
    group = groups[:nog]
  except Exception:
    group = re.split(',', options.number)
    group = [g.strip() for g in group]
    nog = len(group)
  try:
    nos = int(options.count)
    each = int(math.ceil(nos / nog))
    group1 = [a for a in group for j in range(each)]
    sample = [a + options.sep + str(j + 1) for a in group for j in range(each)]
  except Exception:
    def add(x, y):
      return x + y
    each = re.split(',', options.count)
    group1 = [[m] * int(n) for m, n in zip(group, each)]
    group1 = reduce(add, group1)
    sample = [list(range(1, int(sl) + 1)) for sl in each]
    sample = reduce(add, sample)
    sample = [j + options.sep + str(h) for j, h in zip(group1, sample)]
    nos = len(sample)

  mapp = pd.DataFrame({"#SampleID": sample[:nos], "Group": group1[:nos]})
  mapp.insert(2, "Description", sample[:nos])
  mapp.to_csv(options.out + '/' + 'sample-metadata.tsv', sep='\t', index=False)
  for spl in enumerate(sample[:nos]):
    id_ds[pre_id[spl[0]]] = spl[1]


fout = open(options.out + '/' + 'manifest.txt', 'w')
nfile = 0
ff = 0
#id_sets = []
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
#            id_sets.append(pre_id)
            nfile += 1
        except KeyError:
          ff += 1
fout.close()

if options.meta:
  with open(options.meta, 'r') as mp, open(options.out + '/' + 'sample-metadata.tsv', 'w') as outfile:
    for line in enumerate(mp):
      li = re.split('\t', line[1].strip())
      li = [l.strip() for l in li]
      fc = len(li) - 1
      if line[0] == 0:
        li[0] = "#SampleID"
        li[fc] = "Description"
        outfile.write('\t'.join(li) + '\n')
      else:
        li[0] = li[fc]
        outfile.write('\t'.join(li) + '\n')


emm = "\n    %s fastq files were writed" % (nfile)
emm1 = "\n    %s fastq files were filterd as the ids of samples were not found in mapping file" % (ff)
print(emm)
print(emm1)

# open(options.out + '/' + 'sample-metadata.tsv', 'w') as outmeta
#      li[0] = "#SampleID"
#      li[il] = "Description"
