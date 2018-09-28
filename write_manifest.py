# coding:utf-8
import argparse
import re
import sys
import os

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

id_ds = {}
with open(options.meta, 'r') as meta, open(options.out + '/' + 'sample-metadata.tsv', 'w') as outmeta:
  for line in enumerate(meta):
    li = re.split('\t', re.sub('\n$', '', line[1]))
    il = len(li) - 1
    if line[0] == 0:
      li[0] = "#SampleID"
      li[il] = "Description"
      outmeta.write('\t'.join(li) + '\n')
    else:
      id_ds[li[0]] = li[il]
      li[0] = li[il]
      outmeta.write('\t'.join(li) + '\n')

print("########################################################################################writting manifest\n\nThe regular expression for matching sample ID is %s, you should change -s if no fastq files were writed" % (options.sp))
print("\nThe regular expression for matching reverse sequences is %s" % (options.fp))
print("\nThe regular expression for matching forward sequence is %s, you need to change -r  and -f if the direction of sequence is wrong in manifest" % (options.rp))
sp = re.compile(options.sp)
fp = re.compile(options.fp)
rp = re.compile(options.rp)

fout = open(options.out + '/' + 'manifest.txt', 'w')
nfile = 0
ff = 0
fout.write('sample-id,absolute-filepath,direction\n')
for root, dirs, files in os.walk(options.input):
  if len(dirs) == 0:
    root = os.path.abspath(root)
    for fl in files:
      if re.search(sp, fl):
        try:
          info = "%s/%s" % (root, fl)
          sn = id_ds[re.search(sp, fl).group(1)]
          if re.search(rp, fl):
            fout.write("%s,%s,reverse\n" % (sn, info))
            nfile += 1
          elif re.search(fp, fl):
            fout.write("%s,%s,forward\n" % (sn, info))
            nfile += 1
        except KeyError:
          ff += 1

fout.close()
emm = "\n    %s fastq file were writed" % (nfile)
emm1 = "\n    %s fastq file were filterd as the ids of samples were not found in mapping file" % (ff)
print(emm)
print(emm1)
