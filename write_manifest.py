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
               help=r'The regular expression of sample ID in file names. You must use \(\) to contain sample ID expression')
p.add_argument('-f', '--forward-file-pattern', dest='fp', default=r'\.1\.fq$', metavar='<regular expression>',
               help='The regular expression representing forward sequence in file names')
p.add_argument('-r', '--reverse-file-pattern', dest='rp', default=r'\.2\.fq$', metavar='<regular expression>',
               help='The regular expression representing reverse sequence in file names')
p.add_argument('-c', '--ignore-case', dest='ignore_case', action='store_true', help='Ignore case for sample names')
p.add_argument('-o', '--output', dest='out', metavar='<directory>', default='./',
               help='The path of output files')
options = p.parse_args()


if not os.path.exists(options.out):
    os.makedirs(options.out)


print("Assert no duplicated sample names: please make sure no duplicated sample names and did not use capital and small letter to distinguish samples if error happened\n")
id_ds = {}
with open(options.meta, 'r') as meta:
    for line in enumerate(meta):
        li = re.split('\t', line[1].strip())
        if line[0] > 0:
            first_element = li[0].strip().lower() if options.ignore_case else li[0].strip()
            if first_element in id_ds.keys():
                print("Duplicated sample names detected:")
                print(first_element)
            id_ds[first_element] = li[len(li) - 1].strip()


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
    if len(files) == 0:
        print("{} is a empty directory".format(root))
    # if len(dirs) == 0:
    else:
        root = os.path.abspath(root)
        for fl in files:
            if re.search(sp, fl):
                info = "%s/%s" % (root, fl)
                pre_id = re.search(sp, fl).group(1)
                pre_id = pre_id.lower() if options.ignore_case else pre_id
                try:
                    sn = id_ds[pre_id]
                    if re.search(rp, fl):
                        fout.write("%s,%s,reverse\n" % (sn, info))
                        # id_sets.append(pre_id)
                        nfile += 1
                    elif re.search(fp, fl):
                        fout.write("%s,%s,forward\n" % (sn, info))
                        id_sets.append(pre_id)
                        nfile += 1
                    else:
                        print("No direction pattern matched this file: {}".format(fl))
                except KeyError:
                    print("{} was found among fastq files, but not found in pre-map file".format(pre_id))
                    ff += 1
            else:
                print("File {} was filtered out".format(fl))

fout.close()
# print("\nThe following samples were :")

with open(options.meta, 'r') as mp, open(options.out + '/' + 'sample-metadata.tsv', 'w') as outfile:
    miss = 0
    for line in enumerate(mp):
        li = re.split('\t', line[1].strip())
        li = [e.strip() for e in li]
        fc = len(li) - 1
        if line[0] == 0:
            li[0] = "#SampleID"
            li[fc] = "Description"
            outfile.write('\t'.join(li) + '\n')
        else:
            pivot = li[0].lower() if options.ignore_case else li[0]
            if pivot in id_sets:
                li[0] = li[fc]
                outfile.write('\t'.join(li) + '\n')
            else:
                # print('    ' + li[0].lower())
                print("{} was found in pre-map file, but not found among fastq files".format(li[0]))
                miss += 1

print("\nIn summary:")
emm = "\n    %s fastq files were wrote into manifest file" % (nfile)
emm1 = "\n    %s fastq files were filterd as the ids of samples were not found in pre-map file" % (ff)
print(emm)
print(emm1)
print('\n    %s samples were removed from pre-map file, as the ids of samples were not found in forward fastq files\' names' % (miss))

# open(options.out + '/' + 'sample-metadata.tsv', 'w') as outmeta
#      li[0] = "#SampleID"
#      li[il] = "Description"
