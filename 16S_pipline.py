# coding:utf-8
import argparse
import re
import sys
import os

# argument:
p = argparse.ArgumentParser(
    description="This script work well when analyze the paired-end sequences derived from MeiJi company. Sample usage: python ~/github/Bayegy/16S_pipeline.py -i rawData/eachsample/ -m pre_map.txt -c Group1,Group2,Group3 --classifier ~/database_16S/338-806/gg_13_8_99_338_806_classifier.qza --ref-seqs ~/database_16S/338-806/gg_13_5_97_338_806_ref_seqs.qza")
p.add_argument('-i', '--input', dest='input', metavar='<directory>', default=False,
               help='Directory of demuxed fastq files. sub directory is allowed when store fastq files')
p.add_argument('-m', '--map', dest='map', metavar='<path>', default=False,
               help='Sample metadata file. The sample ID must in first and last column. Two columns do not have to be same, but the IDs in first column must be respond to --sample-id-patern')
p.add_argument('-c', '--category', dest='group', metavar='<str>', default=False,
               help='column name of categories in metadata. You can specify more than one category, and the categories must be seprated by commas')
p.add_argument('--run-in-parallel', dest='parallel', action='store_true',
               help="When supplied, OTU clustering is operated one time for each category, otherwise OTU clustering is operated only one time for all categories")
p.add_argument('-s', '--sample-id-patern', dest='sp', default=r'raw\.split\.\(.+\)\.[12]\.fq$', metavar='<regular expression>',
               help='The regular expression of sample ID in file names. You must use \(\) to contain sample ID expression')
p.add_argument('-f', '--forward-file-pattern', dest='fp', default=r'\.1\.fq$', metavar='<regular expression>',
               help='The regular expression representing forward sequence in file names')
p.add_argument('-r', '--reverse-file-pattern', dest='rp', default=r'\.2\.fq$', metavar='<regular expression>',
               help='The regular expression representing reverse sequence in file names')
p.add_argument('-d', '--sample-depth', dest='depth', default='auto', metavar='<int or str>',
               help='Depth for subsampleing. If "auto" is supplied, the min OTU frequency of samples will be caculated and apply to this parameter')
p.add_argument('--classifier', dest='classifier', default='~/database_16S/338-806/gg_13_8_99_338_806_classifier.qza', metavar='<path>',
               help='Path of the classifier for alignment and assigning taxonomy')
p.add_argument('--ref-seqs', dest='ref', default='~/database_16S/338-806/gg_13_5_97_338_806_ref_seqs.qza', metavar='<path>',
               help='Path of the reference sequences for close reference alignment')
p.add_argument('-o', '--outdir', dest='outdir', metavar='<directory>', default='./',
               help='specify the output directory')

options = p.parse_args()

if not os.path.exists(options.outdir):
  os.makedirs(options.outdir)


options.sp = re.sub('\(', '\\\(', options.sp)
options.sp = re.sub('\)', '\\\)', options.sp)

scriptpath = sys.path[0]
if options.parallel:
  for c in re.split(",", options.group.strip()):
    os.makedirs(options.outdir + '/' + c + "_Results")
    os.system("Rscript %s/clean_na_of_inputs.R -m %s --group %s  -o %s" %
              (scriptpath, options.map, c, options.outdir))
    os.system("python %s/write_manifest.py -i %s -m %s/cleaned_map.txt -o %s/%s_data -f %s -r %s -s %s" %
              (scriptpath, options.input, options.outdir, options.outdir, c, options.fp, options.rp, options.sp))
    os.system('''cd %s/%s_Results&&sed -i -e '1s/%s/Group/g' ../%s_data/sample-metadata.tsv&&\
bash %s/body_16S_pipeline.sh ../%s_data/sample-metadata.tsv %s 1000 Group %s %s ../%s_data/manifest.txt  none''' %
              (options.outdir, c, c, c, scriptpath, c, options.depth, options.classifier, options.ref, c))
else:
  os.makedirs(options.outdir + '/Results')
  os.system("python %s/write_manifest.py -i %s -m %s -o %s/data -f %s -r %s -s %s" %
            (scriptpath, options.input, options.map, options.outdir, options.fp, options.rp, options.sp))
  os.system('''cd %s/Results&&\
bash %s/body_16S_pipeline.sh ../data/sample-metadata.tsv %s 1000 %s %s %s ../data/manifest.txt  none''' %
            (options.outdir, scriptpath, options.depth, options.group, options.classifier, options.ref))
