# coding:utf-8
import argparse
import re
import sys
import os

# argument:
p = argparse.ArgumentParser(description="This script is used to plot phylogentic tree of top abundant otus")
p.add_argument('-i', '--input', dest='input', metavar='<path>', default=False,
               help='Path of demuxed PE files')
p.add_argument('-p', '--pre-map', dest='map', metavar='<path>', default=False,
               help='raw sample metadata files. The sample ID must in first and last column. Two columns do not have to be same, but the IDs in first column must be respond to --sample-id-patern')
p.add_argument('-c', '--category', dest='group', metavar='<str>', default=False,
               help='column name of group in metadata. You can specify more than one category, and the categories must be seprated by commas')
p.add_argument('-o', '--outdir', dest='outdir', metavar='<directory>', default='./',
               help='specify the output directory')

options = p.parse_args()

if not os.path.exists(options.outdir):
    os.makedirs(options.outdir)

scriptpath = sys.path[0]
print(scriptpath)
for c in re.split(",", options.group.strip()):
    os.makedirs(options.outdir + '/' + c + "_Results")
    os.system("Rscript %s/clean_na_of_inputs.R -m %s --group %s  -o %s" %
              (scriptpath, options.map, c, options.outdir))
    os.system("python %s/write_manifest.py -i %s -m %s/cleaned_map.txt -o %s/%s_data" %
              (scriptpath, options.input, options.outdir, options.outdir, c))
    os.system('''cd %s/%s_Results&&sed -i -e '1s/%s/Group/g' ../%s_data/sample-metadata.tsv&&\
bash %s/16S_pipeline.V9.sh ../%s_data/sample-metadata.tsv 15000 1000 Group ~/database_16S/338-806/gg_13_8_99_338_806_classifier.qza ~/database_16S/338-806/gg_13_5_97_338_806_ref_seqs.qza ../%s_data/manifest.txt  none''' %
              (options.outdir, c, c, c, scriptpath, c, c))
