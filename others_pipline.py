# coding:utf-8
import argparse
import re
import sys
import os

# argument:
p = argparse.ArgumentParser(
    description="This script work well when analyze the paired-end sequences derived from MeiJi company. python ~/github/Bayegy/16S_pipline.py -i 郑超君老师-ME201808031003-MJ-M-20180807049-68份样本-原始数据/rawData/eachsample/ --classifier ../../database_16S/GG/338-806/gg_13_8_99_338_806_classifier.qza --classifier-type gg -m pre_map.txt -c Group1,Group2,Group3 -o Results --run-in-parallel")
p.add_argument('-i', '--input', dest='input', metavar='<directory>', default=False,
               help='Directory of demuxed fastq files. sub directory is allowed when storing fastq files. This parameter is required if --manifest is not supplied. Supply this parameter together with -s -f -r. This option will be ignored if --manifest is specifed')
p.add_argument('--manifest', dest='manifest', metavar='<path>', default=False,
               help='The manifest file to import demuxed fastq files, required if -i is not specified. And -i, -s, -f, -r will be ignored if this parameter is specifed')
p.add_argument('-m', '--map', dest='map', metavar='<path>', default=False,
               help='Sample metadata file. The sample ID must in first and last column, and two columns do not have to be same when -i is supplied, but the IDs in first column must be respond to --sample-id-patern')
p.add_argument('-c', '--category', dest='group', metavar='<str>', default='auto',
               help='column name of categories in metadata. You can specify more than one category, and the categories must be seprated by commas')
p.add_argument('--run-in-parallel', dest='parallel', action='store_true',
               help="When supplied, OTU clustering is operated for each category, otherwise OTU clustering is operated only one time for all categories. This parameter is only used when --manifest is not supplied")
p.add_argument('-s', '--sample-id-patern', dest='sp', default=r'raw\.split\.(.+)\.[12]\.fq$', metavar='<regular expression>',
               help="The regular expression of sample ID in file names. You must use () to contain sample ID expression, and to prevent unexpected error, using '' to contain the expression is necessary. Supply this parameter together with -i -f -r. This option will be ignored if --manifest is specifed")
p.add_argument('-f', '--forward-file-pattern', dest='fp', default=r'\.1\.fq$', metavar='<regular expression>',
               help='The regular expression representing forward sequence in file names,Supply this parameter together with -i -s -r. This option will be ignored if --manifest is specifed')
p.add_argument('-r', '--reverse-file-pattern', dest='rp', default=r'\.2\.fq$', metavar='<regular expression>',
               help='The regular expression representing reverse sequence in file names,Supply this parameter together with -i -s -f. This option will be ignored if --manifest is specifed')
p.add_argument('-d', '--sample-depth', dest='depth', default='auto', metavar='<int or str>',
               help='Depth for subsampleing. If "auto" is supplied, the min OTU frequency of samples will be caculated and apply to this parameter')
p.add_argument('-e', '--exclude', dest='exclude', default='none', metavar='<str>',
               help='The variables excluded from RDA and correlation heatmap analysis, if more than one variable is expected, you should use comma as delimeters. if more than one run is expected, use semicolon as delimeters. example: "Var1,Var2;Var1,Var3". default is none, meaning no variable will be excluded.')
p.add_argument('--classifier', dest='classifier', default='/home/admin1/database_16S/GG/338-806/gg_13_8_99_338_806_classifier.qza', metavar='<path>',
               help='Path of the classifier for alignment and assigning taxonomy')
p.add_argument('--ref-seqs', dest='ref', default='/home/admin1/database_16S/GG/338-806/gg_13_5_97_338_806_ref_seqs.qza', metavar='<path>',
               help='Path of the reference sequences for close reference alignment')
p.add_argument('--classifier-type', dest='type', default='gg', metavar='<str>',
               help='Specify the type of classifier, either silva or gg')
p.add_argument('-n', '--rar-filename', dest='fname', default='Result_AmpliconSequencing', metavar='<str>',
               help='The name of final results.rar')
p.add_argument('-o', '--outdir', dest='outdir', metavar='<directory>', default='./',
               help='specify the output directory')

options = p.parse_args()


if not ((options.input or options.manifest) and options.map and options.group):
    print("You must specify either -i or --manifest, also -m and -c")
    sys.exit()

if not os.path.exists(options.outdir):
    os.makedirs(options.outdir)
if options.manifest:
    options.manifest = os.path.abspath(options.manifest)
options.map = os.path.abspath(options.map)
options.classifier = os.path.abspath(options.classifier)
options.ref = os.path.abspath(options.ref)
options.outdir = os.path.abspath(options.outdir)
scriptpath = sys.path[0]


def parse_map(mapping):
    os.system('''sed -i -e 's/ //g' %s''' % (mapping))
    with open(mapping, 'r') as mapin:
        l1 = mapin.readline()
        gps = re.findall('Group\d*', l1)
    return ','.join(gps)


def rar(name, target):
    os.system('''rar a %s %s''' % (name, target))


if os.path.isfile(options.map):
    cgp = parse_map(options.map)
    options.group = cgp if options.group == "auto" else options.group
    if not options.manifest:
        if options.parallel:
            for c in re.split(",", options.group.strip()):
                print('##############################################################Start to generate results according to the samples of %s' % (c))
                if not os.path.exists(options.outdir + '/' + c + "_Results"):
                    os.makedirs(options.outdir + '/' + c + "_Results")
                os.system("Rscript %s/clean_na_of_inputs.R -m %s --group %s  -o %s" %
                          (scriptpath, options.map, c, options.outdir))
                os.system("python %s/write_manifest.py -i %s -m %s/cleaned_map.txt -o %s/%s_data -f '%s' -r '%s' -s '%s'" %
                          (scriptpath, options.input, options.outdir, options.outdir, c, options.fp, options.rp, options.sp))
                os.system('''cd %s/%s_Results&&sed -i -e '1s/%s/Group/g' ../%s_data/sample-metadata.tsv&&\
      bash %s/others_pipeline.sh ../%s_data/sample-metadata.tsv %s 1000 Group %s %s ../%s_data/manifest.txt  "%s" %s''' %
                          (options.outdir, c, c, c, scriptpath, c, options.depth, options.classifier, options.ref, c, options.exclude, options.type))
                if not os.path.exists(options.outdir + '/' + 'All_Results_Summary'):
                    os.makedirs(options.outdir + '/' + 'All_Results_Summary')
                os.system('''rm -r %s/All_Results_Summary/%s_Result_AmpliconSequencing; mv %s/%s_Results/Result_AmpliconSequencing %s/All_Results_Summary/%s_Result_AmpliconSequencing''' %
                          (options.outdir, c, options.outdir, c, options.outdir, c))

        else:
            if not os.path.exists(options.outdir + '/Results'):
                os.makedirs(options.outdir + '/Results')
            os.system("python %s/write_manifest.py -i %s -m %s -o %s/data -f '%s' -r '%s' -s '%s'" %
                      (scriptpath, options.input, options.map, options.outdir, options.fp, options.rp, options.sp))
            os.system('''cd %s/Results&&\
      bash %s/others_pipeline.sh ../data/sample-metadata.tsv %s 1000 %s %s %s ../data/manifest.txt  "%s" %s&&\
      rar a %s.rar Result_AmpliconSequencing''' % (options.outdir, scriptpath, options.depth, options.group, options.classifier, options.ref, options.exclude, options.type, options.fname))
    else:
        os.system('''cd %s&&\
      bash %s/others_pipeline.sh %s %s 1000 %s %s %s %s  "%s" %s&&\
      rar a %s.rar Result_AmpliconSequencing''' % (options.outdir, scriptpath, options.map, options.depth, options.group, options.classifier, options.ref, options.manifest, options.exclude, options.type, options.fname))
else:
    mapdir = os.listdir(options.map)
    excs = re.split('::', options.exclude)
    gp_set = re.split('::', options.group)
    if not os.path.exists(options.outdir + '/Results_Result_AmpliconSequencing'):
        os.makedirs(options.outdir + '/Results_Result_AmpliconSequencing')

    for mp in mapdir:
        pmp = options.map + '/' + mp
        if os.isfile(pmp):
            cgp = parse_map(pmp)
            if cgp:
                od = re.search('\d+', mp).group()
                odn = int(od) - 1
                exc = excs[odn] if len(excs) > 1 else options.exclude
                cgp = cgp if options.group == "auto" else gp_set[odn]
                if not os.path.exists(options.outdir + '/Results' + od):
                    os.makedirs(options.outdir + '/Results' + od)
                os.system("python %s/write_manifest.py -i %s -m %s -o %s/data%s -f '%s' -r '%s' -s '%s'" %
                          (scriptpath, options.input, pmp, options.outdir, od, options.fp, options.rp, options.sp))
                os.system('''cd %s/Results%s&&\
          bash %s/others_pipeline.sh ../data%s/sample-metadata.tsv %s 1000 %s %s %s ../data%s/manifest.txt  "%s" %s&&\
          mv Result_AmpliconSequencing %sResult_AmpliconSequencing&&\
          rar a %s%s.rar %sResult_AmpliconSequencing&&\
          mv %s%s.rar ../Results_Result_AmpliconSequencing/''' %
                          (options.outdir, od, scriptpath, od, options.depth, cgp, options.classifier, options.ref, od, exc, options.type,
                           od, od, options.fname, od, od, options.fname))
