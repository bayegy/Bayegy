import argparse
import re
import sys
import os
# argument:
p = argparse.ArgumentParser()
p.add_argument('-u', '--unite-taxonomy', dest='unite', metavar='txt',
               help='unite taxonomy')
p.add_argument('-s', '--silva-taxonomy', dest='silva', metavar='txt',
               help='silva taxonomy')
p.add_argument('-c', '--clean_unclassified', dest='clean', action='store_true',
               help='Supply this to clean the unclassified levels')
p.add_argument('-o', '--out-path', dest='out', metavar='path',
               help='out path')


options = p.parse_args()


unites = []

with open(options.unite, 'r') as unite:
    for line in unite:
        li = re.split('\t *', line)
        unites.append(li[1])

levels = ['k__', 'p__', 'c__', 'o__', 'f__', 'g__', 's__']


def clean_tax(tax):
    for ambi in ['Incertae', 'Ambiguous', 'uncultured', 'unclassified', 'unidentified']:
        tax = re.sub(';[^;]*' + ambi + '[^\n]*', '', tax)
    return(tax)


def merge_tax(tax):
    if re.search('D_3__Fungi', tax):
        hind_tax = []
        tax = re.search('D_3__Fungi.*', tax).group()
        tax = re.sub('unidentified', '', tax)
        tax = re.sub('D_\d+__ *[;|\n].*', '', tax)
        tax = re.sub('D_\d+__', '', tax)
        tax = re.split(';', tax)
        tax = [l.strip() for l in tax]
        for ta in reversed(tax):
            ta = re.sub(' ', '_', ta)
            hind_tax = [ta] + hind_tax
            for ref in unites:
                if re.search('__' + ta + ' *[;|\n]', ref):
                    pre_tax = re.search('(^.*__' + ta + ')' + ' *[;|\n]', ref).group(1)
                    del hind_tax[0]
                    if len(hind_tax) > 0:
                        hind_tax = hind_tax + ['unidentified', 'unidentified',
                                               'unidentified', 'unidentified', 'unidentified', 'unidentified']
                        n_l = 7 - len(re.split(';', pre_tax))
                        f_hind_tax = []
                        for ite in zip(levels[-n_l:], hind_tax[:n_l]):
                            f_hind_tax.append(ite[0] + ite[1])
                        return(pre_tax + ';' + ';'.join(f_hind_tax) + '\n')
                    else:
                        return(pre_tax + '\n')
    else:
        return(tax)


with open(options.silva, 'r') as silva, open(options.out, 'w') as out:
    for line in silva:
        li = re.split('\t', line)
        taxonomy = merge_tax(li[1])
        if options.clean:
            taxonomy = clean_tax(taxonomy)
        out.write(li[0] + '\t' + taxonomy)
