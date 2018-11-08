import argparse
import re
import sys
import os
# argument:
p = argparse.ArgumentParser()
p.add_argument('-u', '--unite-taxonomy', dest='unite', metavar='txt',
               help='its taxonomy')
p.add_argument('-s', '--silva-taxonomy', dest='silva', metavar='txt',
               help='silva taxonomy')
p.add_argument('-o', '--out-path', dest='out', metavar='path',
               help='out path')


options = p.parse_args()


phylums = []
classes = []

with open(options.unite, 'r') as unite:
    for line in unite:
        li = re.search('[^\t]+$', line).group().strip()
        phylums.append(re.search('p__([^;]*)', li).group(1).strip())
        classes.append(re.search('c__([^;]*)', li).group(1).strip())

phylums = list(set(phylums))
classes = list(set(classes))


def find_level(tar, level='p'):
    if level == 'p':
        for phylum in phylums:
            if re.search('__' + phylum + ' *;', tar):
                return('p__' + phylum)
        return('')

    elif level == 'c':
        for clas in classes:
            if re.search('__' + clas + ' *;', tar):
                ll = re.search('__' + clas + ' *;.*', tar).group().strip()
                for ambi in ['Incertae', 'Ambiguous', 'uncultured', 'unclassified', 'unidentified']:
                    ll = re.sub('[^;]*' + ambi + '.*',
                                'unidentified;unidentified;unidentified;unidentified;unidentified', ll)
                ll = re.split(';', ll)[0:5]
                ll = [re.sub('^.*__', '', l.strip()) for l in ll]
                tll = ''
                for cp in zip(['c__', 'o__', 'f__', 'g__', 's__'], ll):
                    tll = tll + ';' + cp[0] + cp[1]
                tll = re.sub(' ', '_', tll)
                return(tll)
        return('')


with open(options.silva, 'r') as silva, open(options.out, 'w') as out:
    for line in silva:
        line = re.split('\t', line)
        tax = line[1].strip()
        if re.search('D_3__Fungi', tax):
            tax = re.search('D_3__Fungi.*', tax).group()
            found_phlum = find_level(tax, level='p')
            if found_phlum:
                found_class = find_level(tax, level='c')
                out.write(line[0] + '\t' + 'k__Fungi;' + found_phlum + found_class + '\n')
            else:
                out.write(line[0] + '\t' + 'k__Fungi' + '\n')
        else:
            out.write(line[0] + '\t' + 'Unassigned' + '\n')
