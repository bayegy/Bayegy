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
p.add_argument('-r', '--ref-seqs', dest='ref', default=False,
               help='Supply this to clean the refrence sequences not belong to Fungi')
p.add_argument('-o', '--out-dir', dest='out', metavar='directory', default='./',
               help='out directory')


options = p.parse_args()


if not os.path.exists(options.out):
    os.makedirs(options.out)


unites = []

with open(options.unite, 'r') as unite:
    for line in unite:
        li = re.split('\t *', line)
#        li = re.sub('Incertae', 'unidentified', li[1])
        unites.append(li[1])


notspecies = ['Incertae', 'Ambiguous', 'uncultured', 'unclassified', 'unidentified']
levels = ['k__', 'p__', 'c__', 'o__', 'f__', 'g__', 's__']


def checklevel(level):
    for ncls in notspecies:
        if not level.find(ncls) == -1:
            return(False)
    return(True)


def clean_tax(tax):
    #    for ambi in ['Incertae', 'Ambiguous', 'uncultured', 'unclassified', 'unidentified']:
    for ambi in ['Ambiguous', 'uncultured', 'unclassified']:
        tax = re.sub(';[^;]*' + ambi + '[^\n]*', '', tax)
    return(tax)


def adjust_species(tax):
    if re.search('_sp\.', tax):
        spe = re.search('[^;]*_sp\.[^;\n]*', tax).group()
        lev = re.search('^[a-z]__', spe).group()
        if not lev == 's__':
            spe = re.sub('^[a-z]__', 's__', spe)
            tax = re.sub('[^;]*_sp\..*\n', '', tax)
            exn = 7 - len(re.split(';', tax))
            incert = ''
            for lvl in levels[:-1][-exn:]:
                incert = incert + lvl + 'Incertae_Sedis' + ';'
            return(tax + incert + spe + '\n')
    return(tax)


def merge_tax(tax):
    if re.search('D_3__Fungi', tax):
        hind_tax = []
        tax = re.search('D_3__Fungi.*', tax).group()
#        tax = re.sub('unidentified', '', tax)
        tax = re.sub('; *D_\d+__ *[;\n].*', '', tax)
        tax = re.sub('D_\d+__', '', tax)
        tax = re.split(';', tax)
        tax = [l.strip() for l in tax]
        for ta in reversed(tax):
            if checklevel(ta):
                ta = re.sub(' ', '_', ta)
                hind_tax = [ta] + hind_tax
                for ref in unites:
                    if re.search('__' + ta + ' *[;\n]', ref):
                        pre_tax = re.search('(^.*__' + ta + ')' + ' *[;\n]', ref).group(1)
                        del hind_tax[0]
                        if len(hind_tax) > 0:
                            hind_tax = hind_tax + ['unclassified', 'unclassified',
                                                   'unclassified', 'unclassified', 'unclassified', 'unclassified']
                            n_l = 7 - len(re.split(';', pre_tax))
                            f_hind_tax = []
                            for ite in zip(levels[-n_l:], hind_tax[:n_l]):
                                f_hind_tax.append(ite[0] + ite[1])
                            return(pre_tax + ';' + ';'.join(f_hind_tax) + '\n')
                        else:
                            return(pre_tax + '\n')
    else:
        if not options.ref:
            #   return(tax)
            return('Unassigned' + '\n')
        else:
            return(False)


s_otuid = []
with open(options.silva, 'r') as silva, open(options.out + '/' + 'transformed_silva_taxonomy.txt', 'w') as out:
    for line in silva:
        li = re.split('\t', line)
        taxonomy = merge_tax(li[1])
        if taxonomy:
            taxonomy = adjust_species(taxonomy)
            if options.clean:
                taxonomy = clean_tax(taxonomy)
            out.write(li[0].strip() + '\t' + taxonomy)
            s_otuid.append(li[0].strip())

if options.ref:
    with open('%s/filtered_ref_seqs.fasta' % (options.out), 'w') as fout, open(options.ref, 'r') as ref:
        sn = []
        for ln, line in enumerate(ref):
            line = re.sub('^>', '', line.strip())
            if ln in sn:
                fout.write(line + "\n")

            if line in s_otuid:
                fout.write('>' + line + "\n")
                sn.append(ln + 1)
