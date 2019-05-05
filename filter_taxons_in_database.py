import argparse
import re
import sys
import os
# argument:
p = argparse.ArgumentParser()
p.add_argument('-r', '--ref-seqs', dest='ref', default=False,
               help='reference sequences, (.fasta)')
p.add_argument('-t', '--taxonomy', dest='taxonomy', metavar='txt',
               help='taxonomy, (tabular table)')
p.add_argument('-e', '--exclude', dest='ex', metavar='str',
               help='One or more search terms that indicate which taxa should be excluded from the database. If providing more than one term,terms should be delimited by commas.')
p.add_argument('-c', '--clean_unclassified', dest='clean', action='store_true',
               help='Supply this to clean the unclassified levels')
p.add_argument('-o', '--out-dir', dest='out', metavar='directory', default='./',
               help='out directory')


options = p.parse_args()


if not os.path.exists(options.out):
    os.makedirs(options.out)


def clean_tax(tax):
    #    for ambi in ['Incertae', 'Ambiguous', 'uncultured', 'unclassified', 'unidentified']:
    for ambi in ['Ambiguous', 'uncultured', 'unclassified']:
        tax = re.sub(';[^;]*' + ambi + '[^\n]*', '', tax)
    return(tax)


def is_clean(tax):
    ex = re.split(',', options.ex.strip(','))
    for e in ex:
        if not tax.find(e) == -1:
            return False
    return True


s_otuid = []
with open(options.taxonomy, 'r') as taxonomy, open(options.out + '/' + 'filtered_taxonomy.txt', 'w') as out:
    for line in taxonomy:
        li = re.split('\t', line)
        tax = li[1]
        if is_clean(tax):
            if options.clean:
                tax = clean_tax(tax)
            out.write(li[0].strip() + '\t' + tax)
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
