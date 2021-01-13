import argparse
import re
import os
import pandas as pd

# argument:
p = argparse.ArgumentParser(
    description="Summarize abundance table at every taxonomy levels")
p.add_argument('-i', '--input', dest='input', metavar='<path>',
               help='Path of the abundance table with taxonomy description as last column')
p.add_argument('--unspecified', action="store_true",
               help='pass this flag to append Unspecified at missing levels')
p.add_argument('--output-abs', dest='outabs', metavar='<path pattern>', default=None,
               help='output absolute abundance table, {} in the path will be replaced by 1,2,3...')
p.add_argument('--output-rel', dest='outrel', metavar='<path pattern>', default=None,
               help='output relative abundance table, {} in the path will be replaced by 1,2,3...')
options = p.parse_args()

out_dirs = []

if options.outabs:
    out_dirs.append(os.path.dirname(options.outabs))
if options.outrel:
    out_dirs.append(os.path.dirname(options.outrel))

for out_dir in out_dirs:
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)


def remove_prefix(taxon):
    return re.sub('^.+__', '', taxon)


def process_taxonomy(taxa, unspecified=False):
    taxa = [re.split(' *; *', t.strip()) for t in taxa]
    taxa = [[t for t in taxon if t and not t.endswith("__")] for taxon in taxa]
    taxa_len = [len(t) for t in taxa]
    total_levels = max(taxa_len)
    extended_taxa = []
    for taxon in taxa:
        levels = len(taxon)
        if levels < total_levels:
            stuffer = "unclassified"
            if unspecified:
                stuffer = taxon[-1]
                stuffer = remove_prefix(stuffer)
                stuffer = 'Unspecified_' + stuffer
            taxon = taxon + ([stuffer] * (total_levels - levels))
        extended_taxa.append(taxon)
    return total_levels, extended_taxa


def main():
    skipfirst = False
    with open(options.input) as fh:
        line = fh.readline()
        if len(line.split('\t')) < 2:
            skipfirst = True
    df = pd.read_csv(options.input, skiprows=1 if skipfirst else None, index_col=0, sep="\t")
    # print(df)
    df = df.loc[df.iloc[:, -1].notna(), :]
    taxonomy = df.iloc[:, -1]
    df.drop(df.columns[-1], axis=1, inplace=True)
    total_levels, extended_taxa = process_taxonomy(taxonomy, unspecified=options.unspecified)
    # print(total_levels)
    for ln in range(1, total_levels + 1):
        current_levels = [t[:ln] for t in extended_taxa]
        level_dict = {"unclassified": ["unclassified"]}
        merged_levels = []
        for lv in current_levels:
            m_lv = ";".join(lv)
            if m_lv.endswith("unclassified"):
                merged_levels.append("unclassified")
            else:
                merged_levels.append(m_lv)
                lv = [(tn if tn.startswith('Unspecified_') else remove_prefix(tn)) for tn in lv]
                level_dict[m_lv] = lv
        level_df = df.copy()
        level_df['detail'] = merged_levels
        level_df = level_df.groupby('detail').sum()
        layers = []
        for t in level_df.index:
            data = level_dict[t]
            if len(data) > 2:
                data = data[-2:]
            if data[-1] in layers:
                layers.append('_'.join(data))
            else:
                layers.append(data[-1])
        if options.outabs:
            abs_df = level_df.copy()
            abs_df['Tax_detail'] = abs_df.index
            abs_df.insert(0, 'Taxonomy', layers)
            abs_df.to_csv(options.outabs.format(ln), index=False, sep="\t")

        if options.outrel:
            rel_df = level_df.apply(lambda x: x / x.sum(), axis=0)
            rel_df['Tax_detail'] = rel_df.index
            rel_df.insert(0, 'Taxonomy', layers)
            rel_df.to_csv(options.outrel.format(ln), index=False, sep="\t")


if __name__ == '__main__':
    main()
