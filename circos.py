
import json
import numpy as np
import pandas as pd
import re
import os
import shutil


class Circos(object):
    def __init__(self, table, number=10, mapping_file=None, category=None, by_group_mean=False, outpath="./", prefix=""):
        self.__base_path = re.sub('[^/]+$', "", __file__)
        self.__read_conf__()
        self.__numder = number
        self.__prefix = prefix
        self.__read_colors__()
        self.__outpath = outpath + "/"
        self.__outconf = self.__outpath + "circos_conf/"
        self.__init_data__(table, mapping_file, category, by_group_mean)
        self.__number_of_otu, self.__number_of_sample = self.data.shape
        if not os.path.exists(self.__outpath):
            os.makedirs(self.__outpath)
        shutil.copytree(self.__base_path + "/circos_config", self.__outpath + 'circos_conf')
        otu_col_index, sam_col_index = self.generate_span(self.data.shape)
        self.__otu_col = self.__colors[otu_col_index[0]:otu_col_index[1]]
        self.__sam_col = self.__colors[sam_col_index[0]:sam_col_index[1]]
        self.__rev_sam_name, self.__rev_otu_name = (self.data.columns.values.reshape((-1, 1))[::-1, :],
                                                    self.data.index.values.reshape((-1, 1)))[::-1, :]
        self.__rev_sam_col = self.__sam_col[::-1, :]
        self.__rev_otu_col = self.__otu_col[::-1, :]

    def __read_conf__(self):
        with open(self.__base_path + "/circos_config/path.conf") as f:
            js = json.load(f)
            self._circos = js['circos_path']
            self._etc = js['circos_etc']

    def __init_data__(self, table, mapping_file, category, by_group_mean):
        otu = self.read_tsv(table)
        otu = otu if not otu.dtypes[-1] == np.dtype("O") else otu.drop(otu.columns[-1], axis=1)
        if mapping_file and category:
            mapf = self.read_tsv(mapping_file)[category].dropna().sort_values(by=category)
            otu = otu.filter(items=mapf.index, axis=1)
        else:
            otu = otu.sort_index(axis=1)
        otu = otu.iloc[otu.sum(axis=1).values.argsort()[::-1][:self.__numder]]
        otu = otu.groupby(mapf, axis=1).mean() if by_group_mean else otu
        self.data = otu.applymap(lambda x: int(np.round(x)))

    def __read_colors__(self):
        colors = []
        with open(self._etc + "/colors.brewer.conf") as color_file:
            for line in color_file:
                line = line.strip()
                if line and not(line.startswith("#")):
                    colors.append(re.sub('=.+', "", line).strip())
        self.__colors = np.array(colors).reshape((-1, 1))
        # print(self.__colors)

    def generate_span(self, number_list: []) -> []:
        """Please sapply iterable number list"""
        step_sum = []
        current_sum = 0
        for e in number_list:
            current_sum += e
            step_sum.append(current_sum)
        flat = list(np.array(step_sum) - np.array(number_list))
        return [[m, n] for m, n in zip(flat, step_sum)]

    def rep_each(self, x, each) -> np.array:
        return np.array([[i] * each for i in x]).reshape((-1, 1))

    def read_tsv(self, path):
        return pd.read_csv(path, sep="\t", index_col=0)

    def write_conf(self, df, name):
        df.to_csv(self.__outconf + name, sep=" ", header=False, index=False)

    def write_karyotype(self):
        all_sum = pd.concat([self.data.sum(axis=0), self.data.sum(axis=1)], axis=0)
        kar_data = pd.DataFrame(np.zeros((all_sum.shape[0], 7)))
        kar_data[0] = "chr"
        kar_data[1] = "-"
        kar_data[2] = all_sum.index
        kar_data[3] = all_sum.index
        kar_data[4] = 0
        kar_data[5] = all_sum.values
        kar_data[6] = np.vstack((self.__sam_col, self.__otu_col))
        self.write_conf(kar_data, "karyotype.txt")

    def write_highlight(self):

        data = self.data
        data = data.loc[data.index[::-1], data.columns[::-1]]
        pre = "fill_color="
        self.otu = np.array(data.apply(self.generate_span, axis=1).values.tolist())
        otu = pd.DataFrame(np.hstack((self.__rev_otu_name.repeat(self.__number_of_sample, axis=0),
                                      self.otu.reshape((-1, 2)), np.char.add(pre, np.tile(self.__rev_sam_col, (self.__number_of_otu, 1))))))
        self.write_conf(otu, "highlight_spec.txt")

        self.sample = np.array(data.apply(self.generate_span, axis=0).values.T.tolist())
        sample = pd.DataFrame(np.hstack((self.__rev_sam_name.repeat(self.__number_of_otu, axis=0),
                                         self.sample.reshape((-1, 2)), np.char.add(pre, np.tile(self.__rev_otu_col, (self.__number_of_sample, 1))))))
        self.write_conf(sample, "highlight_site.txt")

        self.write_conf(pd.DataFrame(np.vstack((sample.values, otu.values))), "highlight_all.txt")

    def write_links(self):

        link = pd.DataFrame(
            np.hstack((np.tile(self.__rev_sam_name, (self.__number_of_otu, 1)),
                       self.sample.reshape((-1, 2), order="F"),
                       self.__rev_otu_name.repeat(self.__number_of_sample, axis=0),
                       self.otu.reshape((-1, 2)),
                       np.char.add("color=", self.__rev_otu_col.repeat(self.__number_of_sample, axis=0))
                       ))
        )
        self.write_conf(link, "links.txt")

    def __init__path__(self):
        with open(self.__outconf + "image.generic.conf", 'r', encoding='utf-8') as ci:
            out = ci.read() % (self.__outpath, self.__prefix + 'circos.png')

        with open(self.__outconf + "image.generic.conf", 'w', encoding='utf-8') as co:
            co.write(out)

    def plot_circos(self):
        self.write_karyotype()
        self.write_highlight()
        self.write_links()
        self.__init__path__()
        os.system(self._circos + " -conf " + self.__outconf + "circos.conf")


if __name__ == '__main__':
    import argparse
    p = argparse.ArgumentParser(
        description="This script is used to plot RDA of species. The numeric enviroment factors must be encluded in maping file. The categories will be filterd before RDA")
    p.add_argument('-i', '--input', dest='input', metavar='<path>',
                   help='Taxonomic count data file')
    p.add_argument('-o', '--output', dest='output', metavar='<directory>', default='./',
                   help='Given an output directory')
    p.add_argument('-m', '--map', dest='map', metavar='<path>',
                   help='Sample metadata file')
    p.add_argument('-g', '--group', dest='group', metavar='<str>',
                   help='Column name in sample-metadata file')
    p.add_argument('-n', '--number', dest='number', metavar='<int>', default='10',
                   help='Specify how many species to be display, defaulf is 20')
    p.add_argument('-b', '--by-groupMean', dest='by', metavar='<bool>', default=False,
                   help='Pass True to use group mean to plot circos')
    p.add_argument('-p', '--prefix', dest='prefix', metavar='<int>', default="",
                   help='The prefix of output files, default if null')
    options = p.parse_args()
    c = Circos(table=options.input, number=int(options.number), mapping_file=options.map,
               category=options.group, by_group_mean=bool(options.by), outpath=options.output, prefix=options.prefix)
    c.plot_circos()
