import argparse
import re
import sys
import os
import pandas as pd
from lxml import html
# argument:
p = argparse.ArgumentParser(
    description="This script is used to form manifest of 16S demuxed sequences for the use of qiime2")
p.add_argument('-i', '--input', dest='input', metavar='<path>',
               help='Path of input Table, either txt or html')
p.add_argument('-t', '--type', dest='type', metavar='<str>', default="auto",
               help='The type of input file')
p.add_argument('-k', '--key-word', dest='key', metavar='<str>', default="auto",
               help='The key-word supposed to be replaced by this table in target html')
p.add_argument('-o', '--output', dest='out', metavar='<directory>', default='./test.html',
               help='The path of target html to which this table will input')
options = p.parse_args()


file_type = re.search('[^\.]+$', options.input).group() if options.type == "auto" else options.type


if file_type == 'txt' or file_type == 'csv':
    out_df = pd.read_csv(options.input, sep='\t')
elif file_type == "html":
    out_df = pd.read_html(options.input)[0]
else:
    f = open(options.input, 'r')
    raw_html = f.read()
    tree = html.fromstring(raw_html)
    jason = tree.xpath("//script[@id='data']/text()")[0]
    json = pd.read_json(jason, typ='series')
    columns = [i[0] for i in json['columns']]
    out_df = pd.DataFrame(json['data'])
    out_df.columns = columns

table_html = out_df.to_html(index=False)

with open(options.out, 'r') as f1:
    target = f1.read()
    target = target.replace(options.key, '\n' + table_html + '\n')
with open(options.out, 'w') as f2:
    f2.write(target)
