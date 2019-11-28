#!/usr/bin/env python3

import re
import os
import argparse
import shutil


p = argparse.ArgumentParser(description="Change the file name suffix of all files under a specific directory. ")

p.add_argument('directory', help="directory to change suffix", default=False, metavar='<str>')
p.add_argument('--old', '-o', help="old suffix", default='txt', metavar='<str>')
p.add_argument('--new', '-n', help="new suffix, default is xls", default='xls', metavar='<str>')
p.add_argument('--skip', '-s', help="the paths contain the spcified strings (seprated by commas) will be ignored",
               default=False, metavar='<str>')


options = p.parse_args()


def contain(path, skips):
    for s in skips:
        if not path.find(s) == -1:
            return True
    return False


for root, dirs, files in os.walk(options.directory):
    for file in files:
        path = os.path.join(root, file)
        pattern1 = re.compile(r'\.{}$'.format(options.old))
        if (not options.skip or not contain(path, options.skip.split(','))) and re.search(pattern1, path):
            new_path = re.sub(pattern1, '.{}'.format(options.new), path)
            shutil.move(path, new_path)
