#!/usr/bin/evn python3
# ######################################### import ################################################
import argparse
import os
import sys
from datetime import datetime
import pandas as pd
import numpy as np
import seaborn as sns
import json
import collections
# ########################################### ___ #################################################
__doc__ = '用于生成分组配色方案，输出文件为json格式'
__author__ = 'Liu Jiang'
__mail__ = 'liujiang9201@163.com'
__date__ = '2019年04月23日 星期二 20时23分40秒'
__version__ = '1.0.0'
# ########################################### main ##################################################


def Report(level, info):
    date_now = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    if level == "ERROR":
        sys.stderr.write("{0} - {1} - ERROR - {2}\n".format(date_now, os.path.basename(__file__), info))
        sys.exit(1)
    elif level == "INFO":
        sys.stdout.write("{0} - {1} - INFO - {2}\n".format(date_now, os.path.basename(__file__), info))
    elif level == "DEBUG":
        sys.stdout.write("{0} - {1} - DEBUG - {2}\n".format(date_now, os.path.basename(__file__), info))
        sys.exit(1)
    return()


def CheckFile(file):
    if os.path.exists(file):
        return(os.path.abspath(file))
    else:
        info = "输入文件 {0} 不存在！".format(file)
        Report("ERROR", info)


def CheckDir(dir):
    if not os.path.exists(dir):
        tmp = os.system("mkdir -p {0}".format(dir))
        if tmp == 0:
            info = "创建输出目录 {0}".format(dir)
            Report("INFO", info)
            return(dir)
        else:
            Report("ERROR", "创建输出目录失败")


def GetGroup(file, col):
    in_data = pd.read_table(file)
    # groups_name = list(in_data.columns)[col[0] - 1:col[1]]
    groups_name = col
    groups_dic = {}
    groups_list = []
    for group in groups_name:
        tmp_all = [i[0] for i in np.array(in_data[[group]]).tolist()]
        tmp_set = [i for i in set(tmp_all) if not pd.isnull(i)]
        # tmp_set.sort(key=tmp_all.index)
        tmp_set.sort(key=str.lower)
        groups_list += tmp_set
        groups_dic.setdefault(group, tmp_set)
    return(groups_dic, list(set(groups_list)), groups_name)


def main():
    bin = os.path.dirname(os.path.abspath(__file__))
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     epilog='author:\t{0}\nmail:\t{1}\ndate:\t{2}\nversion:\t{3}'.format(__author__, __mail__, __date__, __version__))
    parser.add_argument('-i', help='input file, like {}/mapping_file.txt'.format(bin),
                        dest='input', type=str, required=True)
    parser.add_argument('-c', help='categories names seprated by commas, like Group1,Group2',
                        dest='columns', type=str, required=True)
    parser.add_argument('-p', help='presupplied color palette(a file with each line contains a color). if not passed, color palette Accent will be used.',
                        dest='palette', type=str, required=False, default=False)
    parser.add_argument('-o', help='output file in json format', dest='output', type=str, required=True)
    args = parser.parse_args()
    info = "开始运行"
    Report("INFO", info)
    # col = [int(i) for i in args.columns.split(",")]
    col = [i.strip() for i in args.columns.split(",")]
    Report("INFO", "您指定 {} 列为分组信息".format("-".join([str(i) for i in col])))
    # check inout
    Report("INFO", "检查输入文件")
    in_file = CheckFile(args.input)
    # get groups
    Report("INFO", "获取分组信息")
    groups_dic, groups_list, groups_name = GetGroup(in_file, col)
    # get color
    Report("INFO", "获取配色方案")
    if args.palette:
        with open(args.palette, 'r') as f:
            palette = [l.strip() for l in f.read().split('\n')]
    color_list = sns.color_palette('Accent', len(groups_list)).as_hex(
    ) if not args.palette else palette[:len(groups_list)]
    # print(color_list)
    color_dic = dict(zip(groups_list, color_list))
    # appoint color
    Report("INFO", "分配分组配色")
    groups_color_dic = collections.OrderedDict()
    for group in groups_name:
        groups_color_dic.setdefault(group, [color_dic[i] for i in groups_dic[group]])
    # write json
    outfile = os.path.abspath(args.output)
    outdir = os.path.dirname(outfile)
    CheckDir(outdir)
    Report("INFO", "写入json文件：{}".format(outfile))
    try:
        with open(outfile, 'w') as out_file:
            json.dump(groups_color_dic, out_file, indent=4)
    except PermissionError:
        Report("ERROR", "指定目录无写入权限")


if __name__ == "__main__":
    main()
    info = "运行成功"
    Report("INFO", info)
