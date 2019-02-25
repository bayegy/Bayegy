import sys
import re
import os
import numpy as np
import pandas as pd
import pymysql
import time
import getopt
import warnings
warnings.filterwarnings("ignore")

db = pymysql.connect("localhost", 'root', '947366', 'wstdb_project')
cursor = db.cursor()


def add(record):
    sql = "insert into tb_project (项目编号,客户姓名,项目类型,分析时间,样本数量) values ('%s','%s','%s','%s',%s)" % (
        record[0], record[1], record[2], record[3], record[4])
    print(sql)
    cursor.execute(sql)


def export(tm):
    sql = "select * from tb_project where 分析时间 like '%s%%'" % (tm)
    print(sql)
    cursor.execute(sql)
    dt = cursor.fetchall()
    dt = np.array(dt)
    df = pd.DataFrame(dt[:, 1:])
    df.columns = ("项目编号", "客户姓名", "项目类型", "分析时间", "样本数量")
    df = df.append({'项目编号': '合计', '样本数量': df['样本数量'].sum()}, ignore_index=True)
    print(df)
    return df


def parse_project(line):
    line = line.strip()
    name = re.search('^[^\-]+', line).group().replace('老师', '')
    if re.search('ME[0-9]+', line):
        number = re.search('ME[0-9]+', line).group()
    else:
        number = time.strftime('RD%Y%m%d', time.localtime(time.time()))
    if re.search('([0-9]+)(个|份)', line):
        sample = re.search('([0-9]+)(个|份)', line).group(1)
    else:
        sample = '0'
    atime = re.search('[^\-]+$', line).group().replace('.rar', '')
    ptype = re.search('\-([^\-]+)\-分析结果', line).group(1)
    return [number, name, ptype, atime, sample]


try:
    if sys.argv[1] == 'add':
        if re.search('\.txt$', sys.argv[2]):
            with open(sys.argv[2], 'r') as fin:
                for line in fin:
                    add(parse_project(line))
        else:
            add(parse_project(sys.argv[2]))
        db.commit()
    elif sys.argv[1] == 'export':
        export(sys.argv[2]).to_csv(sys.argv[2] + '分析项目统计.xls', index=False, sep='\t', na_rep="--")

except Exception as e:
    print(e)
    db.rollback()
finally:
    db.close()
