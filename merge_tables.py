#!/usr/bin/env python3
import sys
import pandas as pd
import io
# from functools import reduce

sp, *srcs, tgt = sys.argv

left, *rights = [pd.read_csv(s, sep="\t", index_col=0) for s in srcs]

df_merge = left.join(rights)
if tgt == '-':
    bf = io.StringIO()
    df_merge.to_csv(bf, sep="\t")
    print(bf.getvalue())
else:
    df_merge.to_csv(tgt, sep="\t")
