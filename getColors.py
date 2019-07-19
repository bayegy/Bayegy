import os
import json
import numpy as np
import pandas as pd
# import pdb


def get_colors(category):
    try:
        colors_plan_path = os.environ['COLORS_PLAN_PATH']
        with open(colors_plan_path) as f:
            js = json.load(f)
            return js[category]
    except KeyError:
        return []


def get_lefse_colors(category, mapping_file, lda_file, return_dict=False):
    colors = get_colors(category)
    if colors:
        df = pd.read_csv(mapping_file, sep="\t")
        gps = df[category]
        gps = list(set(gps[gps.notna()]))
        gps.sort(key=str.lower)
        gps_colors = dict()
        gps_sig = set()
        for gp, color in zip(gps, colors):
            gps_colors[gp] = color
        with open(lda_file) as f:
            for line in f:
                li = line.strip().split('\t')
                sgp = li[2].strip()
                if sgp:
                    gps_sig.add(sgp)
        gps_sig = list(gps_sig)
        gps_sig.sort()
        return {gp: gps_colors[gp] for gp in gps_sig} if return_dict else [gps_colors[gp] for gp in gps_sig]
    else:
        return {} if return_dict else []
