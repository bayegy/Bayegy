import os
import json


def get_colors(category):
    try:
        colors_plan_path = os.environ['COLORS_PLAN_PATH']
        with open(colors_plan_path) as f:
            js = json.load(f)
            return js[category]
    except KeyError:
        return False
