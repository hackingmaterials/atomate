# coding: utf-8

from __future__ import division, print_function, unicode_literals, \
    absolute_import

import six

__author__ = 'Anubhav Jain'
__email__ = 'ajain@lbl.gov'

def get_calc_dir(d, fw_spec):
    calc_dir = None
    if "calc_dir" in d:
        calc_dir = d["calc_dir"]

    elif d.get("calc_loc"):
        if isinstance(d["calc_loc"], six.string_types):
            for doc in reversed(fw_spec["calc_locs"]):
                if doc["name"] == d["calc_loc_name"]:
                    calc_dir = doc["path"]
                    break
        else:
            calc_dir = fw_spec["calc_locs"][-1]["path"]

    else:
        raise ValueError("Must specify either calc_dir or calc_loc!")
    return calc_dir
