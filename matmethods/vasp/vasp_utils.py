# coding: utf-8

from __future__ import division, print_function, unicode_literals, \
    absolute_import

import six

__author__ = 'Anubhav Jain, Kiran Mathew'
__email__ = 'ajain@lbl.gov, kmathew@lbl.gov'


def get_calc_key(d, fw_spec, key):
    val = None
    calc_loc_key = "path" if key == "calc_dir" else key
    if key in d:
        val = d[key]
    elif d.get("calc_loc"):
        if isinstance(d["calc_loc"], six.string_types):
            for doc in reversed(fw_spec["calc_locs"]):
                if doc["name"] == d["calc_loc"]:
                    val = doc[calc_loc_key]
                    break
        else:
            val = fw_spec["calc_locs"][-1][calc_loc_key]
    else:
        return None
    return val
