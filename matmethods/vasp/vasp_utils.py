# coding: utf-8

from __future__ import division, print_function, unicode_literals, \
    absolute_import

import six

__author__ = 'Anubhav Jain'
__email__ = 'ajain@lbl.gov'

def get_vasp_dir(d, fw_spec):
    vasp_dir = None
    if "vasp_dir" in d:
        vasp_dir = d["vasp_dir"]

    elif d.get("vasp_loc"):
        if isinstance(d["vasp_loc"], six.string_types):
            for doc in reversed(fw_spec["calc_locs"]):
                if doc["name"] == d["vasp_loc_name"]:
                    vasp_dir = doc["path"]
                    break
        else:
            vasp_dir = fw_spec["calc_locs"][-1]["path"]

    else:
        raise ValueError("Must specify either vasp_dir or vasp_loc!")
    return vasp_dir
