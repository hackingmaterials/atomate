__author__ = 'Anubhav Jain <ajain@lbl.gov>'


def get_vasp_dir(d, fw_spec):

    if "vasp_dir" in d:
        vasp_dir = d["vasp_dir"]

    elif d.get("vasp_loc"):
        if isinstance(d["vasp_loc"], basestring):
            for doc in reversed(fw_spec["vasp_locs"]):
                if doc["name"] == d["vasp_loc_name"]:
                    vasp_dir = doc["path"]
                    break
        else:
            vasp_dir = fw_spec["vasp_locs"][-1]["path"]

    else:
        raise ValueError("Must specify either vasp_dir or vasp_loc!")

    return vasp_dir

