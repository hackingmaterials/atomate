# coding: utf-8



from atomate.feff.workflows.core import get_wf_xas

__author__ = 'Kiran Mathew'
__email__ = 'kmathew@lbl.gov'


def wf_Xanes_K_edge(structure, c=None):

    c = c or {}
    feff_cmd = c.get("FEFF_CMD", "feff_mpi")
    db_file = c.get("DB_FILE", None)
    metadata = c.get("METADATA", {})
    absorbing_atom = c.get("ABSORBING_ATOM")

    user_tag_settings = {"RPATH": -1,
                         "SCF": "7 0 30 0.2 3",
                         "FMS": "9 0",
                         "LDOS": "-30.0 30.0 0.1",
                         "RECIPROCAL": "",
                         "EDGE": "K",
                         "COREHOLE": "RPA"}

    metadata.update({"absorbing_atom_idx": str(absorbing_atom)})

    wf = get_wf_xas(absorbing_atom, structure, edge="K", feff_cmd=feff_cmd, db_file=db_file,
                    metadata=metadata, user_tag_settings=user_tag_settings, use_primitive=False)
    return wf
