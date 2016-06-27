from __future__ import absolute_import

from pymatgen import MPRester

__author__ = 'Anubhav Jain <ajain@lbl.gov>'

# TODO: consider using AbstractStructureFilter / Transmuters framework

class SubmissionFilter:

    NO_POTCARS = ['Po', 'At', 'Rn', 'Fr', 'Ra', 'Am', 'Cm', 'Bk', 'Cf', 'Es',
                  'Fm', 'Md', 'No', 'Lr']

    def __init__(self, is_valid=True, potcar_exists=True,
                 max_natoms=200, is_ordered=True, not_in_MP=True,
                 MAPI_KEY=None, require_bandstructure=False):

        self.is_valid = is_valid
        self.potcar_exists = potcar_exists
        self.max_natoms = max_natoms
        self.is_ordered = is_ordered
        self.not_in_MP = not_in_MP
        self.MAPI_KEY = MAPI_KEY
        self.require_bandstructure = require_bandstructure

    def filter(self, s):
        failures = []

        if self.is_valid:
            if not s.is_valid():
                failures.append("IS_VALID=False")

        if self.potcar_exists:
            elements = s.composition.elements
            if set(elements).intersection(set(self.NO_POTCARS)):
                failures.append("POTCAR_EXISTS=False")

        if self.max_natoms:
            if s.num_sites > self.max_natoms:
                failures.append("MAX_NATOMS=Exceeded")

        if self.is_ordered:
            if not s.is_ordered:
                failures.append("IS_ORDERED=False")

        if self.not_in_MP:
            mpr = MPRester(self.MAPI_KEY)
            mpids = mpr.find_structure(s)
            if mpids:
                if self.require_bandstructure:
                    for mpid in mpids:
                        try:
                            bs = mpr.get_bandstructure_by_material_id(mpid)
                            if bs:
                                failures.append("NOT_IN_MP=False ({})".format(mpid))
                        except:
                            pass
                else:
                    failures.append("NOT_IN_MP=False ({})".format(mpids[0]))

        if not failures:
            return s, None

        else:
            return None, failures