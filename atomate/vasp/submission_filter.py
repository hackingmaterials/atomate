# coding: utf-8


from monty.json import MSONable, MontyDecoder

from pymatgen import MPRester
from pymatgen.alchemy.filters import AbstractStructureFilter

__author__ = 'Anubhav Jain <ajain@lbl.gov>, Kiran Mathew <kmathew@lbl.gov>'


class SubmissionFilter(AbstractStructureFilter):
    NO_POTCARS = ['Po', 'At', 'Rn', 'Fr', 'Ra', 'Am', 'Cm', 'Bk',
                  'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr']

    def __init__(self, is_valid=True, potcar_exists=True, max_natoms=200, is_ordered=True,
                 not_in_MP=True, MAPI_KEY=None, require_bandstructure=False):
        """
        Initialize a submission filter for checking that structures are valid for calculations.

        Args:
            is_valid (bool): If true, checks structure validity
            potcar_exists (bool): If true, ensures all elements have VASP PAW_PBE POTCAR
            max_natoms (int): If not None, ensures structure has <=max_natoms atoms
            is_ordered (bool): If true, ensures structure is ordered
            not_in_MP (bool): If true, ensures structure not in MP
            MAPI_KEY (str): For MP checks, your MAPI key if not previously set as config var
            require_bandstructure (bool): For MP checks, require a band structure calc
        """
        self.is_valid = is_valid
        self.potcar_exists = potcar_exists
        self.max_natoms = max_natoms
        self.is_ordered = is_ordered
        self.not_in_MP = not_in_MP
        self.MAPI_KEY = MAPI_KEY
        self.require_bandstructure = require_bandstructure

    def test(self, structure):
        failures = []

        if self.is_valid:
            if not structure.is_valid():
                failures.append("IS_VALID=False")

        if self.potcar_exists:
            elements = structure.composition.elements
            if set(elements).intersection(set(self.NO_POTCARS)):
                failures.append("POTCAR_EXISTS=False")

        if self.max_natoms:
            if structure.num_sites > self.max_natoms:
                failures.append("MAX_NATOMS=Exceeded")

        if self.is_ordered:
            if not structure.is_ordered:
                failures.append("IS_ORDERED=False")

        if self.not_in_MP:
            mpr = MPRester(self.MAPI_KEY)
            mpids = mpr.find_structure(structure)
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
        return True if not failures else False

    def as_dict(self):
        return MSONable.as_dict(self)

    @classmethod
    def from_dict(cls, d):
        decoded = {k: MontyDecoder().process_decoded(v) for k, v in d.items() if
                   not k.startswith("@")}
        return cls(**decoded)
