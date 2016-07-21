# coding: utf-8

from __future__ import absolute_import, division, print_function, unicode_literals

from tqdm import tqdm

from matgendb.util import get_database
from pymatgen import MPRester, Structure
from pymatgen.entries.computed_entries import ComputedEntry

__author__ = 'Anubhav Jain <ajain@lbl.gov>'


class MaterialsEhullBuilder:
    def __init__(self, materials_write, mapi_key=None, update_all=False):
        """
        Starting with an existing materials collection, adds stability information and
        The Materials Project ID.
        Args:
            materials_write: mongodb collection for materials (write access needed)
            mapi_key: (str) Materials API key (if MAPI_KEY env. var. not set)
            update_all: (bool) - if true, updates all docs. If false, only updates
                docs w/o a stability key
        """
        self._materials = materials_write
        self.mpr = MPRester(api_key=mapi_key)
        self.update_all = update_all

    def run(self):
        print("MaterialsEhullBuilder starting...")
        self._build_indexes()

        q = {"thermo.energy": {"$exists": True}}
        if not self.update_all:
            q["stability"] = {"$exists": False}

        mats = [m for m in self._materials.find(q, {"calc_settings": 1, "structure": 1,
                                                    "thermo.energy": 1, "material_id": 1})]
        pbar = tqdm(mats)
        for m in pbar:
            pbar.set_description("Processing materials_id: {}".format(m['material_id']))
            try:
                params = {}
                for x in ["is_hubbard", "hubbards", "potcar_spec"]:
                    params[x] = m["calc_settings"][x]

                structure = Structure.from_dict(m["structure"])
                energy = m["thermo"]["energy"]
                my_entry = ComputedEntry(structure.composition, energy, parameters=params)
                self._materials.update_one({"material_id": m["material_id"]},
                                           {"$set": {"stability": self.mpr.get_stability([my_entry])[0]}})

                mpids = self.mpr.find_structure(structure)
                self._materials.update_one({"material_id": m["material_id"]}, {"$set": {"mpids": mpids}})

            except:
                import traceback
                print("<---")
                print("There was an error processing material_id: {}".format(m))
                traceback.print_exc()
                print("--->")

        print("MaterialsEhullBuilder finished processing.")


    def reset(self):
        self._materials.update_many({}, {"$unset": {"stability": 1}})
        self._build_indexes()

    def _build_indexes(self):
        self._materials.create_index("stability.e_above_hull")


    @staticmethod
    def from_db_file(db_file, m="materials", **kwargs):
        """
        Get a MaterialsEhullBuilder using only a db file
        Args:
            db_file: (str) path to db file
            m: (str) name of "materials" collection
            **kwargs: other parameters to feed into the builder, e.g. mapi_key
        """
        db_write = get_database(db_file, admin=True)
        return MaterialsEhullBuilder(db_write[m], **kwargs)
