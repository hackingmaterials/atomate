from tqdm import tqdm

from pymatgen import MPRester, Structure
from pymatgen.entries.computed_entries import ComputedEntry

__author__ = 'Anubhav Jain <ajain@lbl.gov>'


class MaterialsEhullBuilder:
    def __init__(self, materials_write, mapi_key=None, update_all=True):
        """
        Starting with an existing materials collection, adds stability information.
        Args:
            materials_write: mongodb collection for materials (write access needed)
            mapi_key: (str) Materials API key (if MAPI_KEY env. var. not set)
            update_all: (bool) - if true, updates all docs. If false, only updates docs w/o a stability key
        """
        self._materials = materials_write
        self.mpr = MPRester(mapi_key)
        self.update_all = update_all

    def run(self):
        print("MaterialsEhullBuilder starting...")
        q = {"thermo.energy": {"$exists": True}}
        if not self.update_all:
            q["stability"] = {"$exists": False}

        mats = [m for m in self._materials.find(q, {"calc_settings": 1, "structure": 1, "thermo.energy": 1, "material_id": 1})]
        pbar = tqdm(mats)
        for m in pbar:
            pbar.set_description("Processing materials_id: {}".format(m['material_id']))
            try:
                params = {}
                for x in ["is_hubbard", "hubbards", "potcar_spec"]:
                    params[x] = m["calc_settings"][x]

                composition = Structure.from_dict(m["structure"]).composition
                energy = m["thermo"]["energy"]

                my_entry = ComputedEntry(composition, energy, parameters=params)
                self._materials.update_one({"material_id": m["material_id"]}, {"$set": {"stability": self.mpr.get_stability([my_entry])[0]}})

            except:
                import traceback
                print("<---")
                print("There was an error processing material_id: {}".format(m))
                traceback.print_exc()
                print("--->")

        print("MaterialsEhullBuilder finished processing.")
