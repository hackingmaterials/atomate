from pymatgen import MPRester, Structure
from pymatgen.entries.computed_entries import ComputedEntry

__author__ = 'Anubhav Jain <ajain@lbl.gov>'


class MaterialsEhullBuilder:
    def __init__(self, materials_write, mapi_key=None, update_all=True):
        self._materials = materials_write
        self.mpr = MPRester(mapi_key)
        self.update_all = update_all

    def run(self):
        q = {"thermo.energy": {"$exists": True}}
        if not self.update_all:
            q["stability"] = {"$exists": False}

        for m in self._materials.find(q, {"calc_settings": 1, "structure": 1, "thermo.energy": 1, "material_id": 1}):
            print("Processing material_id: {}".format(m["material_id"]))
            params = {}
            for x in ["is_hubbard", "hubbards", "potcar_spec"]:
                params[x] = m["calc_settings"][x]

            composition = Structure.from_dict(m["structure"]).composition
            energy = m["thermo"]["energy"]

            my_entry = ComputedEntry(composition, energy, parameters=params)
            self._materials.update_one({"material_id": m["material_id"]}, {"$set": {"stability": self.mpr.get_stability([my_entry])[0]}})

        print("MaterialsEhullBuilder finished processing.")