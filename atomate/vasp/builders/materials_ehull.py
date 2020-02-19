# coding: utf-8


from tqdm import tqdm

from atomate.utils.utils import get_database

from pymatgen import MPRester, Structure
from pymatgen.entries.computed_entries import ComputedEntry

from atomate.utils.utils import get_logger
from atomate.vasp.builders.base import AbstractBuilder

logger = get_logger(__name__)

__author__ = 'Anubhav Jain <ajain@lbl.gov>'


class MaterialsEhullBuilder(AbstractBuilder):
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
        logger.info("MaterialsEhullBuilder starting...")
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

                # TODO: @computron This only calculates Ehull with respect to Materials Project.
                # It should also account for the current database's results. -computron
                self._materials.update_one({"material_id": m["material_id"]},
                                           {"$set": {"stability": self.mpr.get_stability([my_entry])[0]}})

                # TODO: @computron: also add additional properties like inverse hull energy?

                # TODO: @computron it's better to use PD tool or reaction energy calculator
                # Otherwise the compatibility schemes might have issues...one strategy might be
                # use MP only to retrieve entries but compute the PD locally -computron
                for el, elx in my_entry.composition.items():
                    entries = self.mpr.get_entries(el.symbol, compatible_only=True)
                    min_e = min(entries, key=lambda x: x.energy_per_atom).energy_per_atom
                    energy -= elx * min_e
                self._materials.update_one({"material_id": m["material_id"]},
                                           {"$set": {"thermo.formation_energy_per_atom": energy / structure.num_sites}})

                mpids = self.mpr.find_structure(structure)
                self._materials.update_one({"material_id": m["material_id"]}, {"$set": {"mpids": mpids}})

            except:
                import traceback
                logger.exception("<---")
                logger.exception("There was an error processing material_id: {}".format(m))
                logger.exception(traceback.format_exc())
                logger.exception("--->")

        logger.info("MaterialsEhullBuilder finished processing.")

    def reset(self):
        logger.info("Resetting MaterialsEhullBuilder")
        self._materials.update_many({}, {"$unset": {"stability": 1}})
        self._build_indexes()
        logger.info("Finished resetting MaterialsEhullBuilder")

    def _build_indexes(self):
        self._materials.create_index("stability.e_above_hull")

    @classmethod
    def from_file(cls, db_file, m="materials", **kwargs):
        """
        Get a MaterialsEhullBuilder using only a db file
        Args:
            db_file: (str) path to db file
            m: (str) name of "materials" collection
            **kwargs: other parameters to feed into the builder, e.g. mapi_key
        """
        db_write = get_database(db_file, admin=True)
        return cls(db_write[m], **kwargs)
