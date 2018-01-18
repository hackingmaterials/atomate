from tqdm import tqdm

from atomate.utils.utils import get_database

from pymatgen import Structure
from pymatgen.analysis.structure_analyzer import get_dimensionality

from atomate.utils.utils import get_logger
from atomate.vasp.builders.base import AbstractBuilder

logger = get_logger(__name__)

__author__ = 'Anubhav Jain <ajain@lbl.gov>'


class MaterialsDescriptorBuilder(AbstractBuilder):
    def __init__(self, materials_write, update_all=False):
        """
        Starting with an existing materials collection, adds some compositional and structural
        descriptors.
        
        Args:
            materials_write: mongodb collection for materials (write access needed)
            update_all: (bool) - if true, updates all docs. If false, updates incrementally
        """
        self._materials = materials_write
        self.update_all = update_all

    def run(self):
        logger.info("MaterialsDescriptorBuilder starting...")
        self._build_indexes()

        q = {}
        if not self.update_all:
            q["descriptors.density"] = {"$exists": False}

        mats = [m for m in self._materials.find(q, {"structure": 1, "material_id": 1})]

        pbar = tqdm(mats)
        for m in pbar:
            pbar.set_description("Processing materials_id: {}".format(m['material_id']))
            struct = Structure.from_dict(m["structure"])
            d = {"descriptors": {}}
            d["descriptors"]["dimensionality"] = get_dimensionality(struct)
            d["descriptors"]["density"] = struct.density
            d["descriptors"]["nsites"] = len(struct)
            d["descriptors"]["volume"] = struct.volume

            self._materials.update_one({"material_id": m["material_id"]}, {"$set": d})

    def reset(self):
        logger.info("Resetting MaterialsDescriptorBuilder")
        self._materials.update_many({}, {"$unset": {"descriptors": 1}})
        self._build_indexes()
        logger.info("Finished resetting MaterialsDescriptorBuilder")

    def _build_indexes(self):
        for i in ["descriptors.dimensionality", "descriptors.density", "descriptors.nsites"]:
            self._materials.create_index(i)

    @classmethod
    def from_file(cls, db_file, m="materials", **kwargs):
        """
        Get a MaterialsDescriptorBuilder using only a db file.

        Args:
            db_file: (str) path to db file
            m: (str) name of "materials" collection
            **kwargs: other parameters to feed into the builder
        """
        db_write = get_database(db_file, admin=True)
        return cls(db_write[m], **kwargs)