# coding: utf-8


from tqdm import tqdm

from atomate.utils.utils import get_logger , get_database

from pymatgen import Structure
from pymatgen.analysis.structure_matcher import StructureMatcher, ElementComparator
from pymatgen.electronic_structure.boltztrap import BoltztrapAnalyzer

from atomate.vasp.builders.base import AbstractBuilder

logger = get_logger(__name__)

__author__ = 'Anubhav Jain <ajain@lbl.gov>'


class BoltztrapMaterialsBuilder(AbstractBuilder):
    def __init__(self, materials_write, boltztrap_read):
        """
        Update materials collection based on boltztrap collection.

        Args:
            materials_write (pymongo.collection): mongodb collection for materials (write access needed)
            boltztrap_read (pymongo.collection): mongodb collection for boltztrap (suggest read-only for safety)
        """
        self._materials = materials_write
        self._boltztrap = boltztrap_read

    def run(self):
        logger.info("BoltztrapMaterialsBuilder starting...")
        logger.info("Initializing list of all new boltztrap ids to process ...")
        previous_oids = []
        for m in self._materials.find({}, {"_boltztrapbuilder.all_object_ids": 1}):
            if '_boltztrapbuilder' in m:
                previous_oids.extend(m["_boltztrapbuilder"]["all_object_ids"])

        if not previous_oids:
            self._build_indexes()

        all_btrap_ids = [i["_id"] for i in self._boltztrap.find({}, {"_id": 1})]
        new_btrap_ids = [o_id for o_id in all_btrap_ids if o_id not in previous_oids]

        logger.info("There are {} new boltztrap ids to process.".format(len(new_btrap_ids)))

        pbar = tqdm(new_btrap_ids)
        for o_id in pbar:
            pbar.set_description("Processing object_id: {}".format(o_id))
            try:
                doc = self._boltztrap.find_one({"_id": o_id})
                m_id = self._match_material(doc)
                if not m_id:
                    raise ValueError("Cannot find matching material for object_id: {}".format(o_id))
                self._update_material(m_id, doc)
            except:
                import traceback
                logger.exception("<---")
                logger.exception("There was an error processing task_id: {}".format(o_id))
                logger.exception(traceback.format_exc())
                logger.exception("--->")

        logger.info("BoltztrapMaterialsBuilder finished processing.")

    def reset(self):
        logger.info("Resetting BoltztrapMaterialsBuilder")
        self._materials.update_many({}, {"$unset": {"_boltztrapbuilder": 1,
                                                    "transport": 1}})
        self._build_indexes()
        logger.info("Finished resetting BoltztrapMaterialsBuilder")

    def _match_material(self, doc, ltol=0.2, stol=0.3, angle_tol=5):
        """
        Returns the material_id that has the same structure as this doc as
         determined by the structure matcher. Returns None if no match.

        Args:
            doc (dict): a JSON-like document
            ltol (float): StructureMatcher tuning parameter
            stol (float): StructureMatcher tuning parameter
            angle_tol (float): StructureMatcher tuning parameter

        Returns:
            (int) matching material_id or None
        """
        formula = doc["formula_reduced_abc"]
        sgnum = doc["spacegroup"]["number"]

        for m in self._materials.find({"formula_reduced_abc": formula, "sg_number": sgnum},
                                      {"structure": 1, "material_id": 1}):

            m_struct = Structure.from_dict(m["structure"])
            t_struct = Structure.from_dict(doc["structure"])

            sm = StructureMatcher(ltol=ltol, stol=stol, angle_tol=angle_tol,
                                  primitive_cell=True, scale=True,
                                  attempt_supercell=False, allow_subset=False,
                                  comparator=ElementComparator())

            if sm.fit(m_struct, t_struct):
                return m["material_id"]

        return None

    def _update_material(self, m_id, doc):
        """
        Update a material document based on a new task

        Args:
            m_id (int): material_id for material document to update
            doc (dict): a JSON-like Boltztrap document
        """
        bta = BoltztrapAnalyzer.from_dict(doc)
        d = {}
        d["zt"] = bta.get_extreme("zt")
        d["pf"] = bta.get_extreme("power factor")
        d["seebeck"] = bta.get_extreme("seebeck")
        d["conductivity"] = bta.get_extreme("conductivity")
        d["kappa_max"] = bta.get_extreme("kappa")
        d["kappa_min"] = bta.get_extreme("kappa", maximize=False)

        self._materials.update_one({"material_id": m_id}, {"$set": {"transport": d}})
        self._materials.update_one({"material_id": m_id},
                                   {"$push": {"_boltztrapbuilder.all_object_ids": doc["_id"]}})

    def _build_indexes(self):
        """
        Create indexes for faster searching
        """
        for x in ["zt", "pf", "seebeck", "conductivity", "kappa_max", "kappa_min"]:
            self._materials.create_index("transport.{}.best.value".format(x))

    @classmethod
    def from_file(cls, db_file, m="materials", b="boltztrap", **kwargs):
        """
        Get a BoltztrapMaterialsBuilder using only a db file.

        Args:
            db_file (str): path to db file
            m (str): name of "materials" collection
            b (str): name of "boltztrap" collection
            **kwargs: other params to put into BoltztrapMaterialsBuilder
        """
        db_write = get_database(db_file, admin=True)
        try:
            db_read = get_database(db_file, admin=False)
            db_read.collection_names()  # throw error if auth failed
        except:
            print("Warning: could not get read-only database")
            db_read = get_database(db_file, admin=True)

        return cls(db_write[m], db_read[b], **kwargs)
