# coding: utf-8

from __future__ import absolute_import, division, print_function, unicode_literals

import os
from datetime import datetime
from pymongo import ReturnDocument
from tqdm import tqdm

from matgendb.util import get_database

from monty.serialization import loadfn

from pymatgen import Structure
from pymatgen.analysis.structure_matcher import StructureMatcher, ElementComparator

from atomate.utils.utils import get_mongolike, get_logger
from atomate.vasp.builders.base import AbstractBuilder

logger = get_logger(__name__)

__author__ = 'Anubhav Jain <ajain@lbl.gov>'

module_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)))


class TasksMaterialsBuilder(AbstractBuilder):
    def __init__(self, materials_write, counter_write, tasks_read, tasks_prefix="t",
                 materials_prefix="m", query=None):
        """
        Create a materials collection from a tasks collection.

        Args:
            materials_write (pymongo.collection): mongodb collection for materials (write access needed)
            counter_write (pymongo.collection): mongodb collection for counter (write access needed)
            tasks_read (pymongo.collection): mongodb collection for tasks (suggest read-only for safety)
            tasks_prefix (str): a string prefix for tasks, e.g. "t" gives a task_id like "t-132"
            materials_prefix (str): a string prefix to prepend to material_ids
            query (dict): a pymongo query on tasks_read for which tasks to include in the builder
        """
        x = loadfn(os.path.join(module_dir, "tasks_materials_settings.yaml"))
        self.supported_task_labels = x['supported_task_labels']
        self.property_settings = x['property_settings']
        self.indexes = x.get('indexes', [])
        self.properties_root = x.get('properties_root', [])

        self._materials = materials_write
        if self._materials.count() == 0:
            self._build_indexes()

        self._counter = counter_write
        if self._counter.find({"_id": "materialid"}).count() == 0:
            self._counter.insert_one({"_id": "materialid", "c": 0})

        self._tasks = tasks_read
        self._t_prefix = tasks_prefix
        self._m_prefix = materials_prefix
        self.query = query

    def run(self):
        logger.info("MaterialsTaskBuilder starting...")
        logger.info("Initializing list of all new task_ids to process ...")
        previous_task_ids = []
        for m in self._materials.find({}, {"_tasksbuilder.all_task_ids": 1}):
            previous_task_ids.extend(m["_tasksbuilder"]["all_task_ids"])

        q = {}
        q["state"] = "successful"
        q["task_label"] = {"$in": self.supported_task_labels}

        if self.query:
            common_keys = [k for k in q.keys() if k in self.query.keys()]
            if common_keys:
                raise ValueError("User query parameter cannot contain key(s): {}".
                                 format(common_keys))
            q.update(self.query)

        all_task_ids = [self.tid_to_str(t["task_id"]) for t in self._tasks.find(q, {"task_id": 1})]
        task_ids = [t_id for t_id in all_task_ids if t_id not in previous_task_ids]

        logger.info("There are {} new task_ids to process.".format(len(task_ids)))

        pbar = tqdm(task_ids)
        for t_id in pbar:
            pbar.set_description("Processing task_id: {}".format(t_id))
            try:
                taskdoc = self._tasks.find_one({"task_id": self.tid_to_int(t_id)})
                m_id = self._match_material(taskdoc)
                if not m_id:
                    m_id = self._create_new_material(taskdoc)
                self._update_material(m_id, taskdoc)

            except:
                import traceback
                logger.exception("<---")
                logger.exception("There was an error processing task_id: {}".format(t_id))
                logger.exception(traceback.format_exc())
                logger.exception("--->")

        logger.info("TasksMaterialsBuilder finished processing.")

    def reset(self):
        logger.info("Resetting TasksMaterialsBuilder")
        self._materials.delete_many({})
        self._counter.delete_one({"_id": "materialid"})
        self._counter.insert_one({"_id": "materialid", "c": 0})
        self._build_indexes()
        logger.info("Finished resetting TasksMaterialsBuilder")

    def tid_to_str(self, task_id):
        # converts int task_id to string
        return "{}-{}".format(self._t_prefix, task_id)

    @staticmethod
    def tid_to_int(task_id):
        # converts string task_id to int
        return int(task_id.split("-")[1])

    def mid_to_str(self, material_id):
        # converts int material_id to string
        return "{}-{}".format(self._m_prefix, material_id)

    @classmethod
    def from_file(cls, db_file, m="materials", c="counter", t="tasks", **kwargs):
        """
        Get a TaskMaterialsBuilder using only a db file.

        Args:
            db_file (str): path to db file
            m (str): name of "materials" collection
            c (str): name of "counter" collection
            t (str): name of "tasks" collection
            **kwargs: other params to put into TasksMaterialsBuilder
        """
        db_write = get_database(db_file, admin=True)
        try:
            db_read = get_database(db_file, admin=False)
            db_read.collection_names()  # throw error if auth failed
        except:
            logger.warn("Warning: could not get read-only database; using write creds")
            db_read = get_database(db_file, admin=True)
        return cls(db_write[m], db_write[c], db_read[t], **kwargs)

    def _build_indexes(self):
        """
        Create indexes for faster searching
        """
        self._materials.create_index("material_id", unique=True)
        for index in self.indexes:
            self._materials.create_index(index)

    def _match_material(self, taskdoc, ltol=0.2, stol=0.3, angle_tol=5):
        """
        Returns the material_id that has the same structure as this task as
         determined by the structure matcher. Returns None if no match.

        Args:
            taskdoc (dict): a JSON-like task document
            ltol (float): StructureMatcher tuning parameter 
            stol (float): StructureMatcher tuning parameter 
            angle_tol (float): StructureMatcher tuning parameter

        Returns:
            (int) matching material_id or None
        """
        formula = taskdoc["formula_reduced_abc"]

        # handle the "parent structure" option, which is used to intentionally force slightly
        # different structures to contribute to the same "material", e.g. from an ordering scheme
        if "parent_structure" in taskdoc:
            t_struct = Structure.from_dict(taskdoc["parent_structure"]["structure"])
            q = {"formula_reduced_abc": formula, "parent_structure.spacegroup.number": taskdoc[
                "parent_structure"]["spacegroup"]["number"]}
        else:
            sgnum = taskdoc["output"]["spacegroup"]["number"]
            t_struct = Structure.from_dict(taskdoc["output"]["structure"])
            q = {"formula_reduced_abc": formula, "sg_number": sgnum}

        for m in self._materials.find(q, {"parent_structure": 1, "structure": 1, "material_id": 1}):
            s_dict = m["parent_structure"]["structure"] if "parent_structure" in m else m[
                "structure"]
            m_struct = Structure.from_dict(s_dict)
            sm = StructureMatcher(ltol=ltol, stol=stol, angle_tol=angle_tol,
                                  primitive_cell=True, scale=True,
                                  attempt_supercell=False, allow_subset=False,
                                  comparator=ElementComparator())

            if sm.fit(m_struct, t_struct):
                return m["material_id"]

        return None

    def _create_new_material(self, taskdoc):
        """
        Create a new material document.

        Args:
            taskdoc (dict): a JSON-like task document

        Returns:
            (int) - material_id of the new document
        """
        doc = {"created_at": datetime.utcnow()}
        doc["_tasksbuilder"] = {"all_task_ids": [], "prop_metadata":
            {"labels": {}, "task_ids": {}}, "updated_at": datetime.utcnow()}
        doc["spacegroup"] = taskdoc["output"]["spacegroup"]
        doc["sg_symbol"] = doc["spacegroup"]["symbol"]
        doc["sg_number"] = doc["spacegroup"]["number"]
        doc["structure"] = taskdoc["output"]["structure"]
        doc["material_id"] = self.mid_to_str(self._counter.
                                             find_one_and_update({"_id": "materialid"},
                                                                 {"$inc": {"c": 1}},
                                                                 return_document=ReturnDocument.AFTER)[
                                                 "c"])

        for x in ["formula_anonymous", "formula_pretty", "formula_reduced_abc",
                  "elements", "nelements", "chemsys"]:
            doc[x] = taskdoc[x]

        if "parent_structure" in taskdoc:
            doc["parent_structure"] = taskdoc["parent_structure"]
            t_struct = Structure.from_dict(taskdoc["parent_structure"]["structure"])
            doc["parent_structure"]["formula_reduced_abc"] = t_struct.composition.reduced_formula

        self._materials.insert_one(doc)

        return doc["material_id"]

    def _update_material(self, m_id, taskdoc):
        """
        Update a material document based on a new task and using complex logic

        Args:
            m_id (int): material_id for material document to update
            taskdoc (dict): a JSON-like task document
        """

        # For each materials property, figure out what kind of task the data is currently based on
        # as defined by the task label.  This is used to decide if the new taskdoc is a type of
        # calculation that provides higher quality data for that property
        prop_tlabels = self._materials.find_one(
            {"material_id": m_id}, {"_tasksbuilder.prop_metadata.labels": 1})[
            "_tasksbuilder"]["prop_metadata"]["labels"]

        task_label = taskdoc["task_label"]  # task label of new doc that updates this material

        # figure out what materials properties need to be updated based on new task
        for x in self.property_settings:
            for p in x["properties"]:
                # check if this is a valid task for getting the property
                if task_label in x["quality_scores"]:
                    # assert: this is a valid task for the property
                    # but is it the "best" task for that property (highest quality score)?
                    t_quality = x["quality_scores"][task_label]
                    m_quality = x["quality_scores"].get(prop_tlabels.get(p, None), None)
                    # check if this task's quality is better than existing data
                    # 3 possibilities:
                    # i) materials property data not present, so this is best
                    # ii) task quality higher based on task label
                    # iii) task quality equal to materials; use lowest energy task
                    if not m_quality or t_quality > m_quality \
                            or (t_quality == m_quality
                                and taskdoc["output"]["energy_per_atom"] <
                                    self._materials.find_one({"material_id": m_id}, {
                                        "_tasksbuilder": 1})["_tasksbuilder"]["prop_metadata"][
                                        "energies"][p]):

                        # this task has better quality data
                        # figure out where the property data lives in the materials doc and
                        # in the task doc
                        materials_key = "{}.{}".format(x["materials_key"], p) \
                            if x.get("materials_key") else p
                        tasks_key = "{}.{}".format(x["tasks_key"], p) \
                            if x.get("tasks_key") else p

                        # insert property data AND metadata about this task
                        self._materials.\
                            update_one({"material_id": m_id},
                                       {"$set": {materials_key: get_mongolike(taskdoc, tasks_key),
                                                 "_tasksbuilder.prop_metadata.labels.{}".format(p): task_label,
                                                 "_tasksbuilder.prop_metadata.task_ids.{}".format(p): self.tid_to_str(taskdoc["task_id"]),
                                                 "_tasksbuilder.prop_metadata.energies.{}".format(p): taskdoc["output"]["energy_per_atom"],
                                                 "_tasksbuilder.updated_at": datetime.utcnow()}})

                        # copy property to document root if in properties_root
                        # i.e., intentionally duplicate some data to the root level
                        if p in self.properties_root:
                            self._materials.\
                            update_one({"material_id": m_id},
                                       {"$set": {p: get_mongolike(taskdoc, tasks_key)}})

        self._materials.update_one({"material_id": m_id},
                                   {"$push": {"_tasksbuilder.all_task_ids": self.tid_to_str(taskdoc["task_id"])}})
