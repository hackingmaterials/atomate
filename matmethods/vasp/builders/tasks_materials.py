from __future__ import absolute_import, print_function

import os
from datetime import datetime

from pymongo import ReturnDocument
from tqdm import tqdm

from matgendb.util import get_database
from matmethods.utils.utils import get_mongolike
from monty.serialization import loadfn
from pymatgen import Structure
from pymatgen.analysis.structure_matcher import StructureMatcher, ElementComparator

__author__ = 'Anubhav Jain <ajain@lbl.gov>'

module_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)))

# TODO: make this work in parallel for better performance - watch for race conditions w/same formula+spacegroup combo
# TODO: if multiple entries with same quality, choose lowest energy one. quality score can be a tuple then

class TasksMaterialsBuilder:
    def __init__(self, materials_write, counter_write, tasks_read,
                 tasks_prefix="t", materials_prefix="m"):
        """
        Create a materials collection from a tasks collection
        Args:
            materials_write: mongodb collection for materials (write access needed)
            counter_write: mongodb collection for counter (write access needed)
            tasks_read: mongodb collection for tasks (suggest read-only for safety)
        """
        x = loadfn(os.path.join(module_dir, "tasks_materials_settings.yaml"))
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

    def tid_str(self, task_id):
        # converts int material_id to string
        return "{}-{}".format(self._t_prefix, task_id)

    def tid_int(self, task_id):
        # converts string task_id to int
        return int(task_id.replace(self._t_prefix+"-", ""))

    def mid_str(self, material_id):
        # converts int material_id to string
        return "{}-{}".format(self._m_prefix, material_id)

    def run(self):
        print("MaterialsTaskBuilder starting...")
        print("Initializing list of all new task_ids to process ...")
        previous_task_ids = []
        for m in self._materials.find({}, {"_tmbuilder.all_task_ids": 1}):
            previous_task_ids.extend(m["_tmbuilder"]["all_task_ids"])

        all_task_ids = [self.tid_str(t["task_id"]) for t in self._tasks.find({"state": "successful"}, {"task_id": 1})]
        task_ids = [t_id for t_id in all_task_ids if t_id not in previous_task_ids]

        print("There are {} new task_ids to process.".format(len(task_ids)))

        pbar = tqdm(task_ids)
        for t_id in pbar:
            pbar.set_description("Processing task_id: {}".format(t_id))
            try:
                taskdoc = self._tasks.find_one({"task_id": self.tid_int(t_id)})
                self._preprocess_taskdoc(taskdoc)  # TODO: move pre-process to separate builder

                m_id = self._match_material(taskdoc)
                if not m_id:
                    m_id = self._create_new_material(taskdoc)
                self._update_material(m_id, taskdoc)

            except:
                import traceback
                print("<---")
                print("There was an error processing task_id: {}".format(t_id))
                traceback.print_exc()
                print("--->")

        print("TasksMaterialsBuilder finished processing.")

    def reset(self):
        self._materials.delete_many({})
        self._counter.delete_one({"_id": "materialid"})
        self._counter.insert_one({"_id": "materialid", "c": 0})
        self._build_indexes()

    @staticmethod
    def _preprocess_taskdoc(taskdoc):
        """
        Preprocess a task doc, usually as a way to handle backwards
        incompatibilities and refactorings that accumulate over time

        Args:
            taskdoc: A JSON-like task document
        """

        # cast spacegroup number to int type
        taskdoc["output"]["spacegroup"]["number"] = \
            int(taskdoc["output"]["spacegroup"]["number"])

    def _match_material(self, taskdoc):
        """
        Returns the material_id that has the same structure as this task as
         determined by the structure matcher. Returns None if no match.

        Args:
            taskdoc: a JSON-like task document

        Returns:
            (int) matching material_id or None
        """
        formula = taskdoc["formula_reduced_abc"]
        sgnum = taskdoc["output"]["spacegroup"]["number"]

        for m in self._materials.find({"formula_reduced_abc": formula,
                                       "sg_number": sgnum},
                                      {"structure": 1, "material_id": 1}):

            m_struct = Structure.from_dict(m["structure"])
            t_struct = Structure.from_dict(taskdoc["output"]["structure"])

            sm = StructureMatcher(ltol=0.2, stol=0.3, angle_tol=5,
                                  primitive_cell=True, scale=True,
                                  attempt_supercell=False, allow_subset=False,
                                  comparator=ElementComparator())

            if sm.fit(m_struct, t_struct):
                return m["material_id"]

        return None

    def _create_new_material(self, taskdoc):
        """
        Create a new material document
        Args:
            taskdoc: a JSON-like task document

        Returns:
            (int) - material_id of the new document
        """
        doc = {"created_at": datetime.utcnow()}
        doc["_tmbuilder"] = {"all_task_ids": [], "prop_metadata":
            {"labels": {}, "task_ids": {}}, "updated_at": datetime.utcnow()}
        doc["spacegroup"] = taskdoc["output"]["spacegroup"]
        doc["sg_symbol"] = doc["spacegroup"]["symbol"]
        doc["sg_number"] = doc["spacegroup"]["number"]
        doc["structure"] = taskdoc["output"]["structure"]
        doc["material_id"] = self.mid_str(self._counter.find_one_and_update(
                        {"_id": "materialid"}, {"$inc": {"c": 1}},
                        return_document=ReturnDocument.AFTER)["c"])
        for x in ["formula_anonymous", "formula_pretty", "formula_reduced_abc",
                  "nelements"]:
            doc[x] = taskdoc[x]

        self._materials.insert_one(doc)

        return doc["material_id"]

    def _update_material(self, m_id, taskdoc):
        """
        Update a material document based on a new task

        Args:
            m_id: (int) material_id for material document to update
            taskdoc: a JSON-like task document
        """
        # get list of labels for each existing property in material
        # this is used to decide if the taskdoc has higher quality data
        x = self._materials.find_one({"material_id": m_id},
                                            {"_tmbuilder.prop_metadata.labels":
                                                1})
        m_labels = x["_tmbuilder"]["prop_metadata"]["labels"]

        task_label = taskdoc["task_label"]
        # figure out what properties need to be updated
        for x in self.property_settings:
            for p in x["properties"]:
                # check if this is a valid task for getting the property
                if task_label in x["quality_scores"]:
                    t_quality = x["quality_scores"][task_label]
                    m_label = m_labels.get(p, None)
                    m_quality = x["quality_scores"].get(m_label, None)
                    # check if this task's quality is better than existing data
                    if not m_quality or t_quality > m_quality:
                        # insert task's properties into material
                        materials_key = "{}.{}".format(x["materials_key"], p) \
                            if x.get("materials_key") else p
                        tasks_key = "{}.{}".format(x["tasks_key"], p) \
                            if x.get("tasks_key") else p

                        self._materials.\
                            update_one({"material_id": m_id},
                                       {"$set": {materials_key: get_mongolike(taskdoc, tasks_key),
                                                 "_tmbuilder.prop_metadata.labels.{}".format(p): task_label,
                                                 "_tmbuilder.prop_metadata.task_ids.{}".format(p): self.tid_str(taskdoc["task_id"]),
                                                 "_tmbuilder.updated_at": datetime.utcnow()}})

                        # copy property to document root if in properties_root
                        if p in self.properties_root:
                            self._materials.\
                            update_one({"material_id": m_id},
                                       {"$set": {p: get_mongolike(taskdoc,
                                                                  tasks_key)}})

        
        self._materials.update_one({"material_id": m_id},
                                   {"$push": {"_tmbuilder.all_task_ids":
                                                  self.tid_str(taskdoc["task_id"])}})

    @staticmethod
    def from_db_file(db_file, m="materials", c="counter", t="tasks", **kwargs):
        """
        Get a TaskMaterialsBuilder using only a db file
        Args:
            db_file: (str) path to db file
            m: (str) name of "materials" collection
            c:  (str) name of "counter" collection
            t:  (str) name of "tasks" collection
            **kwargs: other params to put into TasksMaterialsBuilder
        """
        db_write = get_database(db_file, admin=True)
        try:
            db_read = get_database(db_file, admin=False)
        except:
            print("Warning: could not get read-only database")
            db_read = get_database(db_file, admin=True)

        return TasksMaterialsBuilder(db_write[m], db_write[c], db_read[t], **kwargs)


    def _build_indexes(self):
        """
        Create indexes for faster searching
        """
        self._materials.create_index("material_id", unique=True)
        for index in self.indexes:
            self._materials.create_index(index)
