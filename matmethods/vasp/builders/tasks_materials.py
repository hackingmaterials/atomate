import os

from pymongo import ReturnDocument

from monty.serialization import loadfn
from pymatgen import Structure
from pymatgen.analysis.structure_matcher import StructureMatcher, \
    ElementComparator

__author__ = 'Anubhav Jain <ajain@lbl.gov>'

module_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)))

# TODO: !! add a real counter collection
# TODO: make this work in parallel for better performance - watch for race conditions w/same formula+spacegroup combo
# TODO: support nested keys in settings for tasks_key (probably some easy way to do this via monty or something)

class TasksMaterialsBuilder:
    def __init__(self, tasks_read, materials_write, counter_write):
        self._tasks = tasks_read
        self._materials = materials_write
        self._counter = counter_write

        if self._counter.find({"_id": "materialid"}).count() == 0:
            self._counter.insert_one({"_id": "materialid", "c": 0})

        x = loadfn(os.path.join(module_dir, "tasks_materials_settings.yaml"))
        self.property_settings = x['property_settings']
        # TODO: add option to give prefix to task-ids when building materials

    def run(self):

        print("Initiating list of all previously processed task_ids ...")
        previous_task_ids = []
        for m in self._materials.find({}, {"_tmbuilder.all_task_ids": 1}):
            previous_task_ids.extend(m["_tmbuilder"]["all_task_ids"])

        if not previous_task_ids:
            # TODO: make this part of the YAML settings file
            self._materials.create_index("_tmbuilder.all_task_ids")
            self._materials.create_index("material_id", unique=True)

        print("Initiating list of all successful task_ids ...")
        all_task_ids = [t["task_id"] for t in self._tasks.find({"state": "successful"}, {"task_id": 1})]

        print("Determining list of new task_ids ...")
        task_ids = [t_id for t_id in all_task_ids if t_id not in previous_task_ids]
        print("There are {} new task_ids to process".format(len(task_ids)))

        for t_id in task_ids:
            print("Processing task id: {}".format(t_id))
            taskdoc = self._tasks.find_one({"task_id": t_id})
            self._preprocess_taskdoc(taskdoc)

            m_id = self._match_material(taskdoc)
            if not m_id:
                m_id = self._create_new_material(taskdoc)

            self._update_material(m_id, taskdoc)

    def reset(self):
        self._materials.delete_many({})
        self._counter.delete_one({"_id": "materialid"})
        self._counter.insert_one({"_id": "materialid", "c": 0})

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
                                       "spacegroup.number": sgnum},
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
        doc = {}
        doc["_tmbuilder"] = {"all_task_ids": [], "prop_metadata":
            {"labels": {}, "task_ids": {}}}
        doc["formula_reduced_abc"] = taskdoc["formula_reduced_abc"]
        doc["spacegroup"] = taskdoc["output"]["spacegroup"]
        doc["structure"] = taskdoc["output"]["structure"]
        doc["material_id"] = self._counter.find_one_and_update(
                        {"_id": "materialid"}, {"$inc": {"c": 1}},
                        return_document=ReturnDocument.AFTER)["c"]
        self._materials.insert_one(doc)

        return doc["material_id"]

    def _update_material(self, m_id, taskdoc):
        # add task_id to list of all_task_ids
        self._materials.update_one({"material_id": m_id},
                                   {"$push": {"_tmbuilder.all_task_ids":
                                                  taskdoc["task_id"]}})

        # get list of labels for each property
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
                        self._materials.\
                            update_one({"material_id": m_id},
                                       {"$set": {materials_key: taskdoc[x["tasks_key"]][p],
                                                 "_tmbuilder.prop_metadata.labels.{}".format(p): task_label,
                                                 "_tmbuilder.prop_metadata.task_ids.{}".format(p): taskdoc["task_id"]}})