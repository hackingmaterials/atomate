from pymatgen import Structure
from pymatgen.analysis.structure_matcher import StructureMatcher, \
    ElementComparator
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

__author__ = 'Anubhav Jain <ajain@lbl.gov>'

class TaskMaterialsBuilder:
    COUNTER = 0

    def __init__(self, tasks_read, materials_write):
        self._tasks = tasks_read
        self._materials = materials_write

        self.labels_qualities = {}
        self.labels_qualities["structure optimization"] = {"bandgap": 1, "energy_per_atom": 1}
        self.labels_qualities["static"] = {"bandgap": 2, "energy_per_atom": 2}
        self.labels_qualities["nscf uniform"] = {"bandgap": 3, "energy_per_atom": 0}
        self.labels_qualities["nscf line"] = {"bandgap": 4, "energy_per_atom": 0}

        # TODO: add option to give prefix to task-ids when building materials

    def run(self):

        print("Initiating list of all previously processed task_ids ...")
        previous_task_ids = []
        for m in self._materials.find({}, {"all_sources.task_ids": 1}):
            previous_task_ids.extend(m["all_sources"]["task_ids"])

        if not previous_task_ids:
            self._materials.create_index("all_sources.task_ids")

        print("Initiating list of all successful task_ids ...")
        all_task_ids = [t["task_id"] for t in self._tasks.find({"state": "successful"}, {"task_id": 1})]

        print("Determining list of new task_ids ...")

        task_ids = [t_id for t_id in all_task_ids if t_id not in previous_task_ids]
        print("There are {} new task_ids to process".format(len(task_ids)))

        for t_id in task_ids:
            taskdoc = self._tasks.find_one({"task_id": t_id})
            self._preprocess_taskdoc(taskdoc)
            m_id = self._match_material(taskdoc)

            if m_id:
                self._add_task_to_material(m_id, taskdoc)
            else:
                m_id = self._add_new_material(taskdoc)

            self._refresh_material(m_id, taskdoc)

    def _preprocess_taskdoc(self, taskdoc):

        # handles backward incompatibilities
        taskdoc["output"]["spacegroup"]["number"] = int(taskdoc["output"]["spacegroup"]["number"])

    def _match_material(self, taskdoc):
        formula = taskdoc["formula_reduced_abc"]
        sgnum = taskdoc["output"]["spacegroup"]["number"]

        for m in self._materials.find({"formula_reduced_abc": formula, "spacegroup.number": sgnum}, {"structure": 1, "materials_id": 1}):
            m_struct = Structure.from_dict(m["structure"])
            t_struct = Structure.from_dict(taskdoc["output"]["structure"])

            sm = StructureMatcher(ltol=0.2, stol=0.3, angle_tol=5, primitive_cell=True,
                 scale=True, attempt_supercell=False, allow_subset=False,
                 comparator=ElementComparator())

            if sm.fit(m_struct, t_struct):
                return m["materials_id"]

        return None

    def _add_task_to_material(self, m_id, taskdoc):
        self._materials.update_one({"materials_id": m_id}, {"$push": {"all_sources.task_ids": taskdoc["task_id"]}})


    def _add_new_material(self, taskdoc):
        self.COUNTER += 1
        doc = {}
        doc["all_sources"] = {"task_ids": [taskdoc["task_id"]]}
        doc["formula_reduced_abc"] = taskdoc["formula_reduced_abc"]
        doc["spacegroup"] = taskdoc["output"]["spacegroup"]
        doc["structure"] = taskdoc["output"]["structure"]
        doc["props_labels"] = {}
        doc["materials_id"] = self.COUNTER
        self._materials.insert_one(doc)

        return self.COUNTER

    def _refresh_material(self, m_id, taskdoc):
        pl_material = self._materials.find_one({"materials_id": m_id}, {"props_labels": 1})["props_labels"]
        if taskdoc["task_label"] in self.labels_qualities:
            pq_task = self.labels_qualities[taskdoc["task_label"]]
        else:
            print("WARNING: don't know how to handle task label: {}".format(taskdoc["task_label"]))
            pq_task = {}

        for p, q in pq_task.iteritems():
            if p not in pl_material or q > self.labels_qualities[pl_material[p]][p]:
                self._materials.update_one({"materials_id": m_id}, {"$set": {p: taskdoc["output"][p], "props_labels.{}".format(p): taskdoc["task_label"]}})


