"""
To use this file, first modify db.json located in this directory with details for your atomate
output database.
"""
import os

from atomate.vasp.builders.bandgap_estimation import BandgapEstimationBuilder
from atomate.vasp.builders.dielectric import DielectricBuilder
from atomate.vasp.builders.fix_tasks import FixTasksBuilder
from atomate.vasp.builders.materials_descriptor import MaterialsDescriptorBuilder
from atomate.vasp.builders.tags import TagsBuilder
from atomate.vasp.builders.tasks_materials import TasksMaterialsBuilder

__author__ = 'Anubhav Jain <ajain@lbl.gov>'

module_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)))

if __name__ == "__main__":

    dbfile = os.path.join(module_dir, "db.json")  # make sure to modify w/your db details

    build_sequence = [FixTasksBuilder, TasksMaterialsBuilder, TagsBuilder,
                      MaterialsDescriptorBuilder, BandgapEstimationBuilder, DielectricBuilder]
    for cls in build_sequence:
        b = cls.from_file(dbfile)
        # b.reset()  # uncomment if you want to start from a builder from scratch
        b.run()