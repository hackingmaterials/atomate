# coding: utf-8

from __future__ import absolute_import, division, print_function, unicode_literals

from atomate.vasp.builders.boltztrap_materials import BoltztrapMaterialsBuilder
from atomate.vasp.builders.fix_tasks import FixTasksBuilder
from atomate.vasp.builders.materials_ehull import MaterialsEhullBuilder
from atomate.vasp.builders.tasks_materials import TasksMaterialsBuilder

__author__ = 'Anubhav Jain <ajain@lbl.gov>'

"""
This is an example build script - just set the vars below and run
"""

DB_FILE = "PATH/TO/DB.JSON"  # path to valid db.json file
MAPI_KEY = None  # set this str if you don't have MAPI_KEY env var set

if __name__ == "__main__":
    FixTasksBuilder.from_file(DB_FILE).run()
    TasksMaterialsBuilder.from_file(DB_FILE).run()
    MaterialsEhullBuilder.from_file(DB_FILE, mapi_key=MAPI_KEY).run()
    BoltztrapMaterialsBuilder.from_file(DB_FILE).run()
    # add more Builders as desired using the same format...
    # use .reset() method of a Builder to start over
