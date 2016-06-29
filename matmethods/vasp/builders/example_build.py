from matmethods.vasp.builders.boltztrap_materials import \
    BoltztrapMaterialsBuilder
from matmethods.vasp.builders.fix_tasks import FixTasksBuilder
from matmethods.vasp.builders.materials_ehull import MaterialsEhullBuilder
from matmethods.vasp.builders.tasks_materials import TasksMaterialsBuilder

__author__ = 'Anubhav Jain <ajain@lbl.gov>'

"""
This is an example build script - just set the vars below and run
"""

DB_FILE = "PATH/TO/DB.JSON"  # path to valid db.json file
MAPI_KEY = None  # set this str if you don't have MAPI_KEY env var set

if __name__ == "__main__":
    FixTasksBuilder.from_db_file(DB_FILE).run()
    TasksMaterialsBuilder.from_db_file(DB_FILE).run()
    MaterialsEhullBuilder.from_db_file(DB_FILE, mapi_key=MAPI_KEY).run()
    BoltztrapMaterialsBuilder.from_db_file(DB_FILE).run()