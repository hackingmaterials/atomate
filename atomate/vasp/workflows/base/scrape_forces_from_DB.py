#!/usr/bin/env python3


import numpy as np
from pymongo import MongoClient


# MongoDB database settings/details
host = "brmyvxghlvvkgkx-mongodb.services.clever-cloud.com"
port = 27017
database_name = "brmyvxghlvvkgkx"
username = "usxrmfzwh9ed1plwi027"
password = "JJGXekEfnNYtXQSBYHv8"

# Connect to the database
db = MongoClient(host, port)[database_name]
db.authenticate(username, password)
tasks = db.tasks #connects to the 'tasks' collection

"""
This module defines IO ShengBTE for making CONTROL file (i.e. a class that can read/write the ‘control’ file, f90ml)

TO-DO: Test this thing with csld_main
"""

def scrape_forces_to_txt(task_id):
    # Input a Firetask task_id (int)
    # Output a .txt with all the atomic forces

    forces_list = tasks.find_one({'task_id': task_id})['calcs_reversed'][0]['output']['ionic_steps'][0]['forces']
    # print(forces_list)
    # print(len(forces_list))

    num_atoms = len(forces_list)
    forces = np.empty((num_atoms,3))
    for atom in range(num_atoms):
        forces[atom,:] = forces_list[atom][:]
    # print(forces)
    np.savetxt('forces.txt', forces, fmt='%.6f')


def main():
    scrape_forces_to_txt(2)


if __name__ == '__main__':
    main()

# def get_structure_from_task_id(task_id):
#     # Input a Firetask task_id (int)
#     # Return the corresponding Pymatgen structure (Structure)
#     try:
#         return Structure.from_dict(tasks.find_one({"task_id": task_id})["input"]["structure"])
#     except:
#         print("No structure with task id '{}' was found.".format(task_id))
#
#
#
# def get_structures_from_formulaV2(formula, task_label=None):
#     # Required:
#     #   Input a Pymatgen 'pretty formula' (String)
#     # Optional:
#     #   Input a 'task_label' or substring of a 'task_label' (String)
#     # Return the corresponding Pymatgen Structures (list)
#
#     structs = []
#     try:
#         query = {"formula_pretty": formula}
#         if task_label is not None and isinstance(task_label, str):
#             # query["$and"] += [{"task_label": task_label}]
#             query["task_label"] = {"$regex": task_label}
#         dicts = list(tasks.find(query, ["input.structure"]))
#
#         # dicts = list(tasks.find({"$and": [{"formula_pretty": formula},
#         #                                   {"task_label": task_label}]},
#         #                                 ["input.structure"]))
#         for dict in dicts:
#             struct = Structure.from_dict(dict["input"]["structure"])
#             structs += [struct]
#         return structs
#     except:
#         if isinstance(formula, str):
#             print("No structures with the formula '{}' found.".format(formula))
#         else:
#             print("Entered formula is not a string.")
