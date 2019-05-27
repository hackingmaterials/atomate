from fireworks import Firework, FWAction, Workflow, FiretaskBase
from atomate.vasp.firetasks.glue_tasks import pass_vasp_result
from atomate.vasp.firetasks.run_calc import RunVaspCustodian
from atomate.vasp.fireworks.core import OptimizeFW
from atomate.vasp.config import VASP_CMD, DB_FILE

from atomate.vasp.workflows.base.approx_neb_wf_parts import *

#TODO: Write approx_neb_wf_description
def approx_neb_wf(
    structure,
    working_ion,
    insert_coords,
    n_images,
    vasp_input_set=None,
    override_default_vasp_params=None,
    vasp_cmd=VASP_CMD,
    db_file=DB_FILE,
):
    approx_neb_params = override_default_vasp_params or {
        "user_incar_settings": {
            "EDIFF": 0.0005,
            "EDIFFG": -0.05,
            "IBRION": 1,
            "ISIF": 3,
            "ISMEAR": 0,
            "LDAU": False,
            "NSW": 400,
            "ADDGRID": True,
            "ISYM": 1,
            "NELMIN": 4,
        }
    }
    host_lattice_fw = OptimizeFW(
        structure,
        name="approx neb optimize host lattice",
        vasp_input_set=vasp_input_set,
        override_default_vasp_params=approx_neb_params.copy(),
        vasp_cmd=vasp_cmd,
        db_file=db_file,
        vasptodb_kwargs={"parse_chgcar": True, "parse_aeccar": True},
    )

    pass_host_lattice_fw = pass_vasp_result(
        pass_dict={
            "chgcar_file_path": "???",
        }
    )  # best way to pass chgcar from host lattice to PathfinderFW???

    if "user_incar_settings" not in approx_neb_params.keys():
        approx_neb_params = {"user_incar_settings": {}}
    approx_neb_params["user_incar_settings"]["ISIF"] = 2
    approx_neb_params["user_incar_settings"]["ISYM"] = 0

     # firework of single firetask to pass output structure...

    insert_working_ion_fws = []
    for coord in insert_coords:
        insert_working_ion_fws.append(
            InsertSitesFW(
                insert_specie=working_ion,
                insert_coords=coord,
                vasp_input_set=vasp_input_set,
                override_default_vasp_params=approx_neb_params.copy(),
                vasp_cmd=vasp_cmd,
                db_file=db_file,
                parents=host_lattice_fw,
            )
        )
        # spec={}
        # list of fireworks where one ion is inserted - assume two coords provided for now, fws can run independently

    pathfinder_fws = PathFinderFW(
        ep1_struct="???",
        ep2_struct="???",
        n_images=n_images,
        chgcar="???",
        vasp_input_set=vasp_input_set,
        override_default_vasp_params=approx_neb_params.copy(),
        vasp_cmd=vasp_cmd,
        db_file=db_file,
        parents=[host_lattice_fw] + insert_working_ion_fws,
    )
        # kwarg for working ion index? default = 0 in selective dynamics function
        # list of fireworks for all images

    wf = Workflow(
        [host_lattice_fw] + [pass_host_lattice_fw] + insert_working_ion_fws
    )
    #TODO: modify workflow to remove undesirable custodian handlers
    #TODO: need unique id for workflow?
    return wf