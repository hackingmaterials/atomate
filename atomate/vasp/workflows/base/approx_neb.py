from fireworks import Firework, Workflow
from atomate.vasp.fireworks.core import OptimizeFW
from atomate.vasp.config import VASP_CMD, DB_FILE
from uuid import uuid4

from atomate.vasp.firetasks.approx_neb_tasks import HostLatticeToDb
from atomate.vasp.fireworks.approx_neb import InsertSitesFW, ApproxNEBLaunchFW

# TODO: Write approx_neb_wf_description
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
            "NSW": 200,
            "ADDGRID": True,
            "ISYM": 1,
            "NELMIN": 4,
        }
    }

    host_lattice_fw = OptimizeFW(
        structure,
        name="approx neb host lattice optimization",
        vasp_input_set=vasp_input_set,
        override_default_vasp_params=approx_neb_params.copy(),
        vasp_cmd=vasp_cmd,
        db_file=db_file,
        vasptodb_kwargs={
            "parse_chgcar": True,
            "parse_aeccar": True,
            "additional_fields": {
                "approx_neb.py": {"calc_type": "host_lattice", "wf_uuids": []}
            },
            "task_fields_to_push": {"host_lattice_task_id": "task_id"},
        },
    )
    # TODO: check if passed host_lattice_fw task_id correctly

    wf_uuid = str(uuid4())
    # TODO: how to check that fw_spec has host_lattice_task_id
    initialize_db_fw = Firework(
        tasks=[HostLatticeToDb(db_file=DB_FILE, approx_neb_wf_uuid=wf_uuid)],
        parents=host_lattice_fw,
    )

    if "user_incar_settings" not in approx_neb_params.keys():
        approx_neb_params = {"user_incar_settings": {}}
    approx_neb_params["user_incar_settings"]["ISIF"] = 2
    approx_neb_params["user_incar_settings"]["ISYM"] = 0

    # firework of single firetask to pass output structure...

    insert_working_ion_fws = []
    for coord in insert_coords:
        insert_working_ion_fws.append(
            InsertSitesFW(
                approx_neb_wf_uuid=wf_uuid,
                insert_specie=working_ion,
                insert_coords=coord,
                db_file=db_file,
                parents=initialize_db_fw,
            )
        )

    stable_site_fws = []
    for fw in insert_working_ion_fws:
        stable_site_fws.append(
            ApproxNEBLaunchFW(
                calc_type="stable_site", approx_neb_wf_uuid=wf_uuid, parents=fw
            )
        )
    # pathfinder_fws = PathFinderFW(
    #    ep1_struct="???",
    #    ep2_struct="???",
    #    n_images=n_images,
    #    chgcar="???",
    #    vasp_input_set=vasp_input_set,
    #    override_default_vasp_params=approx_neb_params.copy(),
    #    vasp_cmd=vasp_cmd,
    #    db_file=db_file,
    #    parents=[host_lattice_fw] + insert_working_ion_fws,
    # )
    # list of fireworks for all images

    wf = Workflow(
        [host_lattice_fw]
        + [initialize_db_fw]
        + insert_working_ion_fws
        + stable_site_fws
    )
    # TODO: modify workflow to remove undesirable custodian handlers

    return wf


# TODO: store WI in approx_neb.py collection
