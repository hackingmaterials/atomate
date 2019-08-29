from fireworks import Firework, Workflow
from copy import deepcopy
from atomate.vasp.config import VASP_CMD, DB_FILE
from atomate.vasp.powerups import use_custodian, add_tags, add_additional_fields_to_taskdocs
from custodian.vasp.handlers import (
    VaspErrorHandler,
    MeshSymmetryErrorHandler,
    PotimErrorHandler,
    FrozenJobErrorHandler,
    NonConvergingErrorHandler,
    PositiveEnergyErrorHandler,
    StdErrHandler,
    WalltimeHandler
)
from uuid import uuid4

from atomate.vasp.fireworks.approx_neb import (
    HostLatticeFW,
    ApproxNEBLaunchFW,
    StableSiteFW
)
from atomate.vasp.fireworks.approx_neb_dynamic import EvaluatePathFW

# TODO: Write approx_neb_wf_description

def approx_neb_wf(
    structure,
    working_ion,
    insert_coords,
    insert_coords_combinations,
    n_images,
    vasp_input_set=None,
    override_default_vasp_params=None,
    handler_group = None,
    selective_dynamics_scheme="fix_two_atom",
    launch_mode="all",
    vasp_cmd=VASP_CMD,
    db_file=DB_FILE,
    wall_time=None,
    additional_fields = None,
    tags = None,
    name="ApproxNEB",
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
            "LAECHG": True
        }
    }
    # TODO: Add LASPH: True
    # ToDo: Add "LAECHG": True, to all or just host lattice?
    handler_group = handler_group or [
            VaspErrorHandler(),
            MeshSymmetryErrorHandler(),
            NonConvergingErrorHandler(),
            PotimErrorHandler(),
            PositiveEnergyErrorHandler(),
            FrozenJobErrorHandler(),
            StdErrHandler(),
            WalltimeHandler(wall_time=wall_time)
        ]

    wf_uuid = str(uuid4())
    additional_fields = deepcopy(additional_fields)

    host_lattice_fw = HostLatticeFW(
        structure=structure,
        approx_neb_wf_uuid=wf_uuid,
        db_file=db_file,
        vasp_input_set=vasp_input_set,
        vasp_cmd=vasp_cmd,
        override_default_vasp_params=deepcopy(approx_neb_params),
        additional_fields = additional_fields,
        tags = tags
    )

    # modifies incar settings needed for stable site and image structure relaxations
    if "user_incar_settings" not in approx_neb_params.keys():
        approx_neb_params = {"user_incar_settings": {}}
    approx_neb_params["user_incar_settings"]["ISIF"] = 2
    approx_neb_params["user_incar_settings"]["ISYM"] = 0
    approx_neb_params["user_incar_settings"]["LDAU"] = False

    stable_site_fws = []
    for n, coord in enumerate(insert_coords):
        stable_site_fws.append(
            StableSiteFW(
                approx_neb_wf_uuid=wf_uuid,
                insert_specie=working_ion,
                insert_coords=coord,
                stable_sites_index=n,
                db_file=db_file,
                override_default_vasp_params=approx_neb_params,
                parents=host_lattice_fw,
            )
        )

    evaluate_path_fws = []
    for stable_sites_combo in insert_coords_combinations:
        if isinstance(stable_sites_combo, (str)):
            combo = stable_sites_combo.split("+")
            if len(combo) == 2:
                c = [int(combo[0]), int(combo[-1])]
            else:
                raise ValueError("string format in insert_coords_combinations is incorrect")

        evaluate_path_fws.append(
            EvaluatePathFW(
                approx_neb_wf_uuid=wf_uuid,
                stable_sites_combo = stable_sites_combo,
                mobile_specie=working_ion,
                n_images=n_images,
                selective_dynamics_scheme=selective_dynamics_scheme,
                launch_mode=launch_mode,
                vasp_cmd=vasp_cmd,
                db_file=db_file,
                override_default_vasp_params=approx_neb_params,
                handler_group=handler_group,
                parents=[stable_site_fws[c[0]],stable_site_fws[c[1]]],
                add_additional_fields=additional_fields,
                add_tags=tags
            )
        )

    wf = Workflow([host_lattice_fw] + stable_site_fws + evaluate_path_fws)

    wf = use_custodian(
        wf,
        custodian_params={
            "handler_group": handler_group
        },
    )
    if isinstance(tags,(list)):
        wf = add_tags(wf, tags)
    wf = add_additional_fields_to_taskdocs(wf, update_dict=additional_fields)
    wf.metadata.update({"approx_neb_wf_uuid":wf_uuid})
    wf.name = name

    return wf