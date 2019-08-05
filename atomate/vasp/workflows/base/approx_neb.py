from fireworks import Firework, Workflow
from copy import deepcopy
from atomate.vasp.config import VASP_CMD, DB_FILE
from atomate.vasp.powerups import use_custodian
from custodian.vasp.handlers import (
    VaspErrorHandler,
    MeshSymmetryErrorHandler,
    PotimErrorHandler,
    FrozenJobErrorHandler,
    NonConvergingErrorHandler,
    PositiveEnergyErrorHandler,
    StdErrHandler,
)
from uuid import uuid4

from atomate.vasp.fireworks.approx_neb import (
    HostLatticeFW,
    InsertSitesFW,
    ApproxNEBLaunchFW,
    PathFinderFW,
    GetImagesFW,
)

# ToDo: add way to provide tags/additional fields to docs in approx_neb collection
# TODO: Write approx_neb_wf_description
def approx_neb_wf(
    structure,
    working_ion,
    insert_coords,
    n_images,
    vasp_input_set=None,
    override_default_vasp_params=None,
    selective_dynamics_scheme="fix_two_atom",
    vasp_cmd=VASP_CMD,
    db_file=DB_FILE,
    name="Approx NEB",
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
            "LAECHG": True
        }
    }
    # TODO: Add LASPH: True
    # ToDo: Add "LAECHG": True, to all or just host lattice?

    wf_uuid = str(uuid4())

    host_lattice_fw = HostLatticeFW(
        structure=structure,
        approx_neb_wf_uuid=wf_uuid,
        db_file=db_file,
        vasp_input_set=vasp_input_set,
        vasp_cmd=vasp_cmd,
        override_default_vasp_params=deepcopy(approx_neb_params),
    )

    # modifies incar settings needed for end point and image structure relaxations
    if "user_incar_settings" not in approx_neb_params.keys():
        approx_neb_params = {"user_incar_settings": {}}
    approx_neb_params["user_incar_settings"]["ISIF"] = 2
    approx_neb_params["user_incar_settings"]["ISYM"] = 0
    approx_neb_params["user_incar_settings"]["LDAU"] = False

    if len(insert_coords) != 2:
        raise ValueError("insert_coords must specify exactly two sites")

    insert_working_ion_fws = []
    for n,coord in enumerate(insert_coords):
        insert_working_ion_fws.append(
            InsertSitesFW(
                approx_neb_wf_uuid=wf_uuid,
                insert_specie=working_ion,
                insert_coords=coord,
                stable_sites_index=n,
                db_file=db_file,
                parents=host_lattice_fw,
            )
        )

    stable_site_fws = []
    for fw in insert_working_ion_fws:
        stable_site_fws.append(
            ApproxNEBLaunchFW(
                calc_type="stable_site",
                approx_neb_wf_uuid=wf_uuid,
                db_file=db_file,
                override_default_vasp_params=approx_neb_params,
                parents=fw,
            )
        )

    pathfinder_fw = PathFinderFW(
        approx_neb_wf_uuid=wf_uuid,
        n_images=n_images,
        db_file=db_file,
        parents=stable_site_fws,
    )

    get_images_fw = GetImagesFW(
        approx_neb_wf_uuid=wf_uuid,
        mobile_specie=working_ion,
        selective_dynamics_scheme=selective_dynamics_scheme,
        parents=pathfinder_fw
    )

    relax_image_fws = []
    for n in range(0,n_images):
        path = "images." + str(n) + ".input_structure"
        relax_image_fws.append(
            ApproxNEBLaunchFW(
                calc_type="image",
                approx_neb_wf_uuid=wf_uuid,
                structure_path=path,
                db_file=db_file,
                override_default_vasp_params=approx_neb_params,
                parents=get_images_fw,
            )
        )

    wf = Workflow(
        [host_lattice_fw]
        + insert_working_ion_fws
        + stable_site_fws
        + [pathfinder_fw]
        + [get_images_fw]
        + relax_image_fws
    )

    wf = use_custodian(
        wf,
        custodian_params={
            "handler_group": [
                VaspErrorHandler(),
                MeshSymmetryErrorHandler(),
                NonConvergingErrorHandler(),
                PotimErrorHandler(),
                PositiveEnergyErrorHandler(),
                FrozenJobErrorHandler(),
                StdErrHandler(),
            ]
        },
    )
    wf.name = name

    return wf


def approx_neb_screening_wf(
    structure,
    working_ion,
    insert_coords,
    n_images,
    insert_coords_combinations,
    vasp_input_set=None,
    override_default_vasp_params=None,
    selective_dynamics_scheme="fix_two_atom",
    vasp_cmd=VASP_CMD,
    db_file=DB_FILE,
    name="Approx NEB",
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
            "LAECHG": True
        }
    }
    # TODO: Add LASPH: True
    # ToDo: Add "LAECHG": True, to all or just host lattice?

    wf_uuid = str(uuid4())

    host_lattice_fw = HostLatticeFW(
        structure=structure,
        approx_neb_wf_uuid=wf_uuid,
        db_file=db_file,
        vasp_input_set=vasp_input_set,
        vasp_cmd=vasp_cmd,
        override_default_vasp_params=deepcopy(approx_neb_params),
    )

    # modifies incar settings needed for end point and image structure relaxations
    if "user_incar_settings" not in approx_neb_params.keys():
        approx_neb_params = {"user_incar_settings": {}}
    approx_neb_params["user_incar_settings"]["ISIF"] = 2
    approx_neb_params["user_incar_settings"]["ISYM"] = 0
    approx_neb_params["user_incar_settings"]["LDAU"] = False

    insert_working_ion_fws = []
    for n, coord in enumerate(insert_coords):
        insert_working_ion_fws.append(
            InsertSitesFW(
                approx_neb_wf_uuid=wf_uuid,
                insert_specie=working_ion,
                insert_coords=coord,
                stable_sites_index=n,
                db_file=db_file,
                parents=host_lattice_fw,
            )
        )

    stable_site_fws = []
    for fw in insert_working_ion_fws:
        stable_site_fws.append(
            ApproxNEBLaunchFW(
                calc_type="stable_site",
                approx_neb_wf_uuid=wf_uuid,
                db_file=db_file,
                override_default_vasp_params=approx_neb_params,
                parents=fw,
            )
        )

    pathfinder_fws = []
    for stable_sites_combo in insert_coords_combinations:
        if isinstance(stable_sites_combo, (str)):
            combo = stable_sites_combo.split("+")
            if len(combo) == 2:
                c = [int(combo[0]), int(combo[-1])]
            else:
                ValueError("string format in insert_coords_combinations is incorrect")

        pathfinder_fws.append(
            PathFinderFW(
                approx_neb_wf_uuid=wf_uuid,
                n_images=n_images,
                stable_sites_combo = combo,
                db_file=db_file,
                parents=[stable_site_fws[c[0]],stable_site_fws[c[1]]]
            )
        )

    get_images_fw = GetImagesFW(
        approx_neb_wf_uuid=wf_uuid,
        mobile_specie=working_ion,
        selective_dynamics_scheme=selective_dynamics_scheme,
        parents=pathfinder_fws
    )

    #relax_image_fws = [] #ToDo: modify
    #for n in range(0,n_images):
    #    path = "images." + str(n) + ".input_structure"
    #    relax_image_fws.append(
    #        ApproxNEBLaunchFW(
    #            calc_type="image",
    #            approx_neb_wf_uuid=wf_uuid,
    #            structure_path=path,
    #            db_file=db_file,
    #            override_default_vasp_params=approx_neb_params,
    #            parents=get_images_fw,
    #        )
    #    )

    wf = Workflow(
        [host_lattice_fw]
        + insert_working_ion_fws
        + stable_site_fws
        + pathfinder_fws
        + [get_images_fw]
        #+ relax_image_fws
    )

    wf = use_custodian(
        wf,
        custodian_params={
            "handler_group": [
                VaspErrorHandler(),
                MeshSymmetryErrorHandler(),
                NonConvergingErrorHandler(),
                PotimErrorHandler(),
                PositiveEnergyErrorHandler(),
                FrozenJobErrorHandler(),
                StdErrHandler(),
            ]
        },
    )
    wf.name = name

    return wf