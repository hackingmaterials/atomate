from fireworks import Firework

from atomate.vasp.config import DB_FILE, VASP_CMD
from atomate.vasp.firetasks.approx_neb_dynamic_tasks import GetImageFireworks
from atomate.vasp.firetasks.approx_neb_tasks import AddSelectiveDynamics, PathfinderToDb

__author__ = "Ann Rutt"
__email__ = "acrutt@lbl.gov"


class EvaluatePathFW(Firework):
    def __init__(
        self,
        approx_neb_wf_uuid,
        end_points_combo,
        mobile_specie,
        n_images,
        selective_dynamics_scheme,
        launch_mode="all",
        db_file=DB_FILE,
        vasp_input_set=None,
        vasp_cmd=VASP_CMD,
        override_default_vasp_params=None,
        handler_group=None,
        parents=None,
        add_additional_fields=None,
        add_tags=None,
        **kwargs,
    ):
        r"""
        Applies NEBPathFinder (from pymatgen.analysis.path_finder)
        using the host charge density (chgcar from the task_id stored)
        and output structures stored in the "end_points" field of the
        approx_neb collection. Applies selective dynamics to the
        resulting image structures (interpolated between end point
        structures) and stores in the "images" field of the
        approx_neb collection for use in future calculations.
        Note see AddSelectiveDynamics Firetask for more information
        on "selective_dynamics_scheme". Launches VASP calculations
        using the specified launch_mode and adds task doc fields for
        approx_neb workflow record keeping. Stores relevant outputs
        in the approx_neb collection.

        Args:
        approx_neb_wf_uuid (str): unique id for approx neb
            workflow record keeping
        end_points_combo (str): string must have format of
            "0+1", "0+2", etc. to specify which combination
            of end_points to use for path interpolation.
        mobile_specie (str): specie of site of interest such
            as the working ion (e.g. "Li" if the working ion
            of interest is lithium). Used to perform a check
            on the structures pulled from the approx_neb
            collection.
        n_images: n_images (int): number of images
            interpolated between end point structures
        selective_dynamics_scheme (str): "fix_two_atom"
        launch_mode (str): "all" or "screening"
        db_file (str): path to file containing the database
            credentials.
        vasp_input_set (VaspInputSet class): can use to
            define VASP input parameters.
            See pymatgen.io.vasp.sets module for more
            information. MPRelaxSet() and
            override_default_vasp_params are used if
            vasp_input_set = None.
        vasp_cmd (str): the name of the full executable for running
            VASP.
        override_default_vasp_params (dict): if provided,
            vasp_input_set is disregarded and the Vasp Input
            Set is created by passing
            override_default_vasp_params to MPRelaxSet().
            Allows for easy modification of MPRelaxSet().
            For example, to set ISIF=2 in the INCAR use:
            {"user_incar_settings":{"ISIF":2}}
        handler_group (str or [ErrorHandler]): group of handlers to
            use for RunVaspCustodian firetask. See handler_groups
            dict in the code for the groups and complete list of
            handlers in each group. Alternatively, you can specify a
            list of ErrorHandler objects.
        parents ([Firework]): Parents of this particular Firework.
        add_additional_fields (dict): dict of additional fields to
            add to task docs (by additional_fields of VaspToDb).
        add_tags (list of strings): added to the "tags" field of the
            task docs.
        \*\*kwargs: Other kwargs that are passed to Firework.__init__.
        """
        fw_name = f"hop: {mobile_specie} {end_points_combo}"
        fw_spec = {"tags": ["approx_neb", approx_neb_wf_uuid, "evaluate_path"]}

        t = []
        # apply pathfinder pymatgen function and store outputs in approx_neb collection
        t.append(
            PathfinderToDb(
                db_file=db_file,
                n_images=n_images,
                end_points_combo=end_points_combo,
                approx_neb_wf_uuid=approx_neb_wf_uuid,
            )
        )
        # apply selective dynamics to pathfinder outputs to get images input structures
        t.append(
            AddSelectiveDynamics(
                approx_neb_wf_uuid=approx_neb_wf_uuid,
                pathfinder_key=end_points_combo,
                mobile_specie=mobile_specie,
                selective_dynamics_scheme=selective_dynamics_scheme,
                db_file=db_file,
            )
        )
        # add dynamic firetask that will launch image relaxations as desired
        t.append(
            GetImageFireworks(
                launch_mode=launch_mode,
                images_key=end_points_combo,
                approx_neb_wf_uuid=approx_neb_wf_uuid,
                vasp_cmd=vasp_cmd,
                db_file=db_file,
                vasp_input_set=vasp_input_set,
                override_default_vasp_params=override_default_vasp_params,
                handler_group=handler_group,
                add_additional_fields=add_additional_fields,
                add_tags=add_tags,
            )
        )

        super().__init__(tasks=t, spec=fw_spec, name=fw_name, parents=parents, **kwargs)
