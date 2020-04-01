# coding: utf-8

from __future__ import (
    absolute_import,
    division,
    print_function,
    unicode_literals,
)

import warnings

"""
Defines standardized Fireworks that can be chained easily to perform various
sequences of CP2K calculations.
"""

from fireworks.core.firework import Firework

from pymatgen.core.structure import Structure
from pymatgen.io.cp2k.sets import (
    RelaxSet,
    StaticSet,
    HybridStaticSet,
    HybridRelaxSet,
)

from atomate.cp2k.firetasks.write_inputs import (
    WriteCp2kFromIOSet,
    WriteCp2kFromPrevious,
)
from atomate.cp2k.firetasks.run_calc import RunCp2KCustodian
from atomate.cp2k.firetasks.glue_tasks import (
    UpdateStructureFromPrevCalc,
    CopyCp2kOutputs,
)
from atomate.cp2k.firetasks.parse_outputs import Cp2kToDb
from atomate.common.firetasks.glue_tasks import PassCalcLocs
from atomate.cp2k.utils import optimize_structure_sc_scale


# TODO: almost all Fws follow a similar copy from calc loc --> write input --> run --> pass locs/data.
#  Combine into single core firework maybe?


class BaseFW(Firework):
    def __init__(
        self,
        structure=None,
        name=None,
        cp2k_input_set=None,
        cp2k_input_set_params={},
        cp2k_cmd="cp2k",
        prev_calc_loc=True,
        prev_calc_dir=None,
        db_file=None,
        cp2ktodb_kwargs=None,
        parents=None,
        files_to_copy=[],
        **kwargs
    ):
        """
        Standard calculation Firework. A FW generally follows the following structure:
            (1) Update based on previous calculation
            (2) Copy files from previous calculation
            (3) Write input files as necessary
            (4) Run Cp2k with custodian
            (5) Pass this calculation location to the database (in order to stitch calculations together)
            (6) Run CP2K drone to assimilate results of the calculation

        And this doesn't change, whether its a FW that is part of series of calculations, or a parallel
        FW, or a standalone. Therefore, this BaseFW defines this general structure, and other FWs simply
        specify a slightly specific variation on this. Generally, the FWs that inherit will simply specify
        what the default input set is.

        Args:
            structure (Structure): Input structure. Note that for prev_calc_loc jobs, the structure
                is only used to set the name of the FW and any structure with the same composition
                can be used.
            name (str): Name for the Firework.
            cp2k_input_set (Cp2kInputSet): input set to use. Needed unless a previous calculation's inputs
                parameters will be used.
            cp2k_input_set_params (dict): Dict of Cp2kInputSet kwargs. Remember, overriding default set parameters
                is done inside of here using the "{override_default_params: {...} }" as an entry
            cp2k_cmd (str): Command to run cp2k. Supports env_chk.
            prev_calc_loc (bool or str): If true (default), copies outputs from previous calc. If
                a str value, retrieves a previous calculation output by name. If False/None, will create
                new calculation using the provided structure.
            prev_calc_dir (str): Path to a previous calculation to copy from.
            db_file (str): Path to file specifying db credentials. Supports env_chk.
            parents (Firework): Parents of this particular Firework. FW or list of FWS.
            cp2ktodb_kwargs (dict): kwargs to pass to Cp2kToDb
            \*\*kwargs: Other kwargs that are passed to Firework.__init__.
        """
        t = []

        # Set up the name, the input set, and how cp2k will be pushed to DB
        fw_name = "{}-{}".format(
            structure.composition.reduced_formula if structure else "unknown",
            name,
        )

        cp2k_input_set = cp2k_input_set or StaticSet(
            structure, **cp2k_input_set_params
        )

        cp2ktodb_kwargs = cp2ktodb_kwargs or {}
        if "additional_fields" not in cp2ktodb_kwargs:
            cp2ktodb_kwargs["additional_fields"] = {}
            cp2ktodb_kwargs["additional_fields"]["task_label"] = name

        # if continuing from prev calc, update the structure with the previous result
        # For hybrid, should almost always be true (initialize with gga first)
        if prev_calc_loc:
            t.append(
                CopyCp2kOutputs(
                    files_to_copy=files_to_copy, calc_loc=prev_calc_loc
                )
            )  # TODO START HERE WITH TESTING
            t.append(UpdateStructureFromPrevCalc(prev_calc_loc=prev_calc_loc))

        # if prev calc directory is being REPEATED, copy files
        if prev_calc_dir:
            t.append(
                WriteCp2kFromPrevious(
                    cp2k_input_set_params=cp2k_input_set_params
                )
            )

        # else run based on the IO set
        else:
            t.append(
                WriteCp2kFromIOSet(
                    structure=structure,
                    cp2k_input_set=cp2k_input_set,
                    cp2k_input_params=cp2k_input_set_params,
                )
            )

        t.append(RunCp2KCustodian(cp2k_cmd=cp2k_cmd))
        t.append(PassCalcLocs(name=name))
        t.append(Cp2kToDb(db_file=db_file, **cp2ktodb_kwargs))

        spec = {"cp2k_input_set": cp2k_input_set.as_dict()}
        spec.update(kwargs.get('spec', {}))
        super().__init__(
            t,
            parents=parents,
            name=fw_name,
            spec=spec,
            **kwargs
        )


class StaticFW(BaseFW):
    def __init__(
        self,
        structure=None,
        name="Static",
        cp2k_input_set=None,
        cp2k_input_set_params={},
        cp2k_cmd="cp2k",
        prev_calc_loc=True,
        prev_calc_dir=None,
        db_file=None,
        cp2ktodb_kwargs=None,
        parents=None,
        files_to_copy=[],
        **kwargs
    ):
        """
        Standard static calculation. Can start from input set or from previous calculation.
        """

        cp2k_input_set = cp2k_input_set or StaticSet(
            structure, project_name=name,
            **cp2k_input_set_params
        )

        super().__init__(
            structure=structure,
            name=name,
            cp2k_input_set=cp2k_input_set,
            cp2k_input_set_params=cp2k_input_set_params,
            cp2k_cmd=cp2k_cmd,
            prev_calc_loc=prev_calc_loc,
            prev_calc_dir=prev_calc_dir,
            db_file=db_file,
            cp2ktodb_kwargs=cp2ktodb_kwargs,
            parents=parents,
            files_to_copy=files_to_copy,
            **kwargs
        )


class StaticHybridFW(BaseFW):
    def __init__(
        self,
        structure=None,
        name="HybridStatic",
        cp2k_input_set=None,
        cp2k_input_set_params={},
        cp2k_cmd="cp2k",
        prev_calc_loc=True,
        prev_calc_dir=None,
        db_file=None,
        cp2ktodb_kwargs=None,
        parents=None,
        files_to_copy=[],
        **kwargs
    ):
        """
        Hybrid Static calculation. Presumably restarting from a previous GGA calculation.
        """

        cp2k_input_set = cp2k_input_set or HybridStaticSet(
            structure, project_name=name,
            **cp2k_input_set_params
        )

        super().__init__(
            structure=structure,
            name=name,
            cp2k_input_set=cp2k_input_set,
            cp2k_input_set_params=cp2k_input_set_params,
            cp2k_cmd=cp2k_cmd,
            prev_calc_loc=prev_calc_loc,
            prev_calc_dir=prev_calc_dir,
            db_file=db_file,
            cp2ktodb_kwargs=cp2ktodb_kwargs,
            parents=parents,
            files_to_copy=files_to_copy,
            **kwargs
        )


class RelaxFW(BaseFW):
    def __init__(
        self,
        structure=None,
        name="Relax",
        cp2k_input_set=None,
        cp2k_input_set_params={},
        cp2k_cmd="cp2k",
        prev_calc_loc=True,
        prev_calc_dir=None,
        db_file=None,
        cp2ktodb_kwargs=None,
        parents=None,
        files_to_copy=[],
        **kwargs
    ):
        """
        Geometry Relaxation Calculation.
        """

        cp2k_input_set = cp2k_input_set or RelaxSet(
            structure, project_name=name,
            **cp2k_input_set_params
        )

        super().__init__(
            structure=structure,
            name=name,
            cp2k_input_set=cp2k_input_set,
            cp2k_input_set_params=cp2k_input_set_params,
            cp2k_cmd=cp2k_cmd,
            prev_calc_loc=prev_calc_loc,
            prev_calc_dir=prev_calc_dir,
            db_file=db_file,
            cp2ktodb_kwargs=cp2ktodb_kwargs,
            parents=parents,
            files_to_copy=files_to_copy,
            **kwargs
        )


class RelaxHybridFW(BaseFW):
    def __init__(
        self,
        structure=None,
        name="HybridRelax",
        cp2k_input_set=None,
        cp2k_input_set_params={},
        cp2k_cmd="cp2k",
        prev_calc_loc=True,
        prev_calc_dir=None,
        db_file=None,
        cp2ktodb_kwargs=None,
        parents=None,
        files_to_copy=[],
        **kwargs
    ):
        """
        Hybrid Static calculation. Presumably restarting from a previous GGA calculation.
        """

        cp2k_input_set = cp2k_input_set or HybridRelaxSet(
            structure, project_name=name,
            **cp2k_input_set_params
        )

        super().__init__(
            structure=structure,
            name=name,
            cp2k_input_set=cp2k_input_set,
            cp2k_input_set_params=cp2k_input_set_params,
            cp2k_cmd=cp2k_cmd,
            prev_calc_loc=prev_calc_loc,
            prev_calc_dir=prev_calc_dir,
            db_file=db_file,
            cp2ktodb_kwargs=cp2ktodb_kwargs,
            parents=parents,
            files_to_copy=files_to_copy,
            **kwargs
        )
