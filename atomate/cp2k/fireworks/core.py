# coding: utf-8

from __future__ import absolute_import, division, print_function, \
    unicode_literals

import warnings

"""
Defines standardized Fireworks that can be chained easily to perform various
sequences of CP2K calculations.
"""

from fireworks import Firework

from pymatgen import Structure
from pymatgen.io.cp2k.sets import RelaxSet, StaticSet, HybridStaticSet, HybridRelaxSet

from atomate.cp2k.firetasks.write_inputs import WriteCp2kFromIOSet
from atomate.cp2k.firetasks.run_calc import RunCp2KCustodian


class StaticFW(Firework):

    def __init__(self, structure=None, name="static", cp2k_input_set=None,
                 cp2k_input_set_params={}, cp2k_cmd="vasp", prev_calc_loc=True,
                 prev_calc_dir=None, db_file=None, cp2ktodb_kwargs=None,
                 parents=None, **kwargs):
        """
        Standard static calculation Firework - either from a previous location or from a structure.

        Args:
            structure (Structure): Input structure. Note that for prev_calc_loc jobs, the structure
                is only used to set the name of the FW and any structure with the same composition
                can be used.
            name (str): Name for the Firework.
            vasp_input_set (VaspInputSet): input set to use (for jobs w/no parents)
                Defaults to MPStaticSet() if None.
            vasp_input_set_params (dict): Dict of vasp_input_set kwargs.
            vasp_cmd (str): Command to run vasp.
            prev_calc_loc (bool or str): If true (default), copies outputs from previous calc. If
                a str value, retrieves a previous calculation output by name. If False/None, will create
                new static calculation using the provided structure.
            prev_calc_dir (str): Path to a previous calculation to copy from
            db_file (str): Path to file specifying db credentials.
            parents (Firework): Parents of this particular Firework. FW or list of FWS.
            vasptodb_kwargs (dict): kwargs to pass to VaspToDb
            \*\*kwargs: Other kwargs that are passed to Firework.__init__.
        """
        t = []

        fw_name = "{}-{}".format(structure.composition.reduced_formula if structure else "unknown", name)

        cp2k_input_set = cp2k_input_set or StaticSet(structure)
        t.append(WriteCp2kFromIOSet(structure=structure,
                                    cp2k_input_set=cp2k_input_set,
                                    cp2k_input_set_params=cp2k_input_set_params))

        t.append(RunCp2KCustodian(cp2k_cmd=cp2k_cmd))

        super(StaticFW, self).__init__(t, parents=parents, name=fw_name, **kwargs)
