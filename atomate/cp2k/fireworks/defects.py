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

