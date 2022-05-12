from fireworks.core.launchpad import LaunchPad
from pymatgen.core import Molecule

from atomate.qchem.workflows.base.FF_then_fragment import get_wf_FF_then_fragment

mol = Molecule.from_file("BF4-.xyz")
wf = get_wf_FF_then_fragment(molecule=mol, max_cores=32)
lp = LaunchPad.auto_load()
lp.add_wf(wf)
