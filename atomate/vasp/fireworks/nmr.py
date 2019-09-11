# coding: utf-8

from fireworks import Firework

from pymatgen.io.vasp.sets import MPNMRSet

from atomate.common.firetasks.glue_tasks import PassCalcLocs
from atomate.vasp.firetasks.glue_tasks import CopyVaspOutputs
from atomate.vasp.firetasks.parse_outputs import VaspToDb
from atomate.vasp.firetasks.run_calc import RunVaspCustodian
from atomate.vasp.firetasks.write_inputs import WriteVaspNMRFromPrev, WriteVaspFromIOSet


class NMRFW(Firework):
    def __init__(self,
                 structure=None,
                 mode="cs",
                 isotopes=None,
                 name="nmr tensor",
                 prev_calc_dir=None,
                 vasp_cmd="vasp",
                 copy_vasp_outputs=True,
                 db_file=None,
                 parents=None,
                 **kwargs):
        """
        Firework for NMR tensor calculations

        Args:
            structure (Structure): Input structure. If copy_vasp_outputs, used only to set the
                name of the FW.
            mode (str): the NMR calculation type: cs or efg, default is cs
            isotopes (list): list of isotopes to include, default is to include the
                             lowest mass quadrupolar isotope for all applicable elements
            name (str): Name for the Firework.
            prev_calc_dir (str): Path to a previous calculation to copy from
            vasp_cmd (str): Command to run vasp.
            copy_vasp_outputs (bool): Whether to copy outputs from previous
                run. Defaults to True.
            db_file (str): Path to file specifying db credentials.
            parents (Firework): Parents of this particular Firework.
                FW or list of FWS.
            kwargs: Other kwargs that are passed to Firework.__init__.
        """
        fw_name = "{}-{}".format(structure.composition.reduced_formula if structure else "unknown", name)

        isotopes = isotopes.split() if isinstance(isotopes, str) else isotopes
        t = []
        if prev_calc_dir:
            t.append(CopyVaspOutputs(calc_dir=prev_calc_dir, contcar_to_poscar=True))
            t.append(WriteVaspNMRFromPrev(prev_calc_dir=".", mode=mode, isotopes=isotopes))
        elif parents and copy_vasp_outputs:
            t.append(CopyVaspOutputs(calc_loc=True, contcar_to_poscar=True))
            t.append(WriteVaspNMRFromPrev(prev_calc_dir=".", mode=mode, isotopes=isotopes))
        elif structure:
            vasp_input_set = MPNMRSet(structure, mode=mode, isotopes=isotopes)
            t.append(WriteVaspFromIOSet(structure=structure, vasp_input_set=vasp_input_set))
        else:
            raise ValueError("Must specify structure or previous calculation.")

        t.extend([
            RunVaspCustodian(vasp_cmd=vasp_cmd, auto_npar=">>auto_npar<<"),
            PassCalcLocs(name=name),
            VaspToDb(db_file=db_file, additional_fields={"task_label": name})
        ])
        super(NMRFW, self).__init__(t, parents=parents, name=fw_name, **kwargs)