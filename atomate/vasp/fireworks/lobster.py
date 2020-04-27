import logging
import os
from typing import List, Union

from atomate.common.firetasks.glue_tasks import DeleteFiles, PassCalcLocs, DeleteFilesPrevFolder
from atomate.vasp.config import DB_FILE
from atomate.vasp.firetasks.glue_tasks import CopyVaspOutputs
from atomate.vasp.firetasks.lobster_tasks import WriteLobsterinputfromIO, RunLobster, LobsterRunToDb
from fireworks import Firework
from pymatgen.core.structure import Structure
from custodian.custodian import ErrorHandler, Validator


LOBSTER_CMD = ">>lobster_cmd<<"
VASP_OUTPUT_FILES_without_PPKI = ["OUTCAR", "vasprun.xml", "CHG", "CHGCAR", "CONTCAR", "DOSCAR", "EIGENVAL", "IBZKPT",
                                  "OSZICAR", "PCDAT", "PROCAR", "REPORT", "WAVECAR", "XDATCAR"]

MODULE_DIR = os.path.dirname(os.path.abspath(__file__))
logger = logging.getLogger(__name__)


class LobsterFW(Firework):
    """
    A Firework performs a Lobster calculation based on a previous static VASP calculation with specific configurations
    """

    def __init__(self, structure: Structure = None, name: str = "lobster_calculation", lobster_cmd: str = LOBSTER_CMD,
                 db_file: str = DB_FILE, delete_wavecar: bool = False,
                 delete_wavecar_previous_fw: bool = False,
                 handler_group: Union[List[ErrorHandler], str] = "default",
                 validator_group: Union[List[Validator], str] = "default",
                 calculationtype: str = 'standard',
                 parents: Union[List[Firework], Firework] = None,
                 prev_calc_dir: str = None, prev_calc_loc: bool = True, user_supplied_basis: dict = None,
                 lobsterin_key_dict: dict = None, lobstertodb_kwargs: dict = None,
                 additional_outputs: List[str] = None, **kwargs):
        """

        Args:
            structure (Structure): Structure object, will only be used to name firework
            name (str): name of the firework
            lobster_cmd (str): command to run lobster
            db_file (str): address to db_file
            delete_wavecar (bool): If True, WAVECAR will be deleted
            delete_wavecar_previous_fw (bool): If True, WAVECAR from VASP calc will be deleted
            handler_group (Union[List[ErrorHandler],str])): group of handlers to use. See handler_groups dict in the code for
                the groups and complete list of handlers in each group. Alternatively, you can
                specify a list of ErrorHandler objects.
            validator_group (Union[List[Validator],str]): group of validators to use. See validator_groups dict in the
                code for the groups and complete list of validators in each group. Alternatively, you can
            specify a list of Validator objects.
            calculationtype (str): only 'standard' is fully implemented so far
            parents (Union[List[Firework],Firework]): parent Firework
            prev_calc_dir (str): address to previous vasp calculation
            prev_calc_loc (bool): If true, calc wil be started from previous directory
            user_supplied_basis (dict): the user can supply their own basis functions
            lobsterin_key_dict (dict): the user can supply additional changes to the lobsterin with {"COHPendEnergy":10.0}
            lobstertodb_kwargs (dict): dict that will be saved in the mongodb database
            additional_outputs (List[str]): list of additional files to be stored in the
                results DB. They will be stored as files in gridfs. Examples are:
                "ICOHPLIST.lobster" or "DOSCAR.lobster". Note that the file name
                should be given with the full name and the correct capitalization.
            **kwargs:
        """

        # TODO: make this lobster firework more flexible to allow for FATBAND and other types of calculations

        fw_name = "{}-{}".format(structure.composition.reduced_formula if structure else 'unknown', name)
        
        t = []
        # copies all files from previous VASP calculation;
        if prev_calc_dir:
            t.append(
                CopyVaspOutputs(calc_dir=prev_calc_dir, additional_files=VASP_OUTPUT_FILES_without_PPKI,
                                contcar_to_poscar=False))
            calc_loc = False
        elif parents:
            t.append(
                CopyVaspOutputs(calc_loc=prev_calc_loc, additional_files=VASP_OUTPUT_FILES_without_PPKI,
                                contcar_to_poscar=False))
            calc_loc = True
        else:
            raise ValueError("One must specify a VASP calculation for Lobster run")

        t.append(WriteLobsterinputfromIO(poscar_path='POSCAR', incar_path='INCAR', potcar_path='POTCAR',
                                         option=calculationtype, user_supplied_basis=user_supplied_basis,
                                         user_lobsterin_settings=lobsterin_key_dict))

        # runs lobster
        if delete_wavecar:
            t.append(RunLobster(lobster_cmd=lobster_cmd, gzip_output=True, gzip_WAVECAR=False,
                                handler_group=handler_group, validator_group=validator_group))
        else:
            t.append(RunLobster(lobster_cmd=lobster_cmd, gzip_output=True, gzip_WAVECAR=True,
                                handler_group=handler_group, validator_group=validator_group))

        # task to delete wavecar to avoid storage problems with wavecars -> WAVECARs without symmetry can be huge
        if delete_wavecar:
            t.append(DeleteFiles(files=["WAVECAR", "WAVECAR.gz"]))
        if delete_wavecar_previous_fw:
            t.append(DeleteFilesPrevFolder(files=["WAVECAR", "WAVECAR.gz"], calc_loc=calc_loc))

        # Will save Lobster Calculation in Database
        t.append(LobsterRunToDb(db_file=db_file, calc_loc=calc_loc,
                                additional_fields=lobstertodb_kwargs,
                                additional_outputs=additional_outputs))

        # passes information to next firetask
        t.append(PassCalcLocs(name=name))

        super(LobsterFW, self).__init__(t, parents=parents, name=fw_name, **kwargs)
