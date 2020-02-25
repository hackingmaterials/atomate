from typing import Optional, List, Union

from atomate.common.firetasks.glue_tasks import PassCalcLocs
from atomate.vasp.firetasks.lattice_dynamics import CollectPerturbedStructures, \
    FitForceConstants, ForceConstantsToDB
from atomate.vasp.workflows.base.lattice_dynamics import IMAGINARY_TOL, \
    MAX_N_IMAGINARY, MAX_IMAGINARY_FREQ, FIT_METHOD, MESH_DENSITY
from fireworks import Firework

__author__ = "Alex Ganose"
__email__ = "aganose@lbl.gov"


class FitForceConstantsFW(Firework):

    def __init__(
        self,
        wf_uuid,
        name="fit force constants",
        parents: Union[Firework, List[Firework]] = None,
        db_file: str = None,
        cutoffs: Optional[List[List[float]]] = None,
        imaginary_tol: float = IMAGINARY_TOL,
        max_n_imaginary: int = MAX_N_IMAGINARY,
        max_imaginary_freq: float = MAX_IMAGINARY_FREQ,
        fit_method: str = FIT_METHOD,
        mesh_density: float = MESH_DENSITY,
        additional_fields: dict = None,
        **kwargs
    ):
        """
        Compile perturbed supercell calculations and fit force constants
        using hiPhive.

        Args:
            wf_uuid: Workflow identifier, from which the perturbed supercell
                static calculations will be compiled.
            parents: Parent(s) of this Firework.
            name: Name of this FW.
            db_file: Path to a db file.
            cutoffs: A list of cutoffs to trial. If None, a set of trial cutoffs
                will be generated based on the structure (default).
            imaginary_tol: Tolerance used to decide if a phonon mode is
                imaginary, in THz.
            max_n_imaginary: Maximum number of imaginary modes allowed in the
                final fitted force constant solution. If this criteria is not
                reached by any cutoff combination the Firework will fizzle.
            max_imaginary_freq: Maximum allowed imaginary frequency in the
                final fitted force constant solution. If this criteria is not
                reached by any cutoff combination this FireTask will fizzle.
            fit_method: Method used for fitting force constants. This can be
                any of the values allowed by the hiphive ``Optimizer`` class.
            mesh_density: The density of the q-point mesh used to calculate the
                phonon density of states. See the docstring for the ``mesh``
                argument in Phonopy.init_mesh() for more details.
            additional_fields: Additional fields added to the document, such as
                user-defined tags, name, ids, etc.
            **kwargs: Other kwargs that are passed to Firework.__init__.
        """
        fw_name = "{}-{}".format(wf_uuid, name)
        additional_fields = additional_fields or {}

        collect_structures = CollectPerturbedStructures(db_file, wf_uuid)
        fit_constants = FitForceConstants(
            cutoffs=cutoffs,
            imaginary_tol=imaginary_tol,
            max_n_imaginary=max_n_imaginary,
            max_imaginary_freq=max_imaginary_freq,
            fit_method=fit_method
        )
        fc_to_db = ForceConstantsToDB(
            db_file=db_file,
            wf_uuid=wf_uuid,
            mesh_density=mesh_density,
            additional_fields=additional_fields
        )
        pass_locs = PassCalcLocs(name=name)

        tasks = [collect_structures, fit_constants, fc_to_db, pass_locs]
        super().__init__(tasks, parents=parents, name=fw_name, **kwargs)


class LatticeThermalConductivityFW(Firework):

    def __init__(
            self,
            name="fit force constants",
            parents: Union[Firework, List[Firework]] = None,
            db_file: str = None,
            cutoffs: Optional[List[List[float]]] = None,
            imaginary_tol: float = IMAGINARY_TOL,
            max_n_imaginary: int = MAX_N_IMAGINARY,
            max_imaginary_freq: float = MAX_IMAGINARY_FREQ,
            fit_method: str = FIT_METHOD,
            mesh_density: float = MESH_DENSITY,
            additional_fields: dict = None,
            **kwargs
    ):
        """
        Compile perturbed supercell calculations and fit force constants
        using hiPhive.

        Args:
            wf_uuid: Workflow identifier, from which the perturbed supercell
                static calculations will be compiled.
            parents: Parent(s) of this Firework.
            name: Name of this FW.
            db_file: Path to a db file.
            cutoffs: A list of cutoffs to trial. If None, a set of trial cutoffs
                will be generated based on the structure (default).
            imaginary_tol: Tolerance used to decide if a phonon mode is
                imaginary, in THz.
            max_n_imaginary: Maximum number of imaginary modes allowed in the
                final fitted force constant solution. If this criteria is not
                reached by any cutoff combination the Firework will fizzle.
            max_imaginary_freq: Maximum allowed imaginary frequency in the
                final fitted force constant solution. If this criteria is not
                reached by any cutoff combination this FireTask will fizzle.
            fit_method: Method used for fitting force constants. This can be
                any of the values allowed by the hiphive ``Optimizer`` class.
            mesh_density: The density of the q-point mesh used to calculate the
                phonon density of states. See the docstring for the ``mesh``
                argument in Phonopy.init_mesh() for more details.
            additional_fields: Additional fields added to the document, such as
                user-defined tags, name, ids, etc.
            **kwargs: Other kwargs that are passed to Firework.__init__.
        """
        fw_name = "{}-{}".format(wf_uuid, name)
        additional_fields = additional_fields or {}

        collect_structures = CollectPerturbedStructures(db_file, wf_uuid)
        fit_constants = FitForceConstants(
            cutoffs=cutoffs,
            imaginary_tol=imaginary_tol,
            max_n_imaginary=max_n_imaginary,
            max_imaginary_freq=max_imaginary_freq,
            fit_method=fit_method
        )
        fc_to_db = ForceConstantsToDB(
            db_file=db_file,
            wf_uuid=wf_uuid,
            mesh_density=mesh_density,
            additional_fields=additional_fields
        )
        pass_locs = PassCalcLocs(name=name)

        tasks = [collect_structures, fit_constants, fc_to_db, pass_locs]
        super().__init__(tasks, parents=parents, name=fw_name, **kwargs)
