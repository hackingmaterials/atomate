# coding: utf-8

import os

from atomate.vasp.fireworks.core import OptimizeFW, StaticFW
from fireworks import Workflow, Firework
from atomate.vasp.powerups import (
    add_tags,
    add_additional_fields_to_taskdocs,
    add_wf_metadata,
    add_common_powerups,
)
from atomate.vasp.workflows.base.core import get_wf
from atomate.vasp.firetasks.parse_outputs import (
    MagneticDeformationToDb,
    MagneticOrderingsToDb,
)

from pymatgen.alchemy.materials import TransformedStructure

from atomate.utils.utils import get_logger

logger = get_logger(__name__)

from atomate.vasp.config import VASP_CMD, DB_FILE, ADD_WF_METADATA

from atomate.vasp.workflows.presets.scan import wf_scan_opt
from uuid import uuid4
from pymatgen.io.vasp.sets import MPRelaxSet
from pymatgen.core import Lattice, Structure
from pymatgen.analysis.magnetism.analyzer import (
    CollinearMagneticStructureAnalyzer,
    MagneticStructureEnumerator,
)

__author__ = "Matthew Horton"
__maintainer__ = "Matthew Horton"
__email__ = "mkhorton@lbl.gov"
__status__ = "Production"
__date__ = "March 2017"

__magnetic_deformation_wf_version__ = 1.2
__magnetic_ordering_wf_version__ = 2.0

module_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)))


def get_wf_magnetic_deformation(structure, c=None, vis=None):
    """
    Minimal workflow to obtain magnetic deformation proxy, as
    defined by Bocarsly et al. 2017, doi: 10.1021/acs.chemmater.6b04729

    Args:
        structure: input structure, must be structure with magnetic
    elements, such that pymatgen will initalize ferromagnetic input by
    default -- see MPRelaxSet.yaml for list of default elements
        c: Workflow config dict, in the same format
    as in presets/core.py and elsewhere in atomate
        vis: A VaspInputSet to use for the first FW

    Returns: Workflow
    """

    if not structure.is_ordered:
        raise ValueError(
            "Please obtain an ordered approximation of the input structure."
        )

    structure = structure.get_primitive_structure(use_site_props=True)

    # using a uuid for book-keeping,
    # in a similar way to other workflows
    uuid = str(uuid4())

    c_defaults = {"vasp_cmd": VASP_CMD, "db_file": DB_FILE}
    if c:
        c.update(c_defaults)
    else:
        c = c_defaults

    wf = get_wf(structure, "magnetic_deformation.yaml", common_params=c, vis=vis)

    fw_analysis = Firework(
        MagneticDeformationToDb(
            db_file=DB_FILE, wf_uuid=uuid, to_db=c.get("to_db", True)
        ),
        name="MagneticDeformationToDb",
    )

    wf.append_wf(Workflow.from_Firework(fw_analysis), wf.leaf_fw_ids)

    wf = add_common_powerups(wf, c)

    if c.get("ADD_WF_METADATA", ADD_WF_METADATA):
        wf = add_wf_metadata(wf, structure)

    wf = add_additional_fields_to_taskdocs(
        wf,
        {
            "wf_meta": {
                "wf_uuid": uuid,
                "wf_name": "magnetic_deformation",
                "wf_version": __magnetic_deformation_wf_version__,
            }
        },
    )

    return wf


class MagneticOrderingsWF:
    def __init__(
        self,
        structure,
        default_magmoms=None,
        strategies=("ferromagnetic", "antiferromagnetic"),
        automatic=True,
        truncate_by_symmetry=True,
        transformation_kwargs=None,
        static=False,
    ):
        """
        This workflow will try several different collinear
        magnetic orderings for a given input structure,
        and output a summary to a dedicated database
        collection, magnetic_orderings, in the supplied
        db_file.

        If the input structure has magnetic moments defined, it
        is possible to use these as a hint as to which elements are
        magnetic, otherwise magnetic elements will be guessed
        (this can be changed using default_magmoms kwarg).

        A brief description on how this workflow works:
            1. We make a note of the input structure, and then
               sanitize it (make it ferromagnetic, primitive)
            2. We gather information on which sites are likely
               magnetic, how many unique magnetic sites there
               are (whether one species is in several unique
               environments, e.g. tetrahedral/octahedra as Fe
               in a spinel)
            3. We generate ordered magnetic structures, first
               antiferromagnetic, and then, if appropriate,
               ferrimagnetic_Cr2NiO4 structures either by species or
               by environment -- this makes use of some new
               additions to MagOrderingTransformation to allow
               the spins of different species to be coupled together
               (e.g. all one species spin up, all other species spin
               down, with an overall order parameter of 0.5)
            4. For each ordered structure, we perform a relaxation
               and static calculation. Then an aggregation is performed
               which finds which ordering is the ground state (of
               those attempted in this specific workflow). For
               high-throughput studies, a dedicated builder is
               recommended.
            5. For book-keeping, a record is kept of whether the
               input structure is enumerated by the algorithm or
               not. This is useful when supplying the workflow with
               a magnetic structure obtained by experiment, to measure
               the performance of the workflow.

        Args:
            structure: input structure
            default_magmoms (dict): (optional, defaults provided) dict of
                magnetic elements to their initial magnetic moments in µB,
                generally these are chosen to be high-spin since they can
                relax to a low-spin configuration during a DFT electronic
                configuration
            strategies (tuple): different ordering strategies to use, choose
                from: ferromagnetic, antiferromagnetic,
                antiferromagnetic_by_motif, ferrimagnetic_by_motif and
                ferrimagnetic_by_species (here, "motif", 
                means to use a different ordering parameter for symmetry
                inequivalent sites)
            automatic (bool): if True, will automatically choose sensible
                strategies
            truncate_by_symmetry (bool): if True, will remove very
                unsymmetrical orderings that are likely physically implausible
            transformation_kwargs (dict): keyword arguments to pass to
                MagOrderingTransformation, to change automatic cell size
                limits, etc.
            static (bool): Run only static calcs (no optimization) of 
                different magnetic orderings in a fixed ground state
                geometry.

        """

        self.uuid = str(uuid4())
        self.wf_meta = {
            "wf_uuid": self.uuid,
            "wf_name": self.__class__.__name__,
            "wf_version": __magnetic_ordering_wf_version__,
        }
        self.static = static

        enumerator = MagneticStructureEnumerator(
            structure,
            default_magmoms=default_magmoms,
            strategies=strategies,
            automatic=automatic,
            truncate_by_symmetry=truncate_by_symmetry,
            transformation_kwargs=transformation_kwargs,
        )

        self.sanitized_structure = enumerator.sanitized_structure
        self.ordered_structures = enumerator.ordered_structures
        self.ordered_structure_origins = enumerator.ordered_structure_origins
        self.input_index = enumerator.input_index
        self.input_origin = enumerator.input_origin

    def get_wf(
        self, scan=False, perform_bader=True, num_orderings_hard_limit=16, c=None
    ):
        """
        Retrieve the FireWorks workflow.

        Args:
            scan (bool): if True, use the SCAN functional instead of GGA+U,
                since the SCAN functional has shown to have improved
                performance for magnetic systems in some cases
            perform_bader (bool): if True, make sure the "bader" binary is in
                your path, will use Bader analysis to calculate
                atom-projected magnetic moments
            num_orderings_hard_limit (int): will make sure total number of
                magnetic orderings does not exceed this number even if there
                are extra orderings of equivalent symmetry
            c (dict): additional config dict (as used elsewhere in atomate)

        Returns: FireWorks Workflow

        """

        c_defaults = {"VASP_CMD": VASP_CMD, "DB_FILE": DB_FILE}
        additional_fields = {"relax": not self.static}
        c = c or {}
        for k, v in c_defaults.items():
            if k not in c:
                c[k] = v

        fws = []
        analysis_parents = []

        # trim total number of orderings (useful in high-throughput context)
        # this is somewhat course, better to reduce num_orderings kwarg and/or
        # change enumeration strategies
        ordered_structures = self.ordered_structures
        ordered_structure_origins = self.ordered_structure_origins

        def _add_metadata(structure):
            """
            For book-keeping, store useful metadata with the Structure
            object for later database ingestion including workflow
            version and a UUID for easier querying of all tasks generated
            from the workflow.

            Args:
                structure: Structure

            Returns: TransformedStructure

            """

            # this could be further improved by storing full transformation
            # history, but would require an improved transformation pipeline
            return TransformedStructure(
                structure, other_parameters={"wf_meta": self.wf_meta}
            )

        ordered_structures = [_add_metadata(struct) for struct in ordered_structures]

        if (
            num_orderings_hard_limit
            and len(self.ordered_structures) > num_orderings_hard_limit
        ):
            ordered_structures = self.ordered_structures[0:num_orderings_hard_limit]
            ordered_structure_origins = self.ordered_structure_origins[
                0:num_orderings_hard_limit
            ]
            logger.warning(
                "Number of ordered structures exceeds hard limit, "
                "removing last {} structures.".format(
                    len(self.ordered_structures) - len(ordered_structures)
                )
            )
            # always make sure input structure is included
            if self.input_index and self.input_index > num_orderings_hard_limit:
                ordered_structures.append(self.ordered_structures[self.input_index])
                ordered_structure_origins.append(
                    self.ordered_structure_origins[self.input_index]
                )

        # default incar settings
        user_incar_settings = {"ISYM": 0, "LASPH": True, "EDIFFG": -0.05}
        if scan:
            # currently, using SCAN relaxation as a static calculation also
            # since it is typically high quality enough, but want to make
            # sure we are also writing the AECCAR* files
            user_incar_settings.update({"LAECHG": True})
        user_incar_settings.update(c.get("user_incar_settings", {}))
        c["user_incar_settings"] = user_incar_settings

        for idx, ordered_structure in enumerate(ordered_structures):

            analyzer = CollinearMagneticStructureAnalyzer(ordered_structure)

            name = " ordering {} {} -".format(idx, analyzer.ordering.value)

            if not scan:

                vis = MPRelaxSet(
                    ordered_structure, user_incar_settings=user_incar_settings
                )

                if not self.static:

                    # relax
                    fws.append(
                        OptimizeFW(
                            ordered_structure,
                            vasp_input_set=vis,
                            vasp_cmd=c["VASP_CMD"],
                            db_file=c["DB_FILE"],
                            max_force_threshold=0.05,
                            half_kpts_first_relax=False,
                            name=name + " optimize",
                        )
                    )

                # static
                fws.append(
                    StaticFW(
                        ordered_structure,
                        vasp_cmd=c["VASP_CMD"],
                        db_file=c["DB_FILE"],
                        name=name + " static",
                        prev_calc_loc=True,
                        parents=fws[-1],
                        vasptodb_kwargs={'parse_chgcar': True, 'parse_aeccar': True}
                    )
                )
                
                if not self.static:
                    # so a failed optimize doesn't crash workflow
                    fws[-1].spec["_allow_fizzled_parents"] = True

            elif scan:

                # wf_scan_opt is just a single FireWork so can append it directly
                scan_fws = wf_scan_opt(ordered_structure, c=c).fws
                # change name for consistency with non-SCAN
                new_name = scan_fws[0].name.replace(
                    "structure optimization", name + " optimize"
                )
                scan_fws[0].name = new_name
                scan_fws[0].tasks[-1]["additional_fields"]["task_label"] = new_name
                fws += scan_fws

            analysis_parents.append(fws[-1])

        fw_analysis = Firework(
            MagneticOrderingsToDb(
                db_file=c["DB_FILE"],
                wf_uuid=self.uuid,
                parent_structure=self.sanitized_structure,
                origins=ordered_structure_origins,
                input_index=self.input_index,
                perform_bader=perform_bader,
                scan=scan,
                additional_fields=additional_fields,
            ),
            name="Magnetic Orderings Analysis",
            parents=analysis_parents,
            spec={"_allow_fizzled_parents": True},
        )
        fws.append(fw_analysis)

        formula = self.sanitized_structure.composition.reduced_formula
        wf_name = "{} - magnetic orderings".format(formula)
        if scan:
            wf_name += " - SCAN"
        wf = Workflow(fws, name=wf_name)

        wf = add_additional_fields_to_taskdocs(wf, {"wf_meta": self.wf_meta})

        tag = "magnetic_orderings group: >>{}<<".format(self.uuid)
        wf = add_tags(wf, [tag, ordered_structure_origins])

        return wf


def get_commensurate_orderings(magnetic_structures, energies):
    """Generate supercells for static calculations.

    From the ground state magnetic ordering and first excited state,
    generate supercells with magnetic orderings.

    If the gs is FM, match supercells to the 1st es.
    If the gs is AFM, FM is included by default, 1st es may not be.
    If the gs is FiM, 1st es may not be captured.

    Args:
        magnetic_structures (list): Structures.
        energies (list): Energies per atom.

    Returns:
        matched_structures (list): Commensurate supercells for static 
            calculations.

    TODO:
        * Only consider orderings with |S_i| = ground state
        * Constrain noncollinear magmoms

    """

    # Sort by energies
    ordered_structures = [
        s for _, s in sorted(zip(energies, magnetic_structures), reverse=False)
    ]
    ordered_energies = sorted(energies, reverse=False)

    # Ground state and 1st excited state
    gs_struct = ordered_structures[0]
    es_struct = ordered_structures[1]

    cmsa = CollinearMagneticStructureAnalyzer(
        es_struct, threshold=0.0, threshold_nonmag=1.0, make_primitive=False
    )
    es_ordering = cmsa.ordering.value
    es_struct = cmsa.structure

    cmsa = CollinearMagneticStructureAnalyzer(
        gs_struct, threshold=0.0, threshold_nonmag=1.0, make_primitive=False
    )
    gs_ordering = cmsa.ordering.value
    gs_struct = cmsa.structure

    # FM gs will always be commensurate so we match to the 1st es
    if gs_ordering == "FM":
        enum_struct = es_struct
        fm_moments = np.array(gs_struct.site_properties["magmom"])
        fm_moment = np.mean(fm_moments[np.nonzero(fm_moments)])
        gs_magmoms = [fm_moment for m in es_struct.site_properties["magmom"]]
    elif gs_ordering in ["FiM", "AFM"]:
        enum_struct = gs_struct
        gs_magmoms = [abs(m) for m in gs_struct.site_properties["magmom"]]

    mse = MagneticStructureEnumerator(
        enum_struct,
        strategies=("ferromagnetic", "antiferromagnetic"),
        automatic=False,
        transformation_kwargs={
            "min_cell_size": 1,
            "max_cell_size": 2,
            "check_ordered_symmetry": False,
        },
    )

    # Enumerator bookkeeping
    input_index = mse.input_index
    ordered_structure_origins = mse.ordered_structure_origins

    matched_structures = []

    sm = StructureMatcher(
        primitive_cell=False, attempt_supercell=True, comparator=ElementComparator()
    )

    # Get commensurate supercells
    for s in mse.ordered_structures:
        try:
            s2 = sm.get_s2_like_s1(enum_struct, s)
        except:
            s2 = None
        if s2 is not None:
            # Standardize magnetic structure
            cmsa = CollinearMagneticStructureAnalyzer(
                s2, threshold=0.0, make_primitive=False
            )
            s2 = cmsa.structure
            matched_structures.append(s2)

    # Find the gs ordering in the enumerated supercells
    cmsa = CollinearMagneticStructureAnalyzer(
        gs_struct, threshold=0.0, make_primitive=False
    )

    enum_index = [cmsa.matches_ordering(s) for s in matched_structures].index(True)

    # Enforce all magmom magnitudes to match the gs
    for s in matched_structures:
        ms = s.site_properties["magmom"]
        magmoms = [np.sign(m1) * m2 for m1, m2 in zip(ms, gs_magmoms)]
        s.add_site_property("magmom", magmoms)

    return matched_structures, input_index, ordered_structure_origins


if __name__ == "__main__":

    # for trying workflows

    from fireworks import LaunchPad

    latt = Lattice.cubic(4.17)
    species = ["Ni", "O"]
    coords = [[0.00000, 0.00000, 0.00000], [0.50000, 0.50000, 0.50000]]
    NiO = Structure.from_spacegroup(225, latt, species, coords)

    wf_deformation = get_wf_magnetic_deformation(NiO)

    wf_orderings = MagneticOrderingsWF(NiO).get_wf()

    lpad = LaunchPad.auto_load()
    lpad.add_wf(wf_orderings)
    lpad.add_wf(wf_deformation)
