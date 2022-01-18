import itertools
import os
from uuid import uuid4

import numpy as np
from fireworks import Firework, Workflow
from pymatgen.core import Element, Structure
from pymatgen.io.vasp import Incar, Kpoints, Poscar
from pymatgen.io.vasp.sets import MPStaticSet

from atomate.common.firetasks.glue_tasks import PassCalcLocs
from atomate.utils.utils import get_logger
from atomate.vasp.config import ADD_WF_METADATA, DB_FILE, VASP_CMD
from atomate.vasp.firetasks.glue_tasks import CopyVaspOutputs
from atomate.vasp.firetasks.parse_outputs import HubbardHundLinRespToDb, VaspToDb
from atomate.vasp.firetasks.run_calc import RunVaspCustodian
from atomate.vasp.firetasks.write_inputs import WriteVaspFromIOSet
from atomate.vasp.powerups import (
    add_additional_fields_to_taskdocs,
    add_common_powerups,
    add_wf_metadata,
)

__author__ = "Guy Moore, Martin Siron"
__maintainer__ = "Guy Moore"
__email__ = "gmoore@lbl.gov"
__status__ = "Production"
__date__ = "February 2020"
logger = get_logger(__name__)

__hubbard_hund_linresp_wf_version__ = 0.0

module_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)))


def get_wf_hubbard_hund_linresp(
    structure,
    user_incar_settings=None,
    relax_nonmagnetic=True,
    spin_polarized=True,
    applied_potential_range=(-0.2, 0.2),
    num_evals=9,
    site_indices_perturb=None,
    species_perturb=None,
    find_nearest_sites=True,
    parallel_scheme=0,
    ediff_tight=None,
    c=None,
):
    """
    Compute Hubbard U (and Hund J) on-site interaction values using GGA+U
    linear response method proposed by Cococcioni et. al.
    (DOI: 10.1103/PhysRevB.71.035105)
    and the spin-polarized response formalism developed by Linscott et. al.
    (DOI: 10.1103/PhysRevB.98.235157)

    This workflow relies on the constrained on-site potential functional implemented in VASP,
    with a helpful tutorial found here:
    https://www.vasp.at/wiki/index.php/Calculate_U_for_LSDA%2BU

    Args:
        structure:
        user_incar_settings: user INCAR settings
        relax_nonmagnetic: Restart magnetic SCF runs from
    non-magnetic calculation, using WAVECAR
        spin_polarized: Perform spin-dependent perturbations
        applied_potential_range: Bounds of applied potential
        num_evals: Number of perturbation evaluations
        site_indices_perturb: (must specify if species_perturb=None)
    List of site indices within
    Structure indicating perturbation sites;
        species_perturb: (must specify if site_indices_perturb=None)
    List of names of species (string)
    of sites to perturb; First site of that species
    is selected in the structure
        find_nearest_sites: If set to true and species_perturb != None,
    the closest sites (by the Structure distance matrix) will be selected
    in the response analysis to account for inter-site screening effects
        parallel_scheme: 0 - (default) self-consistent (SCF)
    runs use WAVECAR from non-self consistent (NSCF) run
    at same applied potential; 1 - SCF runs use WAVECAR
    from ground-state (V=0) run.
    While reusing the WAVECAR from NSCF run in SCF run may be more
    efficient (parallel_scheme: 0), the user may also choose to
    remove the dependency between NSCF and SCF runs
    (parallel_scheme: 1)
        ediff_tight: Final energy convergence tolerance,
    if restarting from a previous run
    (if not specified, will default to pymatgen default EDIFF)
        c: Workflow config dict, in the same format
    as in presets/core.py and elsewhere in atomate

    Returns: Workflow
    """

    if not structure.is_ordered:
        raise ValueError(
            "Please obtain an ordered approximation of the input structure."
        )

    if not site_indices_perturb:
        site_indices_perturb = []

    if species_perturb:

        if find_nearest_sites:
            site_indices_perturb = find_closest_sites(structure, species_perturb)
        else:
            for specie_u in species_perturb:
                found_specie = False
                for s in range(len(structure)):
                    site = structure[s]
                    if (Element(str(site.specie)) == Element(specie_u)) and (
                        s not in site_indices_perturb
                    ):
                        found_specie = True
                        break
                if not found_specie:
                    raise ValueError("Could not find specie(s) in structure.")
                site_indices_perturb.append(s)

    elif not site_indices_perturb:
        logger.warning(
            "Sites for computing U value are not specified. "
            "Computing U for first site in structure. "
        )

    site_indices_perturb = list(tuple(site_indices_perturb))
    num_perturb = len(site_indices_perturb)

    sites_perturb = []
    for site_index_perturb in site_indices_perturb:
        site = structure[site_index_perturb]
        sites_perturb.append(site)

    structure.remove_sites(indices=site_indices_perturb)

    for site in sites_perturb:
        structure.insert(
            i=0,
            species=site.specie,
            coords=site.frac_coords,
            properties=site.properties,
        )

    # using a uuid for book-keeping,
    # in a similar way to other workflows
    uuid = str(uuid4())

    c_defaults = {"vasp_cmd": VASP_CMD, "db_file": DB_FILE}
    if c:
        c.update(c_defaults)
    else:
        c = c_defaults

    # Calculate groundstate

    # set user_incar_settings
    if not user_incar_settings:
        user_incar_settings = {}

    # setup VASP input sets
    uis_gs, uis_ldau, val_dict, vis_ldau = init_linresp_input_sets(
        user_incar_settings, structure, num_perturb
    )

    fws = []
    index_fw_gs = [0]

    ediff_default = vis_ldau.incar["EDIFF"]
    if not ediff_tight:
        ediff_tight = 0.1 * ediff_default

    append_linresp_ground_state_fws(
        fws,
        structure,
        num_perturb,
        index_fw_gs,
        uis_gs,
        relax_nonmagnetic,
        ediff_default,
        ediff_tight,
    )

    # generate list of applied on-site potentials in linear response
    applied_potential_value_list = []
    for counter_perturb in range(num_perturb):
        applied_potential_values = np.linspace(
            applied_potential_range[0], applied_potential_range[1], num_evals
        )
        applied_potential_values = np.around(applied_potential_values, decimals=9)

        if 0.0 in applied_potential_values:
            applied_potential_values = list(applied_potential_values)
            applied_potential_values.pop(applied_potential_values.index(0.0))
            applied_potential_values = np.array(applied_potential_values)

        applied_potential_value_list.append(applied_potential_values.copy())

    for counter_perturb in range(num_perturb):

        applied_potential_values = applied_potential_value_list[counter_perturb]

        for v in applied_potential_values:

            append_linresp_perturb_fws(
                v,
                fws,
                structure,
                counter_perturb,
                num_perturb,
                index_fw_gs,
                uis_ldau,
                val_dict,
                spin_polarized,
                relax_nonmagnetic,
                ediff_default,
                ediff_tight,
                parallel_scheme,
            )

    wf = Workflow(fws)

    fw_analysis = Firework(
        HubbardHundLinRespToDb(
            num_perturb=num_perturb,
            spin_polarized=spin_polarized,
            relax_nonmagnetic=relax_nonmagnetic,
            db_file=DB_FILE,
            wf_uuid=uuid,
        ),
        name="HubbardHundLinRespToDb",
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
                "wf_name": "hubbard_hund_linresp",
                "wf_version": __hubbard_hund_linresp_wf_version__,
            }
        },
    )

    return wf


def find_closest_sites(struct, species_perturb):
    """
    Function to find closest cluster of sites of target species of
    sites to perturb
    """

    d_ij = struct.distance_matrix

    n_species = {}
    for site in struct:
        k = str(site.specie)
        if k not in n_species.keys():
            n_species.update({k: 0})
        n_species[k] += 1

    n_s = len(n_species.keys())
    n_config = 1
    for k in n_species:
        n_config *= n_species[k]
    indices = -1 * np.ones(n_s, int)
    dist_min = -1.0
    indices_nearest = indices.copy()
    for k in range(n_config):
        denom = 1
        for j in range(n_s):
            n = n_species[list(n_species.keys())[j]]
            indices[j] = np.mod(k // denom, n)
            denom *= n
        nn_s = [n_species[list(n_species.keys())[j]] for j in range(n_s)]
        iindxs = [indices[j] + int(np.sum(nn_s[0:j])) for j in range(n_s)]
        indxs = []
        for j, jk in enumerate(n_species):
            if jk in species_perturb:
                indxs.append(iindxs[j])
        dist = 0.0
        for x in indxs:
            for y in indxs:
                dist += d_ij[x, y]
        if dist < dist_min or dist_min == -1.0:
            dist_min = dist
            indices_nearest = indxs.copy()

    return indices_nearest


def init_linresp_input_sets(
    user_incar_settings,
    structure,
    num_perturb,
):
    """
    Function to setup VASP input sets for ground-state
    and perturbation Fireworks
    """

    # set LMAXMIX in user_incar_settings
    if "LMAXMIX" not in user_incar_settings.keys():
        lmaxmix_dict = {"p": 2, "d": 4, "f": 6}
        lmm = "p"
        for site in structure:
            block = str(site.specie.block)
            if block == "d" and lmm != "f":
                lmm = "d"
            elif block == "f":
                lmm = "f"
                break
        user_incar_settings.update({"LMAXMIX": lmaxmix_dict[lmm]})

    # ground state user incar settings
    uis_gs = user_incar_settings.copy()
    if "LMAXMIX" not in uis_gs:
        logger.warning("You have not specified LMAXMIX. Defaulting to VASP default.")
    uis_gs.update({"LDAU": False, "LORBIT": 11, "LWAVE": True, "ISPIN": 2})

    uis_ldau = uis_gs.copy()
    uis_ldau.update({"LDAU": True, "LDAUTYPE": 3, "LDAUPRINT": 2})

    # Initialize vasp input set
    vis_params = {"user_incar_settings": uis_ldau.copy()}
    vis_ldau = HubbardHundLinRespSet(
        structure=structure, num_perturb=num_perturb, **vis_params.copy()
    )

    val_dict = {"LDAUL": {}, "LDAUU": {}, "LDAUJ": {}}
    for k in ["LDAUL", "LDAUU", "LDAUJ"]:
        if k == "LDAUL":
            # for LDAUL
            for i in range(num_perturb):
                val_dict[k].update({f"perturb{i}": -1})
            for s in vis_ldau.poscar.site_symbols[num_perturb:]:
                val_dict[k].update({s: -1})
        else:
            # for LDAUU and LDAUJ
            for i in range(num_perturb):
                val_dict[k].update({f"perturb{i}": 0})
            for s in vis_ldau.poscar.site_symbols[num_perturb:]:
                v = 0
                val_dict[k].update({s: v})
        uis_ldau.update({k: val_dict[k].copy()})
    vis_params = {"user_incar_settings": uis_ldau.copy()}
    vis_ldau = HubbardHundLinRespSet(
        structure=structure, num_perturb=num_perturb, **vis_params.copy()
    )

    uis_gs = uis_ldau.copy()

    return uis_gs, uis_ldau, val_dict, vis_ldau


def append_linresp_ground_state_fws(
    fws,
    structure,
    num_perturb,
    index_fw_gs,
    uis_gs,
    relax_nonmagnetic,
    ediff_default,
    ediff_tight,
):
    """
    Appends the ground-state (zero potential) Fireworks to the list of Fireworks
    """

    if relax_nonmagnetic:
        uis_gs.update({"ISPIN": 1, "EDIFF": ediff_default})
    else:
        uis_gs.update({"ISPIN": 2, "EDIFF": ediff_tight})
    vis_params = {"user_incar_settings": uis_gs.copy()}
    vis_gs = HubbardHundLinRespSet(
        structure=structure, num_perturb=num_perturb, **vis_params.copy()
    )
    fws.append(
        HubbardHundLinRespFW(
            structure=structure,
            name="initial_static",
            vasp_input_set=vis_gs,
            vasp_cmd=VASP_CMD,
            db_file=DB_FILE,
        )
    )

    if relax_nonmagnetic:
        index_fw_gs[0] += 1
        uis_gs.update({"ISPIN": 2, "ICHARG": 1, "EDIFF": ediff_tight})
        additional_files = ["WAVECAR", "CHGCAR"]
        vis_params = {"user_incar_settings": uis_gs.copy()}
        vis_gs = HubbardHundLinRespSet(
            structure=structure, num_perturb=num_perturb, **vis_params.copy()
        )
        fws.append(
            HubbardHundLinRespFW(
                structure=structure,
                parents=fws[-1],
                additional_files=additional_files.copy(),
                name="initial_static_magnetic",
                vasp_input_set=vis_gs,
                vasp_cmd=VASP_CMD,
                db_file=DB_FILE,
            )
        )


def append_linresp_perturb_fws(
    applied_pot,
    fws,
    structure,
    counter_perturb,
    num_perturb,
    index_fw_gs,
    uis_ldau,
    val_dict,
    spin_polarized,
    relax_nonmagnetic,
    ediff_default,
    ediff_tight,
    parallel_scheme,
):
    """
    Appends the perturbed (non-zero on-site potential) Fireworks to the list of Fireworks
    """

    v = applied_pot

    sign = "neg" if v < 0 else "pos"
    signs = []

    spin_potential_values = []
    if spin_polarized:
        spin_potential_values = [{"LDAUU": v, "LDAUJ": 0.0}, {"LDAUU": 0.0, "LDAUJ": v}]
        signs = [{"LDAUU": sign, "LDAUJ": ""}, {"LDAUU": "", "LDAUJ": sign}]
    else:
        spin_potential_values = [{"LDAUU": v, "LDAUJ": v}]
        signs = [{"LDAUU": sign, "LDAUJ": sign}]

    block_dict = {"s": 0, "p": 1, "d": 2, "f": 3}

    for spin_pot_dict, sign_dict in zip(spin_potential_values, signs):

        # Update perturbation potential for U and J
        for k in ["LDAUL", "LDAUU", "LDAUJ"]:
            if k == "LDAUL":
                # for LDAUL
                for i in range(num_perturb):
                    if i == counter_perturb:
                        block = str(structure[counter_perturb].specie.block)
                        val_dict[k].update({f"perturb{i}": block_dict[block]})
                    else:
                        val_dict[k].update({f"perturb{i}": -1})
            else:
                # for LDAUU and LDAUJ
                for i in range(num_perturb):
                    if i == counter_perturb:
                        val_dict[k].update({f"perturb{i}": spin_pot_dict[k]})
                    else:
                        val_dict[k].update({f"perturb{i}": 0})
            uis_ldau.update({k: val_dict[k].copy()})

        # Non-SCF runs
        uis_ldau.update({"ISTART": 1, "ICHARG": 11, "ISPIN": 2, "EDIFF": ediff_default})

        vis_params = {"user_incar_settings": uis_ldau.copy()}
        vis_ldau = HubbardHundLinRespSet(
            structure=structure, num_perturb=num_perturb, **vis_params.copy()
        )

        parents = fws[index_fw_gs[0]]

        additional_files = ["WAVECAR", "CHGCAR"]

        fw = HubbardHundLinRespFW(
            structure=structure,
            parents=parents,
            name="nscf_site{}_vup_{}{}_vdn_{}{}".format(
                counter_perturb,
                sign_dict["LDAUU"],
                abs(round(spin_pot_dict["LDAUU"], 6)),
                sign_dict["LDAUJ"],
                abs(round(spin_pot_dict["LDAUJ"], 6)),
            ),
            vasp_input_set=vis_ldau,
            additional_files=additional_files.copy(),
            vasp_cmd=VASP_CMD,
            db_file=DB_FILE,
        )
        fws.append(fw)

        # SCF runs
        uis_ldau.update({"ISTART": 1, "ICHARG": 0})
        if relax_nonmagnetic:
            uis_ldau.update({"ISPIN": 1, "EDIFF": ediff_default})
        else:
            uis_ldau.update({"ISPIN": 2, "EDIFF": ediff_tight})

        vis_params = {"user_incar_settings": uis_ldau.copy()}
        vis_ldau = HubbardHundLinRespSet(
            structure=structure, num_perturb=num_perturb, **vis_params.copy()
        )

        if parallel_scheme == 0:
            parents = fws[-1]
        else:
            parents = fws[index_fw_gs[0]]

        additional_files = ["WAVECAR"]

        fw = HubbardHundLinRespFW(
            structure=structure,
            parents=parents,
            name="scf_site{}_vup_{}{}_vdn_{}{}".format(
                counter_perturb,
                sign_dict["LDAUU"],
                abs(round(spin_pot_dict["LDAUU"], 6)),
                sign_dict["LDAUJ"],
                abs(round(spin_pot_dict["LDAUJ"], 6)),
            ),
            vasp_input_set=vis_ldau,
            additional_files=additional_files.copy(),
            vasp_cmd=VASP_CMD,
            db_file=DB_FILE,
        )
        fws.append(fw)

        # SCF magnetic runs
        if relax_nonmagnetic:
            uis_ldau.update(
                {"ISTART": 1, "ICHARG": 1, "ISPIN": 2, "EDIFF": ediff_tight}
            )

            vis_params = {"user_incar_settings": uis_ldau.copy()}
            vis_ldau = HubbardHundLinRespSet(
                structure=structure, num_perturb=num_perturb, **vis_params.copy()
            )

            parents = fws[-1]
            additional_files = ["WAVECAR", "CHGCAR"]

            fw = HubbardHundLinRespFW(
                structure=structure,
                parents=parents,
                name="scf_magnetic_site{}_vup_{}{}_vdn_{}{}".format(
                    counter_perturb,
                    sign_dict["LDAUU"],
                    abs(round(spin_pot_dict["LDAUU"], 6)),
                    sign_dict["LDAUJ"],
                    abs(round(spin_pot_dict["LDAUJ"], 6)),
                ),
                vasp_input_set=vis_ldau,
                additional_files=additional_files.copy(),
                vasp_cmd=VASP_CMD,
                db_file=DB_FILE,
            )
            fws.append(fw)


class PoscarPerturb(Poscar):
    """
    Derived Poscar class that allows the distinction of individual sites
    in the Structure
    """

    def __init__(self, structure: Structure, num_perturb: int = 1, **kwargs):
        """
        Args:
            structure:
            num_perturb: Number of sites to perturb;
        First n sites are indicated as "separate" species
            **kwargs:
        """
        super().__init__(structure, sort_structure=False, **kwargs)

        self.structure = structure
        self.num_perturb = num_perturb

    @property
    def site_symbols(self):
        """
        Sequence of symbols associated with the Poscar
        """

        if self.num_perturb > 0 and self.num_perturb <= len(self.structure):
            syms = [site.specie.symbol for site in self.structure[self.num_perturb :]]
            syms = [a[0] for a in itertools.groupby(syms)]
            syms_perturb = [
                site.specie.symbol for site in self.structure[0 : self.num_perturb]
            ]
            syms_perturb.extend(syms)
        else:
            raise ValueError("Invalid atom index to perturb")

        return syms_perturb

    @property
    def natoms(self):
        """
        Sequence of number of sites of each type associated with the Poscar
        """

        if self.num_perturb > 0 and self.num_perturb <= len(self.structure):
            syms = [site.specie.symbol for site in self.structure[self.num_perturb :]]
            n_atoms = [len(tuple(a[1])) for a in itertools.groupby(syms)]

            n_atoms_perturb = [1 for i in range(self.num_perturb)]
            n_atoms_perturb.extend(n_atoms)
        else:
            raise ValueError("Invalid atom index to perturb")

        return n_atoms_perturb


class HubbardHundLinRespSet(MPStaticSet):
    """
    VASP input set for Linear Response Hubbard U workflow perturbation Fireworks
    """

    def __init__(
        self,
        structure,
        num_perturb,
        prev_incar=None,
        prev_kpoints=None,
        reciprocal_density=100,
        small_gap_multiply=None,
        **kwargs,
    ):
        """
        Args:
            structure:
            num_perturb: Number of sites to perturb
            prev_incar:
            prev_kpoints:
            reciprocal_density:
            small_gap_multiply:
            **kwargs:
        """

        super().__init__(structure, sort_structure=False, **kwargs)

        self.num_perturb = num_perturb

        if isinstance(prev_kpoints, str):
            prev_kpoints = Kpoints.from_file(prev_kpoints)
        self.prev_kpoints = prev_kpoints

        self.reciprocal_density = reciprocal_density
        self.kwargs = kwargs
        self.small_gap_multiply = small_gap_multiply

    @property
    def incar(self):
        """
        Custom Incar attribute for HubbardHundLinRespSet
        """
        parent_incar = super().incar
        incar = Incar(parent_incar)

        incar.update({"ISYM": -1, "ISMEAR": 0, "LREAL": False, "LASPH": True})
        incar.update({"ISTART": 1})
        # incar.update({"ALGO": "Fast"})
        incar.pop("NSW", None)

        if self.kwargs.get("user_incar_settings")["LDAUU"]:

            incar.update({"LDAUL": self.kwargs.get("user_incar_settings")["LDAUL"]})
            incar.update({"LDAUU": self.kwargs.get("user_incar_settings")["LDAUU"]})
            incar.update({"LDAUJ": self.kwargs.get("user_incar_settings")["LDAUJ"]})

            incar["LDAUL"] = [incar["LDAUL"][key] for key in incar["LDAUL"].keys()]
            incar["LDAUU"] = [incar["LDAUU"][key] for key in incar["LDAUU"].keys()]
            incar["LDAUJ"] = [incar["LDAUJ"][key] for key in incar["LDAUJ"].keys()]

            incar["LDAU"] = self.kwargs.get("user_incar_settings")["LDAU"]
            incar["LDAUTYPE"] = self.kwargs.get("user_incar_settings")["LDAUTYPE"]
            incar["LDAUPRINT"] = self.kwargs.get("user_incar_settings")["LDAUPRINT"]
            incar["LORBIT"] = self.kwargs.get("user_incar_settings")["LORBIT"]

        return incar

    @property
    def poscar(self):
        """
        Custom Poscar for HubbardHundLinRespSet
        """
        poscar = PoscarPerturb(
            structure=super().structure, num_perturb=self.num_perturb
        )
        return poscar

    @property
    def kpoints(self):
        """
        Custom Kpoints for HubbardHundLinRespSet
        """
        kpoints = super().kpoints

        return kpoints


class HubbardHundLinRespFW(Firework):
    def __init__(
        self,
        structure=None,
        name="linresponse_U",
        vasp_input_set=None,
        vasp_input_set_params=None,
        vasp_cmd=VASP_CMD,
        prev_calc_loc=True,
        prev_calc_dir=None,
        db_file=DB_FILE,
        vasptodb_kwargs=None,
        parents=None,
        additional_files=None,
        **kwargs,
    ):
        """
        Standard static calculation Firework - either from a previous location
        or from a structure.

        Args:
            structure (Structure): Input structure. Note that for prev_calc_loc
                jobs, the structure is only used to set the name of the FW and
                any structure with the same composition can be used.
            name (str): Name for the Firework.
            vasp_input_set (VaspInputSet): input set to use (for jobs w/no parents)
                Defaults to MPStaticSet() if None.
            vasp_input_set_params (dict): Dict of vasp_input_set kwargs.
            vasp_cmd (str): Command to run vasp.
            prev_calc_loc (bool or str): If true (default), copies outputs from
                previous calc. If str value, retrieves a previous calculation
                output by name. If False/None, will create
                new static calculation using the provided structure.
            prev_calc_dir (str): Path to a previous calculation to copy from
            db_file (str): Path to file specifying db credentials.
            parents (Firework): Parents of this particular Firework. FW or list
                of FWS.
            vasptodb_kwargs (dict): kwargs to pass to VaspToDb
            kwargs: Other kwargs that are passed to Firework.__init__().
        """
        t = []

        vasp_input_set_params = vasp_input_set_params or {}
        vasptodb_kwargs = vasptodb_kwargs or {}
        if "additional_fields" not in vasptodb_kwargs:
            vasptodb_kwargs["additional_fields"] = {}
        vasptodb_kwargs["additional_fields"]["task_label"] = name

        formula = structure.composition.reduced_formula if structure else "unknown"
        fw_name = f"{formula}-{name}"

        if prev_calc_dir:
            t.append(
                CopyVaspOutputs(
                    calc_dir=prev_calc_dir,
                    additional_files=additional_files,
                    contcar_to_poscar=False,
                )
            )
        elif parents:
            if prev_calc_loc:
                t.append(
                    CopyVaspOutputs(
                        calc_loc=prev_calc_loc,
                        additional_files=additional_files,
                        contcar_to_poscar=False,
                    )
                )

        if structure:
            vasp_input_set = vasp_input_set or HubbardHundLinRespSet(
                structure, **vasp_input_set_params
            )
            t.append(
                WriteVaspFromIOSet(structure=structure, vasp_input_set=vasp_input_set)
            )
        else:
            raise ValueError("Must specify structure")

        t.append(RunVaspCustodian(vasp_cmd=vasp_cmd, auto_npar=">>auto_npar<<"))
        t.append(PassCalcLocs(name=name))
        t.append(VaspToDb(db_file=db_file, **vasptodb_kwargs))
        super().__init__(t, parents=parents, name=fw_name, **kwargs)
