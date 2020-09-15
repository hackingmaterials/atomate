# coding: utf-8

import os
import itertools

from atomate.vasp.fireworks.core import StaticFW
from fireworks import Workflow, Firework
from atomate.vasp.powerups import (
    add_tags,
    add_additional_fields_to_taskdocs,
    add_wf_metadata,
    add_common_powerups,
)
from atomate.vasp.workflows.base.core import get_wf

from atomate.vasp.firetasks.run_calc import RunVaspCustodian
from atomate.common.firetasks.glue_tasks import PassCalcLocs
from atomate.vasp.firetasks.glue_tasks import CopyVaspOutputs
from atomate.vasp.firetasks.parse_outputs import VaspToDb, HubbardHundLinRespToDb
from atomate.vasp.firetasks.write_inputs import WriteVaspFromIOSet

from atomate.utils.utils import get_logger

logger = get_logger(__name__)

from atomate.vasp.config import VASP_CMD, DB_FILE, ADD_WF_METADATA

from atomate.vasp.workflows.presets.scan import wf_scan_opt
from uuid import uuid4
from pymatgen.io.vasp.sets import MPRelaxSet, MPStaticSet
from pymatgen.io.vasp.inputs import Poscar, Incar
from pymatgen.core import Lattice, Structure, Element

import numpy as np

__author__ = "Guy Moore, Martin Siron"
__maintainer__ = "Guy Moore"
__email__ = "gmoore@lbl.gov"
__status__ = "Production"
__date__ = "February 2020"

__hubbard_hund_linresp_wf_version__ = 0.0

module_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)))

def get_wf_hubbard_hund_linresp(structure, uis,
                             relax_nonmagnetic=True,
                             spin_polarized=True,
                             applied_potential_range=[-0.2, 0.2],
                             num_parallel_evals=9, num_total_evals=9,
                             site_indices_u=None, species_u=None,
                             find_nearest_sites=True,
                             parallel_scheme=1,
                             ediff_tight=1.0e-6,
                             c=None):
    """
    Compute Hubbard U on-site interaction values using GGA+U linear response method 
    proposed by Cococcioni et. al. (DOI: 10.1103/PhysRevB.71.035105).

    Args:
        structure:
        spin_polarized: Perform spin-dependent perturbations
        applied_potential_range: Bounds of applied potential 
        num_evals: Number of perturbation evalutaions
        uis: user INCAR settings (must set LMAXMIX, i.e. uis = {"LMAXMIX": 4}) 
        site_indices_u: List of site indices within 
    Structure indicating perturbation sites
        species_u: List of names of species (string) 
    of sites to perturb; First site of that species 
    is selected in the structure
        use_default_uvals: Use the default U values 
    for non-perturbed sites
        parallel_scheme: 0 - (default) self-consistent (SCF) 
    runs use WAVECAR from non-self consistent (NSCF) run
    at same applied potential; 1 - SCF runs use WAVECAR 
    from ground-state (V=0) run
        c: Workflow config dict, in the same format
    as in presets/core.py and elsewhere in atomate

    Returns: Workflow
    """

    if not structure.is_ordered:
        raise ValueError(
            "Please obtain an ordered approximation of the input structure."
        )

    def find_closest_sites(struct, species_u):

        Dij = struct.distance_matrix

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
        indices = -1*np.ones(n_s, int)
        dist_min = -1.0
        indices_nearest = indices.copy()
        for k in range(n_config):
            denom = 1
            for j in range(n_s):
                n = n_species[list(n_species.keys())[j]]
                indices[j] = np.mod(k // denom, n)
                denom *= n
            Ns = [n_species[list(n_species.keys())[j]] for j in range(n_s)]
            iindxs = [indices[j]+int(np.sum(Ns[0:j])) for j in range(n_s)]            
            indxs = []
            for j,jk in enumerate(n_species):
                if jk in species_u:
                    indxs.append(iindxs[j])
            dist = 0.0
            for x in indxs:
                for y in indxs:
                    dist += Dij[x,y]
            if dist < dist_min or dist_min == -1.0:
                dist_min = dist
                indices_nearest = indxs.copy()

        return indices_nearest

    if not site_indices_u:
        site_indices_u = []

    if species_u:

        if find_nearest_sites:
            site_indices_u = find_closest_sites(structure, species_u)
        else:
            for specie_u in species_u:
                foundSpecie = False
                for s in range(len(structure)):
                    site = structure[s]
                    if (Element(str(site.specie)) == Element(specie_u)) and (s not in site_indices_u):
                        foundSpecie = True
                        break
                if not foundSpecie:
                    raise ValueError(
                        "Could not find specie(s) in structure."
                    )
                site_indices_u.append(s)

    elif not site_indices_u:
        logger.warning(
            "Sites for computing U value are not specified. "
            "Computing U for first site in structure. "
        )

    site_indices_u = list(tuple(site_indices_u))
    num_perturb = len(site_indices_u)

    sites_perturb = []
    for site_index_perturb in site_indices_u:
        site = structure[site_index_perturb]
        sites_perturb.append(site)

    structure.remove_sites(indices=site_indices_u)

    for site in sites_perturb:
        structure.insert(i=0, species=site.specie, coords=site.frac_coords,
                         properties=site.properties)

    # using a uuid for book-keeping,
    # in a similar way to other workflows
    uuid = str(uuid4())

    c_defaults = {"vasp_cmd": VASP_CMD, "db_file": DB_FILE}
    if c:
        c.update(c_defaults)
    else:
        c = c_defaults

    # Calculate groundstate

    # ground state user incar settings
    uis_gs = uis.copy() # LMAXMIX
    if "LMAXMIX" not in uis_gs:
        logger.warning(
            "You have not specified LMAXMIX. Defaulting to VASP default."
        )
    uis_gs.update({"LDAU":False, "LORBIT": 11, "LWAVE":True, "ISPIN": 2})

    uis_ldau = uis_gs.copy()
    uis_ldau.update({"LDAU":True, "LDAUTYPE":3, "LDAUPRINT":2})

    # Initialize vasp input set
    vis_params = {"user_incar_settings": uis_ldau.copy()}
    vis_ldau = HubbardHundLinRespSet(structure=structure, num_perturb=num_perturb, **vis_params.copy())

    val_dict = {"LDAUL": {}, "LDAUU": {}, "LDAUJ": {}}
    for k in ["LDAUL", "LDAUU", "LDAUJ"]:
        if (k == "LDAUL"):
            # for LDAUL
            for i in range(num_perturb):
                val_dict[k].update({"perturb"+str(i):-1})
            for s in vis_ldau.poscar.site_symbols:
                l = -1
                val_dict[k].update({s:l})
        else:
            # for LDAUU and LDAUJ
            for i in range(num_perturb):
                val_dict[k].update({"perturb"+str(i):0})
            for s in vis_ldau.poscar.site_symbols:
                v = 0
                val_dict[k].update({s:v})
        uis_ldau.update({k:val_dict[k].copy()})
    vis_params = {"user_incar_settings": uis_ldau.copy()}
    vis_ldau = HubbardHundLinRespSet(structure=structure, num_perturb=num_perturb, **vis_params.copy())
    
    fws = []
    index_fw_gs = 0

    ediff_default = vis_ldau.incar['EDIFF']

    uis_gs = uis_ldau.copy()

    if relax_nonmagnetic:
        uis_gs.update({"ISPIN":1, "EDIFF":ediff_default})
    else:
        uis_gs.update({"ISPIN":2, "EDIFF":ediff_tight})
    vis_params = {"user_incar_settings": uis_gs.copy()}
    vis_gs = HubbardHundLinRespSet(structure=structure, num_perturb=num_perturb, **vis_params.copy())
    fws.append(HubbardHundLinRespFW(structure=structure,
                                 name="initial_static", vasp_input_set=vis_gs,
                                 vasp_cmd=VASP_CMD, db_file=DB_FILE))

    if relax_nonmagnetic:
        index_fw_gs += 1
        uis_gs.update({"ISPIN":2, "ICHARG":1, "EDIFF":ediff_tight})
        additional_files = ["WAVECAR", "CHGCAR"]
        vis_params = {"user_incar_settings": uis_gs.copy()}
        vis_gs = HubbardHundLinRespSet(structure=structure, num_perturb=num_perturb, **vis_params.copy())
        fws.append(HubbardHundLinRespFW(structure=structure, parents=fws[-1],
                                     additional_files=additional_files.copy(),
                                     name="initial_static_magnetic", vasp_input_set=vis_gs,
                                     vasp_cmd=VASP_CMD, db_file=DB_FILE))

    # Determine applied potential range
    if num_parallel_evals != num_total_evals:
        raise ValueError(
            "Different # of parallel & total evaluations not currently implemented."
        )

    applied_potential_value_list = []
    for counter_perturb in range(num_perturb):
        applied_potential_values = np.linspace(applied_potential_range[0],
                                               applied_potential_range[1], num_parallel_evals)
        applied_potential_values = np.around(applied_potential_values, decimals=9)

        if 0.0 in applied_potential_values:
            applied_potential_values = list(applied_potential_values)
            applied_potential_values.pop(applied_potential_values.index(0.0))
            applied_potential_values = np.array(applied_potential_values)

        applied_potential_value_list.append(applied_potential_values.copy())

    block_dict = {"s":0, "p":1, "d":2, "f":3}
    
    for counter_perturb in range(num_perturb):

        applied_potential_values = applied_potential_value_list[counter_perturb]

        for v in applied_potential_values:

            sign = 'neg' if str(v)[0] == '-' else 'pos'
            signs = []

            spin_potential_values = []
            if spin_polarized:
                spin_potential_values = [{"LDAUU": v, "LDAUJ": 0.0}, {"LDAUU": 0.0, "LDAUJ": v}]
                signs = [{"LDAUU": sign, "LDAUJ": ""}, {"LDAUU": "", "LDAUJ": sign}]
            else:
                spin_potential_values = [{"LDAUU": v, "LDAUJ": v}]
                signs = [{"LDAUU": sign, "LDAUJ": sign}]

            for spin_pot_dict, sign_dict in zip(spin_potential_values, signs):

                # Update perturbation potential for U and J
                for k in ["LDAUL", "LDAUU", "LDAUJ"]:
                    if (k == "LDAUL"):
                        # for LDAUL
                        # FIX ME: shouldn't hard code LDAUL = 2
                        for i in range(num_perturb):
                            if i == counter_perturb:
                                val_dict[k].update({"perturb"+str(i):block_dict[str(structure[counter_perturb].specie.block)]})
                            else:
                                val_dict[k].update({"perturb"+str(i):-1})
                    else:
                        # for LDAUU and LDAUJ
                        for i in range(num_perturb):
                            if i == counter_perturb:
                                val_dict[k].update({"perturb"+str(i):spin_pot_dict[k]})
                            else:
                                val_dict[k].update({"perturb"+str(i):0})
                    uis_ldau.update({k:val_dict[k].copy()})

                # Non-SCF runs
                uis_ldau.update({"ISTART":1, "ICHARG":11, "ISPIN":2, "EDIFF":ediff_default})

                vis_params = {"user_incar_settings": uis_ldau.copy()}
                vis_ldau = HubbardHundLinRespSet(structure=structure, num_perturb=num_perturb, **vis_params.copy())

                parents=fws[index_fw_gs]

                additional_files = ["WAVECAR", "CHGCAR"]

                fw = HubbardHundLinRespFW(structure=structure, parents=parents,
                                       name="nscf_site{}_vup_{}{}_vdn_{}{}".format(counter_perturb,
                                            sign_dict["LDAUU"], abs(round(spin_pot_dict["LDAUU"],6)),
                                            sign_dict["LDAUJ"], abs(round(spin_pot_dict["LDAUJ"],6))),
                                       vasp_input_set=vis_ldau,
                                       additional_files=additional_files.copy(),
                                       vasp_cmd=VASP_CMD, db_file=DB_FILE)
                fws.append(fw)

                # SCF runs
                uis_ldau.update({"ISTART":1, "ICHARG":0})
                if relax_nonmagnetic:
                    uis_ldau.update({"ISPIN":1, "EDIFF":ediff_default})
                else:
                    uis_ldau.update({"ISPIN":2, "EDIFF":ediff_tight})

                vis_params = {"user_incar_settings": uis_ldau.copy()}
                vis_ldau = HubbardHundLinRespSet(structure=structure, num_perturb=num_perturb, **vis_params.copy())

                if parallel_scheme == 0:
                    parents=fws[-1]
                else:
                    parents=fws[index_fw_gs]

                additional_files = ["WAVECAR"]

                fw = HubbardHundLinRespFW(structure=structure, parents=parents,
                                       name="scf_site{}_vup_{}{}_vdn_{}{}".format(counter_perturb,
                                            sign_dict["LDAUU"], abs(round(spin_pot_dict["LDAUU"],6)),
                                            sign_dict["LDAUJ"], abs(round(spin_pot_dict["LDAUJ"],6))),
                                       vasp_input_set=vis_ldau,
                                       additional_files=additional_files.copy(),
                                       vasp_cmd=VASP_CMD, db_file=DB_FILE)
                fws.append(fw)

                # SCF magnetic runs
                if relax_nonmagnetic:
                    uis_ldau.update({"ISTART":1, "ICHARG":1, "ISPIN":2, "EDIFF":ediff_tight})

                    vis_params = {"user_incar_settings": uis_ldau.copy()}
                    vis_ldau = HubbardHundLinRespSet(structure=structure, num_perturb=num_perturb, **vis_params.copy())
                    
                    parents=fws[-1]
                    additional_files = ["WAVECAR", "CHGCAR"]

                    fw = HubbardHundLinRespFW(structure=structure, parents=parents,
                                           name="scf_magnetic_site{}_vup_{}{}_vdn_{}{}".format(counter_perturb,
                                                sign_dict["LDAUU"], abs(round(spin_pot_dict["LDAUU"],6)),
                                                sign_dict["LDAUJ"], abs(round(spin_pot_dict["LDAUJ"],6))),
                                           vasp_input_set=vis_ldau,
                                           additional_files=additional_files.copy(),
                                           vasp_cmd=VASP_CMD, db_file=DB_FILE)
                    fws.append(fw)

    wf = Workflow(fws)

    fw_analysis = Firework(
        HubbardHundLinRespToDb(
            num_perturb=num_perturb,
            spin_polarized=spin_polarized,
            relax_nonmagnetic=relax_nonmagnetic,
            db_file=DB_FILE, wf_uuid=uuid
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

class PoscarPerturb(Poscar):
    """
    Derived Poscar class that allows the distinction of individual sites in the Structure
    """

    def __init__(
            self,
            structure: Structure,
            num_perturb: int = 1,
            **kwargs
    ):
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

        if (self.num_perturb > 0 and self.num_perturb <= len(self.structure)):
            syms = [site.specie.symbol for site in self.structure[self.num_perturb:]]
            syms = [a[0] for a in itertools.groupby(syms)]
            syms_perturb = [site.specie.symbol for site in self.structure[0:self.num_perturb]]
            syms_perturb.extend(syms)
        else:
            raise ValueError(
                "Invalid atom index to perturb"
            )

        return syms_perturb

    @property
    def natoms(self):
        """
        Sequence of number of sites of each type associated with the Poscar
        """

        if (self.num_perturb > 0 and self.num_perturb <= len(self.structure)):
            syms = [site.specie.symbol for site in self.structure[self.num_perturb:]]
            n_atoms = [len(tuple(a[1])) for a in itertools.groupby(syms)]

            n_atoms_perturb = [1 for i in range(self.num_perturb)]
            n_atoms_perturb.extend(n_atoms)
        else:
            raise ValueError(
                "Invalid atom index to perturb"
            )

        return n_atoms_perturb


class HubbardHundLinRespSet(MPStaticSet):
    """
    VASP input set for Linear Response Hubbard U workflow perturbation Fireworks
    """
    def __init__(self, structure, num_perturb, prev_incar=None, prev_kpoints=None,
                 reciprocal_density=100, small_gap_multiply=None, **kwargs):
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
        
        incar.update({"ISYM": -1, "ISMEAR": 0, "LREAL":False, "LASPH":True})
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
        poscar = PoscarPerturb(structure=super().structure, num_perturb=self.num_perturb)
        return poscar

    @property
    def kpoints(self):
        """
        Custom Kpoints for HubbardHundLinRespSet
        """
        kpoints = super().kpoints

        return kpoints

class HubbardHundLinRespFW(Firework):
    def __init__(self, structure=None, name="linresponse_U",
                 vasp_input_set=None, vasp_input_set_params=None,
                 vasp_cmd=VASP_CMD, prev_calc_loc=True, prev_calc_dir=None,
                 db_file=DB_FILE, vasptodb_kwargs=None, parents=None,
                 additional_files=None,
                 **kwargs):
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

        vasp_input_set_params = vasp_input_set_params or {}
        vasptodb_kwargs = vasptodb_kwargs or {}
        if "additional_fields" not in vasptodb_kwargs:
            vasptodb_kwargs["additional_fields"] = {}
        vasptodb_kwargs["additional_fields"]["task_label"] = name

        fw_name = "{}-{}".format(structure.composition.reduced_formula if structure else "unknown", name)

        if prev_calc_dir:
            t.append(CopyVaspOutputs(calc_dir=prev_calc_dir, additional_files=additional_files,
                                     contcar_to_poscar=False))
        elif parents:
            if prev_calc_loc:
                t.append(CopyVaspOutputs(calc_loc=prev_calc_loc, additional_files=additional_files,
                                         contcar_to_poscar=False))

        if structure:
            vasp_input_set = vasp_input_set or HubbardHundLinRespSet(structure, **vasp_input_set_params)
            t.append(WriteVaspFromIOSet(structure=structure,
                                        vasp_input_set=vasp_input_set))
        else:
            raise ValueError("Must specify structure")

        t.append(RunVaspCustodian(vasp_cmd=vasp_cmd, auto_npar=">>auto_npar<<"))
        t.append(PassCalcLocs(name=name))
        t.append(
            VaspToDb(db_file=db_file, **vasptodb_kwargs))
        super(HubbardHundLinRespFW, self).__init__(t, parents=parents, name=fw_name, **kwargs)
