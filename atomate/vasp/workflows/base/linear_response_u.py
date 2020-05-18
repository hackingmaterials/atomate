# coding: utf-8

import os

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
from atomate.vasp.firetasks.parse_outputs import VaspToDb, LinearResponseUToDb
from atomate.vasp.firetasks.write_inputs import WriteVaspFromIOSet

from atomate.utils.utils import get_logger

logger = get_logger(__name__)

from atomate.vasp.config import VASP_CMD, DB_FILE, ADD_WF_METADATA

from atomate.vasp.workflows.presets.scan import wf_scan_opt
from uuid import uuid4
from pymatgen.io.vasp.sets import MPRelaxSet, MPStaticSet
from pymatgen.io.vasp.inputs import Poscar, Incar
from pymatgen.core import Lattice, Structure

import numpy as np

__author__ = ""
__maintainer__ = ""
__email__ = ""
__status__ = "Production"
__date__ = "February 2020"

__linear_response_u_wf_version__ = 0.0

module_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)))

def get_wf_linear_response_u(structure,
                             applied_potential_range=[-0.2, 0.2],
                             num_parallel_evals=9, num_total_evals=9,
                             site_indices_u=None, species_u=None,
                             use_default_uvals=False,
                             ground_state_ldau=True, ground_state_dir=None,
                             c=None, vis=None):
    """
    Args:
        structure: 
        c: Workflow config dict, in the same format
    as in presets/core.py and elsewhere in atomate
        vis: A VaspInputSet to use for the first FW

    Returns: Workflow
    """

    if not structure.is_ordered:
        raise ValueError(
            "Please obtain an ordered approximation of the input structure."
        )

    # Reorder structure
    if not site_indices_u:
        site_indices_u = []
    
    if species_u:
        foundSpecie = False
        for s in range(len(structure)):
            site = structure[s]
            if site.specie == Element(target_elem):
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
    uis_gs = {"LDAU":False, "LMAXMIX":4, "LORBIT": 11, "ISPIN": 2}

    uis_ldau = uis_gs.copy()
    uis_ldau.update({"LDAU":True, "LDAUTYPE":3, "LDAUPRINT":2})

    # Load default U values
    vis_params = {"user_incar_settings":{}}
    set_default_uvals = MPRelaxSet(structure=structure, sort_structure=False, **vis_params.copy())
    incar_dict_default_u = set_default_uvals.incar.as_dict()
    sitesym_default_u = set_default_uvals.poscar.site_symbols
    default_uvals = {}
    if 'LDAUU' in incar_dict_default_u.keys():
        uvals = incar_dict_default_u['LDAUU']
        lvals = incar_dict_default_u['LDAUL']
        for sym, u, l in zip(sitesym_default_u, uvals, lvals):
            default_uvals.update({sym:{'LDAUU':u, 'LDAUL':l}})

    # Initialize vasp input set
    vis_params = {"user_incar_settings": uis_ldau.copy()}
    vis_ldau = LinearResponseUSet(structure=structure, **vis_params.copy())

    for k in ["LDAUL", "LDAUU", "LDAUJ"]:
        val_dict = {}
        if (k == "LDAUL"):
            # for LDAUL
            val_dict.update({"perturb":2})             # FIXME: shouldn't hard code LDAUL
            for s in vis_ldau.poscar.site_symbols:
                l = -1
                if use_default_uvals:
                    if s in default_uvals.keys():
                        if k in default_uvals[s].keys():
                            l = default_uvals[s][k]
                val_dict.update({s:l})
            uis_ldau.update({k:val_dict.copy()})
        else:
            # for LDAUU and LDAUJ
            val_dict.update({"perturb":0})
            for s in vis_ldau.poscar.site_symbols:
                v = 0
                if use_default_uvals:
                    if s in default_uvals.keys():
                        if 'LDAUU' in default_uvals[s].keys():
                            v = default_uvals[s]['LDAUU']
                val_dict.update({s:v})
            uis_ldau.update({k:val_dict.copy()})

    if ground_state_ldau:
        uis_gs = uis_ldau.copy()
        vis_params = {"user_incar_settings": uis_gs.copy()}
        vis_gs = LinearResponseUSet(structure=structure, **vis_params.copy())
        fw_gs = LinearResponseUFW(structure=structure, name="initial static", vasp_input_set=vis_gs,
                                  vasp_cmd=VASP_CMD, db_file=DB_FILE)
    else:
        vis_params = {"user_incar_settings": uis_gs.copy()}
        vis_gs = MPStaticSet(structure=structure, sort_structure=False, **vis_params.copy())
        fw_gs = StaticFW(structure=structure, name="initial static", vasp_input_set=vis_gs,
                         vasp_cmd=VASP_CMD, db_file=DB_FILE)

    if ground_state_dir:
        fws = []
    else:
        fws = [fw_gs]

    # Determine applied potential range
    if num_parallel_evals != num_total_evals:
        raise ValueError(
            "Different # of parallel & total evaluations not currently implemented."
        )

    applied_potential_values = np.linspace(applied_potential_range[0],
                                           applied_potential_range[1], num_parallel_evals)
    applied_potential_values = np.around(applied_potential_values, decimals=9)

    if 0.0 in applied_potential_values:
        applied_potential_values = list(applied_potential_values)
        applied_potential_values.pop(applied_potential_values.index(0.0))
        applied_potential_values = np.array(applied_potential_values)

    for v in applied_potential_values:

        sign = 'neg' if str(v)[0] == '-' else 'pos'

        # Update applied potential to U and J
        uis_ldau.update({"ISTART":1})

        # Update perturbation potential for U and J
        for k in ["LDAUU", "LDAUJ"]:
            # for LDAUU and LDAUJ
            val_dict.update({"perturb":v})
            uis_ldau.update({k:val_dict.copy()})

        # Non-SCF runs
        uis_ldau.update({"ICHARG":11})

        vis_params = {"user_incar_settings": uis_ldau.copy()}
        vis_ldau = LinearResponseUSet(structure=structure, **vis_params.copy())

        if ground_state_dir:
            parents = []
        else:
            parents=fws[0]

        fw = LinearResponseUFW(structure=structure, parents=parents,
                               name="nscf_u_eq_{}{}".format(sign, abs(round(v,6))),
                               vasp_input_set=vis_ldau,
                               additional_files=["WAVECAR","CHGCAR"],
                               prev_calc_dir=ground_state_dir,
                               vasp_cmd=VASP_CMD, db_file=DB_FILE)

        fws.append(fw)

        # SCF runs
        uis_ldau.update({"ICHARG":0})

        vis_params = {"user_incar_settings": uis_ldau.copy()}
        vis_ldau = LinearResponseUSet(structure=structure, **vis_params.copy())

        # NOTE: More efficient to reuse WAVECAR or remove dependency of SCF on NSCF?
        if ground_state_dir:
            parents = []
        else:
            parents=fws[0]
            # parents=fws[-1]

        fw = LinearResponseUFW(structure=structure, parents=parents,
                               name="scf_u_eq_{}{}".format(sign, abs(round(v,6))),
                               vasp_input_set=vis_ldau,
                               additional_files=["WAVECAR"],
                               prev_calc_dir=ground_state_dir,
                               vasp_cmd=VASP_CMD, db_file=DB_FILE)
        fws.append(fw)

    wf = Workflow(fws)

    fw_analysis = Firework(
        LinearResponseUToDb(
            db_file=DB_FILE, wf_uuid=uuid
        ),
        name="LinearResponseUToDb",
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
                "wf_name": "linear_response_u",
                "wf_version": __linear_response_u_wf_version__,
            }
        },
    )

    return wf

class PoscarPerturb(Poscar):
    """
    FILL
    """

    def __init__(
            self,
            structure: Structure,
            perturb_index: int = 0,
            comment: str = None,
            selective_dynamics=None,
            true_names: bool = True,
            velocities=None,
            predictor_corrector=None,
            predictor_corrector_preamble=None,
            sort_structure: bool = False,
    ):
        """
        FILL
        """
        # super().__init__(structure=Structure)

        self.perturb_index = perturb_index

        if structure.is_ordered:
            site_properties = {}
            if selective_dynamics:
                site_properties["selective_dynamics"] = selective_dynamics
            if velocities:
                site_properties["velocities"] = velocities
            if predictor_corrector:
                site_properties["predictor_corrector"] = predictor_corrector
            structure = Structure.from_sites(structure)
            self.structure = structure.copy(site_properties=site_properties)
            if sort_structure:
                self.structure = self.structure.get_sorted_structure()
            self.true_names = true_names
            self.comment = structure.formula if comment is None else comment
            self.predictor_corrector_preamble = predictor_corrector_preamble
        else:
            raise ValueError(
                "Structure with partial occupancies cannot be " "converted into POSCAR!"
            )

        self.temperature = -1

    @property
    def site_symbols(self):
        """
        Sequence of symbols associated with the Poscar. Similar to 6th line in
        vasp 5+ POSCAR.
        """

        syms = super().site_symbols

        if (self.perturb_index == 0):
            syms_perturb = []
            if ((syms[0] == syms[1]) & (len(syms) > 1)):
                syms_perturb = [syms[0]]
            syms_perturb.extend(syms)
        else:
            raise ValueError(
                "Invalid atom index to perturb"
            )

        return syms_perturb

    @property
    def sites(self):
        """
        FILL
        """
        sites_array = [site for site in self.structure]
        return sites_array

    @property
    def natoms(self):
        """
        Sequence of number of sites of each type associated with the Poscar.
        Similar to 7th line in vasp 5+ POSCAR or the 6th line in vasp 4 POSCAR.
        """

        if (super().natoms[self.perturb_index] > 1):
            if (self.perturb_index == 0):
                n_atoms = [1]
                n_atoms.extend(super().natoms)
                n_atoms[1] -= 1
            else:
                raise ValueError(
                    "Invalid atom index to perturb"
                )
        else:
            n_atoms = super().natoms

        return n_atoms


class LinearResponseUSet(MPStaticSet):
    """
    FILL
    """
    def __init__(self, structure, prev_incar=None, prev_kpoints=None,
                 lepsilon=False, lcalcpol=False, reciprocal_density=100,
                 small_gap_multiply=None, **kwargs):
        """
        FILL
        """

        super().__init__(structure, sort_structure=False, **kwargs)

        if isinstance(prev_kpoints, str):
            prev_kpoints = Kpoints.from_file(prev_kpoints)
        self.prev_kpoints = prev_kpoints

        self.reciprocal_density = reciprocal_density
        self.kwargs = kwargs
        self.lepsilon = lepsilon
        self.lcalcpol = lcalcpol
        self.small_gap_multiply = small_gap_multiply

    @property
    def incar(self):
        """
        FILL
        """
        parent_incar = super().incar
        settings = dict(self._config_dict["INCAR"])

        structure = self.structure

        incar = Incar(parent_incar)

        settings.pop("LDAUU", None)
        settings.pop("LDAUJ", None)
        settings.pop("LDAUL", None)

        # Note that DFPT calculations MUST unset NSW. NSW = 0 will fail
        # to output ionic.

        settings.pop("NSW", None)
        incar.pop("NSW", None)

        incar.update({"ISYM": -1, "IBRION": -1, "LCHARG": True, "LWAVE": True})
        # "LORBIT": 11, "LVHAR": True, "LAECHG": True

        for k, v in settings.items():
            if k == "MAGMOM":
                mag = []
                for site in structure:
                    if hasattr(site, 'magmom'):
                        mag.append(site.magmom)
                    elif hasattr(site.specie, 'spin'):
                        mag.append(site.specie.spin)
                    elif str(site.specie) in v:
                        mag.append(v.get(str(site.specie)))
                    else:
                        mag.append(v.get(site.specie.symbol, 0.6))
                incar[k] = mag
            elif k.startswith("EDIFF") and k != "EDIFFG":
                if "EDIFF" not in settings and k == "EDIFF_PER_ATOM":
                    incar["EDIFF"] = float(v) * structure.num_sites
                else:
                    incar["EDIFF"] = float(settings["EDIFF"])
            else:
                incar[k] = v

        for k in ["MAGMOM", "NUPDOWN"] + list(self.kwargs.get(
                "user_incar_settings", {}).keys()):
            # For these parameters as well as user specified settings, override
            # the incar settings.
            if parent_incar.get(k, None) is not None:
                incar[k] = parent_incar[k]
            else:
                incar.pop(k, None)

        if incar.get('LDAU'):
            # ensure to have LMAXMIX for GGA+U static run
            if "LMAXMIX" not in incar:
                incar.update({"LMAXMIX": parent_incar["LMAXMIX"]})

        # Compare ediff between previous and staticinputset values,
        # choose the tighter ediff
        incar["EDIFF"] = min(incar.get("EDIFF", 1), parent_incar["EDIFF"])
        
        if self.kwargs.get("user_incar_settings")["LDAUU"]:

            # Need to add another parameter for perturbed atom

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
        FILL
        """
        poscar = PoscarPerturb(structure=super().structure)
        return poscar

    @property
    def kpoints(self):
        """
        FILL
        """
        self._config_dict["KPOINTS"]["reciprocal_density"] = self.reciprocal_density
        kpoints = super().kpoints

        # Prefer to use k-point scheme from previous run
        # except for when lepsilon = True is specified
        if self.prev_kpoints and self.prev_kpoints.style != kpoints.style:
            if (self.prev_kpoints.style == Kpoints.supported_modes.Monkhorst) \
               and (not self.lepsilon):
                k_div = [kp + 1 if kp % 2 == 1 else kp
                         for kp in kpoints.kpts[0]]
                kpoints = Kpoints.monkhorst_automatic(k_div)
            else:
                kpoints = Kpoints.gamma_automatic(kpoints.kpts[0])
        return kpoints

class LinearResponseUFW(Firework):
    def __init__(self, structure=None, name="linresponse_U",
                 vasp_input_set=None, vasp_input_set_params=None,
                 vasp_cmd=VASP_CMD, prev_calc_loc=True, prev_calc_dir=None,
                 db_file=DB_FILE, vasptodb_kwargs=None, parents=None,
                 additional_files=None,
                 is_nscf=False,
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
            vasp_input_set = vasp_input_set or LinearResponseUSet(structure, **vasp_input_set_params)
            t.append(WriteVaspFromIOSet(structure=structure,
                                        vasp_input_set=vasp_input_set))
        else:
            raise ValueError("Must specify structure")

        t.append(RunVaspCustodian(vasp_cmd=vasp_cmd, auto_npar=">>auto_npar<<"))
        t.append(PassCalcLocs(name=name))
        t.append(
            VaspToDb(db_file=db_file, **vasptodb_kwargs))
        super(LinearResponseUFW, self).__init__(t, parents=parents, name=fw_name, **kwargs)
