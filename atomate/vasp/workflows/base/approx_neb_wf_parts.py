from fireworks import Firework, FWAction, Workflow, FiretaskBase
from atomate.common.firetasks.glue_tasks import PassCalcLocs
from atomate.vasp.firetasks.glue_tasks import pass_vasp_result
from atomate.vasp.firetasks.run_calc import RunVaspCustodian
from atomate.vasp.firetasks.write_inputs import WriteVaspFromIOSet
from atomate.vasp.firetasks.parse_outputs import VaspToDb
from atomate.vasp.fireworks.core import OptimizeFW
from atomate.vasp.config import VASP_CMD, DB_FILE
from pymatgen.io.vasp.sets import MPRelaxSet
from pymatgen import Structure


@explicit_serialize
class InsertSites(FiretaskBase):
    """
    Insert sites into an output structure for use of the modified structure.

    Args:
        base_name (str): Name of the file to create.
        format (str): Optional. one of "zip", "tar", "bztar" or "gztar".
        """

    _fw_name = "InsertSites"
    required_params = ["insert_specie", "insert_coords"]
    optional_params = []

    def run_task(self, fw_spec):
        structure = Structure.from_dict(fw_spec["host_lattice_structure"])
        insert_coords = self["insert_coords"]
        insert_specie = self["insert_specie"]

        if structure.site_properties != {}:  # removes site properties to avoid error
            for p in structure.site_properties.keys():
                structure.remove_site_property(p)
        for coords in sorted(insert_coords, reverse=True):
            structure.insert(
                0, insert_specie, coords, **kwargs
            )  # add kwarg for coords_are_cartesian=False
        return structure



class InsertSitesFW(Firework):
    # TODO: Write class description
    def __init__(
        self,
        structure,
        insert_specie,
        insert_coords,
        name="approx neb insert working ion",
        vasp_input_set=None,
        override_default_vasp_params={},
        vasp_cmd=VASP_CMD,
        db_file=DB_FILE,
        parents=None,
        **kwargs
    ):
        # if structure == None and parents == None:
        #   print("ERROR")
        # elif structure == None: #setting structure supercedes parents
        #   connect to database
        #   query for parent using fw_spec['_job_info'][-1]['launch_dir']
        #   get structure...
        #   structure #from parents
        #TODO:Is pass structure needed in this FW? How to ensure pass_dict key matches?
        pass_structure_fw = pass_vasp_result(
            pass_dict={
                "host_lattice_structure": ">>output.structure"
            }
        )
        structure = InsertSites(insert_specie=insert_specie,insert_coords=insert_coords)
        vasp_input_set = vasp_input_set or MPRelaxSet(structure, **override_default_vasp_params)
        t = [pass_structure_fw]
        t.append(WriteVaspFromIOSet(structure=structure, vasp_input_set=vasp_input_set))
        t.append(RunVaspCustodian(vasp_cmd=vasp_cmd, job_type="double_relaxation_run"))
        t.append(PassCalcLocs(name=name))
        t.append(VaspToDb(db_file=db_file, additional_fields={"task_label": name}))


class PathFinderFW(Firework):
    #TODO: Write PathfinderFW
    def add_fix_two_atom_selective_dynamics(structure, fixed_index=0):
        """
        Returns structure with selective dynamics assigned to fix the position of two sites.
        Two sites will be fixed: 1) the site specified by fixed_index and 2) the site positioned furthest from the specified fixed_index site.

        Args:
            structure (Structure): Input structure (e.g. host lattice with one working ion intercalated)
            fixed_index (int): Index of site in structure whose position will be fixed (e.g. working ion site)
        Returns:
            Structure
        """
        sd_structure = structure.copy()
        sd_array = [[True, True, True]] * sd_structure.num_sites
        sd_array[fixed_index] = [False, False, False]
        ref_site = sd_structure.sites[fixed_index]
        distances = [site.distance(ref_site) for site in sd_structure.sites]
        farthest_index = distances.index(max(distances))
        sd_array[farthest_index] = [False, False, False]
        sd_structure.add_site_property("selective_dynamics", sd_array)
        return sd_structure
