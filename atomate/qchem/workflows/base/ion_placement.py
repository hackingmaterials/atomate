# coding: utf-8

from __future__ import absolute_import, division, print_function, \
    unicode_literals

from fireworks import Workflow
from atomate.qchem.fireworks.core import FrequencyFlatteningOptimizeFW, PlaceIonFW, GatherGeomsFW
from atomate.utils.utils import get_logger

__author__ = "Samuel Blau"
__copyright__ = "Copyright 2018, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Samuel Blau"
__email__ = "samblau1@gmail.com"
__status__ = "Alpha"
__date__ = "11/13/18"

logger = get_logger(__name__)


def get_ion_placement_wf(molecule,
                         mulliken=None,
                         ion=None,
                         charges=None,
                         stop_num=None,
                         do_triplets=True,
                         pcm_dielectric=None,
                         do_optimization=True,
                         linked=False,
                         qchem_input_params=None,
                         name="place ions and optimize then gather",
                         db_file=">>db_file<<",
                         test_positions=None,
                         ref_dirs=None,
                         **kwargs):
    """
    """

    fw_name = "place ions and optimize"
    gather_name = "gather geometries"
    if name != "place ions and optimize then gather":
        fw_name = name
        gather_name = "_"+name

    qchem_input_params = qchem_input_params or {}
    if pcm_dielectric != None:
        qchem_input_params["pcm_dielectric"] = pcm_dielectric

    if do_optimization:
        # Optimize the original molecule
        fw1 = FrequencyFlatteningOptimizeFW(
            molecule=molecule,
            name="first FF",
            qchem_cmd=">>qchem_cmd<<",
            max_cores=">>max_cores<<",
            qchem_input_params=qchem_input_params,
            linked=linked,
            db_file=db_file)

        # Place ions around the optimized molecule and optimize the resulting structures
        fw2 = PlaceIonFW(
            ion=ion,
            charges=charges,
            stop_num=stop_num,
            do_triplets=do_triplets,
            linked=linked,
            name=fw_name,
            qchem_input_params=qchem_input_params,
            test_positions=test_positions,
            ref_dirs=ref_dirs,
            parents=fw1)

        # Gather all optimized ion+mol structures and identify the unique minima
        fw3 = GatherGeomsFW(
            prefix="ion_pos_",
            db_file=db_file,
            name=gather_name,
            parents=fw2)

        fws = [fw1, fw2, fw3]

    else:
        # Place ions around the given molecule and optimize the resulting structures
        fw1 = PlaceIonFW(
            molecule=molecule,
            mulliken=mulliken,
            ion=ion,
            charges=charges,
            stop_num=stop_num,
            do_triplets=do_triplets,
            linked=linked,
            name=fw_name,
            qchem_input_params=qchem_input_params,
            test_positions=test_positions,
            ref_dirs=ref_dirs)

        # Gather all optimized ion+mol structures and identify the unique minima
        fw2 = GatherGeomsFW(
            prefix="ion_pos_",
            db_file=db_file,
            name=gather_name,
            parents=fw1)

        fws = [fw1, fw2]

    wfname = "{}:{}".format(molecule.composition.reduced_formula, name)

    return Workflow(fws, name=wfname, **kwargs)
