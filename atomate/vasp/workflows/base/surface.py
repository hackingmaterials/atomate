# coding: utf-8


## for Surface Energy Calculation
from __future__ import division, unicode_literals

"""
If you use this module, please consider citing the following work::

    R. Tran, Z. Xu, B. Radhakrishnan, D. Winston, W. Sun, K. A. Persson,
    S. P. Ong, "Surface Energies of Elemental Crystals", Scientific Data,
    2016, 3:160080, doi: 10.1038/sdata.2016.80.

as well as::

    Sun, W.; Ceder, G. Efficient creation and convergence of surface slabs,
    Surface Science, 2013, 617, 53–59, doi:10.1016/j.susc.2013.05.016.
"""

__author__ = "Richard Tran"
__version__ = "0.2"
__email__ = "rit001@eng.ucsd.edu"
__date__ = "11/30/17"

import os

from pymatgen import Structure

from fireworks import Workflow


class SurfaceWorkflowManager(object):
    """
        Workflow manager with a set of common database and calculation specs for
            all workflows. The workflow will use VASP to ultimately calculate the
            results needed to derive the surface energy and work function. Other
            properties such as the Wulff shape, weighted surface energy and weighted
            work functions can then be derived from these basic properties.

        The user has the option of started their calculations from the conventional
            unit cell which will calculate all oriented unit cells, terminations,
            and reconstructions of slabs, or from an oriented unit cell which will
            calculate all possible terminations for a specific facet or from a
            single slab structure which will just calculate the given slab. The
            use of the oriented unit cell calculations provides us with well
            converged surface quantities.

    .. attribute:: k_product

        The k-point mesh is the reciprocal of the
            lattice LENGTH multiplied by this integer.

    .. attribute:: vasp_cmd

        Command for running vasp.

    .. attribute:: db_file

        File containing the specifications for the
        database to run and insert calculations.

    .. attribute:: scratch_dir

        The directory to run the calculations in.

    .. attribute:: cwd

        Directory to run the workflow and store the final outputs.

    """

    def __init__(self, db_file, cwd=os.getcwd(),
                 scratch_dir="", k_product=50, vasp_cmd="vasp"):

        """
        Initializes the workflow manager with common database and calculations specs

        Args:
            db_file (str): Location of file containing database specs.
            cwd (str): Location of directory to operate the
                workflow and store the final outputs
            scratch_dir (str): if specified, uses this directory as the root
                scratch dir. Supports env_chk.
            k_product (int): Kpoint number * length for a & b directions, also for c
                direction in bulk calculations. Default to 40.
            vasp_cmd (str): Command used to run vasp.
        """

        self.k_product = k_product
        self.vasp_cmd = vasp_cmd
        self.db_file = db_file
        self.scratch_dir = scratch_dir
        self.cwd = cwd

    def from_conventional_unit_cell(self, structure, mmi, mpid="--"):
        """
        Calculate surface properties from a conventional unit cell. This workflow
            will continue running calculations until all possible slabs up to a max
            miller index of mmi has been completed.

        Args:
            structure (Structure): The conventional unit cell.
            mmi (int): Max Miller index.
            mpid (str): Materials Project ID of the conventional unit cell.
        """

        return Workflow([SurfCalcOptimizer(structure, self.scratch_dir,
                                           self.k_product, self.db_file,
                                           self.vasp_cmd, "conventional_unit_cell",
                                           mmi=mmi, cwd=self.cwd, mpid=mpid)])

    def from_oriented_unit_cell(self, structure, miller_index, scale_factor,
                                reconstruction=None, mpid="--"):
        """
        Calculate surface properties from an oriented unit cell. This workflow will run
            calculations on all possible terminations and reconstructions for a specific
            facet (miller index).

        Args:
            structure (Structure): Oriented unit cell structure.
            miller_index ([h, k, l]): Miller index of plane parallel to
                surface (and oriented unit cell).
            scale_factor (array): Final computed scale factor that brings
                the parent cell to the surface cell.
            mpid (str): Materials Project ID of the conventional unit cell.
        """

        return Workflow([SurfCalcOptimizer(structure, self.scratch_dir, self.k_product,
                                           self.db_file, self.vasp_cmd,
                                           "oriented_unit_cell", cwd=self.cwd,
                                           miller_index=miller_index, mpid=mpid,
                                           reconstruction=reconstruction,
                                           scale_factor=scale_factor, **kwargs)])

    def from_slab_cell(self, structure, miller_index, shift, scale_factor,
                       ouc, ssize, vsize, reconstruction=None, mpid="--"):
        """
        Calculates the surface properties of a single slab structure

        Args:
            structure (Structure): Slab structure.
            miller_index ([h, k, l]): Miller index of plane parallel to
                surface (and oriented unit cell).
            ouc (Structure): The oriented_unit_cell from which
                this Slab is created (by scaling in the c-direction).
            shift (float): The shift in the c-direction applied to get the
                termination.
            ssize (float): Minimum slab size in Angstroms or number of hkl planes
            vsize (float): Minimum vacuum size in Angstroms or number of hkl planes
            scale_factor (array): Final computed scale factor that brings
                the parent cell to the surface cell.
            mpid (str): Materials Project ID of the conventional unit cell.
            reconstruction (str): The name of the reconstruction
                (if it is a reconstructed slab).
        """

        return Workflow([SurfCalcOptimizer(structure, self.scratch_dir, self.k_product,
                                           self.db_file, self.vasp_cmd, "slab_cell",
                                           miller_index=miller_index, ouc=ouc, shift=shift,
                                           scale_factor=scale_factor, ssize=ssize,
                                           reconstruction=reconstruction, cwd=self.cwd,
                                           vsize=vsize, mpid=mpid, **kwargs)])


