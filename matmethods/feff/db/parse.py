
#!/usr/bin/env python

"""
TODO: Modify module doc.
"""

from __future__ import division

__author__ = "Kiran Mathew, Chen Zheng"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Chen Zheng"
__email__ = "chz022@ucsd.edu"
__date__ = "10/06/16"

from matgendb.util import get_database


class FEFFDBManager(object):

    def __init__(self,settings_file,admin=False):
        """
        Init the db manager
        Returns: FEFFDBManager object

        """

        self.db = self.get_db(settings_file,admin)


    def get_db(cls, settings_file, admin=True):
        """
        TODO filled in document
        Args:
            settings_file:
            admin:

        Returns:

        """

        return get_database(settings_file,admin=admin)

    def get_eels_param(self,eels_index):
        """
        TODO filled in documents
        Args:
            eels_index:

        Returns:

        """

        sample = self.db.eelsdb.find({'eels_index':eels_index})
        convg_angle = float(sample["meta_data"]["Convergence Semi-angle"].split()[0])
        coll_angle = float(sample["meta_data"]["Collection Semi-angle"].split()[0])
        beam_energy = float(sample["meta_data"]["Incident Beam Energy"].split()[0])
        absorber = sample["meta_data"]["Elemental Edges"].split("_")[0]
        edge = sample["meta_data"]["Elemental Edges"].split("_")[1]
        return convg_angle, coll_angle, beam_energy, absorber, edge

    def get_xafs_param(self,eels_index):
        """
        TODO filled in document
        Args:
            eels_index:

        Returns:

        """

        sample_data = self.get_eels_param(eels_index)
        absorber, edge = sample_data[3],sample_data[4]
        return absorber, edge

    def get_mp_id(self,eels_index):

        """

        Args:
            eels_index:

        Returns: mp_id entry dictionary object correspoinding to eels_index

        """

        sample = self.db.eelsdb.find({'eels_index':eels_index})

        return sample["mp_id"]





