# coding: utf-8

from __future__ import absolute_import, division, print_function, unicode_literals

import six
from abc import ABCMeta, abstractmethod

__author__ = "Kiran Mathew"
__email__ = "kmathew@lbl.gov"


class AbstractBuilder(six.with_metaclass(ABCMeta)):
    """
    Abstract builder class. Defines the contract and must be subclassed by all builders.
    """

    @abstractmethod
    def run(self):
        """
        Run the builder.
        """
        pass

    @abstractmethod
    def reset(self):
        """
        Unset the building.
        """
        pass

    @classmethod
    @abstractmethod
    def from_file(cls, filename):
        """
        Set the builder from a db file
        """
        pass
