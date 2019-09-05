# coding: utf-8

from abc import ABCMeta, abstractmethod

__author__ = "Kiran Mathew"
__email__ = "kmathew@lbl.gov"


class AbstractBuilder(metaclass=ABCMeta):
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
        Set the builder from a config file, e.g., a db file
        """
        pass
