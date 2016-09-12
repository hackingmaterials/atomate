# coding: utf-8

from __future__ import division, print_function, unicode_literals, absolute_import

"""
This module defines tasks for writing FEFF input sets.
"""

from fireworks import FireTaskBase, explicit_serialize

__author__ = 'Kiran Mathew'
__email__ = 'kmathew@lbl.gov'


def load_class(mod, name):
    mod = __import__(mod, globals(), locals(), [name], 0)
    return getattr(mod, name)


@explicit_serialize
class WriteFeffFromIOSet(FireTaskBase):
    """
    """

    required_params = ["absorbing_atom", "structure", "feff_input_set"]
    optional_params = ["radius", "other_params"]

    def run_task(self, fw_spec):
        # if a full object is provided.
        if hasattr(self['feff_input_set'], 'write_input'):
            fis = self['feff_input_set']

        # if inputset String + parameters was provided
        else:
            fis_cls = load_class("pymatgen.io.feff.sets", self["feff_input_set"])
            fis = fis_cls(self["absorbing_atom"], self["structure"], self.get("radius", 10.0),
                          **self.get("other_params", {}))
        print(fis.tags)
        fis.write_input(".")
