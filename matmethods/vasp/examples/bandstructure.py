from matmethods.vasp.workflows.automatic.core import wf_band_structure
from pymatgen import Structure

__author__ = 'Anubhav Jain, Kiran Mathew'
__email__ = 'ajain@lbl.gov, kmathew@lbl.gov'


if __name__ == "__main__":
    struct_si = {u'lattice': {u'a': 3.839943374653261, u'c': 3.83994338,
                              u'b': 3.839943378813096,
                              u'matrix': [[3.32548851, 0.0, 1.91997169],
                                          [1.10849617, 3.13530064, 1.91997169],
                                          [0.0, 0.0, 3.83994338]],
                              u'volume': 40.036809671145996,
                              u'beta': 59.99999995393976,
                              u'gamma': 60.00000000512866,
                              u'alpha': 59.99999998977525}, u'sites': [
        {u'abc': [0.875, 0.875, 0.875],
         u'xyz': [3.879736595, 2.74338806, 6.719900914999999], u'label': u'Si',
         u'species': [{u'occu': 1, u'element': u'Si'}]},
        {u'abc': [0.125, 0.125, 0.125],
         u'xyz': [0.554248085, 0.39191258, 0.959985845], u'label': u'Si',
         u'species': [{u'occu': 1, u'element': u'Si'}]}],
                 u'@class': 'Structure',
                 u'@module': 'pymatgen.core.structure'}

    structure = Structure.from_dict(struct_si)
    wf = wf_band_structure(structure)