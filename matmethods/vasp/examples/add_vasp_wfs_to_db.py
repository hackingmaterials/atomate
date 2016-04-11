from fireworks import LaunchPad
from matmethods.vasp.workflows.automatic.standard import wf_band_structure
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
                 u'@module': 'pymatgen.core.structure'}  # struct_al =
    struct_al = {u'lattice': {u'a': 2.8765103857767107, u'c': 2.87651039,
                              u'b': 2.876510387621576,
                              u'matrix': [[2.49113107, 0.0, 1.43825519],
                                          [0.83037702, 2.3486609, 1.43825519],
                                          [0.0, 0.0, 2.87651039]],
                              u'volume': 16.82995067829534,
                              u'beta': 60.000000066431895,
                              u'gamma': 60.0000001318988,
                              u'alpha': 60.00000008764776}, u'sites': [
        {u'abc': [0.0, 0.0, 0.0], u'xyz': [0.0, 0.0, 0.0], u'label': u'Al',
         u'species': [{u'occu': 1, u'element': u'Al'}]}],
                 u'@class': 'Structure', u'@module': 'pymatgen.core.structure'}
    struct_fe2o3 = {
        u'lattice': {u'a': 5.08889296704748, u'c': 5.088892940235904,
                     u'b': 5.433876369359943,
                     u'matrix': [[-1.63726303, -4.20120808, 2.3592482],
                                 [1.60869519, 4.18168075, 3.07451757],
                                 [-2.87170368, 4.20120808, 0.0]],
                     u'volume': 102.51695391433915,
                     u'beta': 119.99999825324359, u'gamma': 117.92116960555998,
                     u'alpha': 62.078837305387026}, u'sites': [
            {u'abc': [0.35372165, 0.06116494, 0.64627835],
             u'xyz': [-2.336659551680289, 1.484863824071041,
                      1.0225698487615258], u'label': u'Fe',
             u'species': [{u'occu': 1, u'element': u'Fe'}]},
            {u'abc': [0.14627835, 0.43883506, 0.85372165],
             u'xyz': [-1.9851798883197114, 4.807184630928958,
                      1.6943130362384742], u'label': u'Fe',
             u'species': [{u'occu': 1, u'element': u'Fe'}]},
            {u'abc': [0.64628007, 0.93884022, 0.35371993],
             u'xyz': [-0.5636015442146127, 2.6968240512458337,
                      4.41121584365604], u'label': u'Fe',
             u'species': [{u'occu': 1, u'element': u'Fe'}]},
            {u'abc': [0.85371993, 0.56115978, 0.14628007],
             u'xyz': [-0.9151020557853873, -0.6255110062458334,
                      3.7394328113439608], u'label': u'Fe',
             u'species': [{u'occu': 1, u'element': u'Fe'}]},
            {u'abc': [0.55410064, 0.75, 0.05410064],
             u'xyz': [0.14395189275030562, 1.0356565224999998,
                      3.613149115038848], u'label': u'O',
             u'species': [{u'occu': 1, u'element': u'O'}]},
            {u'abc': [0.25, 0.75, 0.44589936],
             u'xyz': [-0.4832851980216448, 3.959274536598829, 2.8957002275],
             u'label': u'O', u'species': [{u'occu': 1, u'element': u'O'}]},
            {u'abc': [0.05402801, 0.25, 0.25],
             u'xyz': [-0.4042101858574703, 1.8687392953416793,
                      0.896094877842082], u'label': u'O',
             u'species': [{u'occu': 1, u'element': u'O'}]},
            {u'abc': [0.75, 0.25, 0.55402801],
             u'xyz': [-2.416777750140077, 0.22210107965832115, 2.5380655425],
             u'label': u'O', u'species': [{u'occu': 1, u'element': u'O'}]},
            {u'abc': [0.44597199, 0.25, 0.94597199],
             u'xyz': [-3.0445508990024526, 3.1460242275, 1.820788007157918],
             u'label': u'O', u'species': [{u'occu': 1, u'element': u'O'}]},
            {u'abc': [0.94589936, 0.75, 0.75],
             u'xyz': [-2.495942419728661, 2.3132465884011713,
                      4.537499539961152], u'label': u'O',
             u'species': [{u'occu': 1, u'element': u'O'}]}],
        u'@class': 'Structure', u'@module': 'pymatgen.core.structure'}

    lp = LaunchPad.auto_load()

    for s in [struct_si, struct_al, struct_fe2o3]:
        wf = wf_band_structure(Structure.from_dict(s))
        lp.add_wf(wf)
