# coding: utf-8

from __future__ import division, print_function, unicode_literals, absolute_import

import logging
import sys
import six

from pymatgen import Composition

__author__ = 'Anubhav Jain, Kiran Mathew'
__email__ = 'ajain@lbl.gov, kmathew@lbl.gov'


def env_chk(val, fw_spec, strict=True, default=None):
    """
    env_chk() is a way to set different values for a property depending
    on the worker machine. For example, you might have slightly different
    executable names or scratch directories on different machines.

    env_chk() works using the principles of the FWorker env in FireWorks.

    This helper method translates string values that look like this:
    ">>ENV_KEY<<"
    to the contents of:
    fw_spec["_fw_env"][ENV_KEY]

    The fw_spec["_fw_env"] is in turn set by the FWorker. For more details,
    see: https://pythonhosted.org/FireWorks/worker_tutorial.html

    Since the fw_env can be set differently for each FireWorker, one can
    use this method to translate a single value into multiple possibilities,
    thus achieving different behavior on different machines.

    Args:
        val: any value, with ">><<" notation reserved for special env lookup
            values
        fw_spec: (dict) fw_spec where one can find the _fw_env keys
        strict (bool): if True, errors if env value cannot be found
        default: if val is None or env cannot be found in non-strict mode,
                 return default
    """
    if val is None:
        return default

    if isinstance(val, six.string_types) and val.startswith(">>") and val.endswith("<<"):
        if strict:
            return fw_spec['_fw_env'][val[2:-2]]
        return fw_spec.get('_fw_env', {}).get(val[2:-2], default)
    return val


def get_calc_loc(target_name, calc_locs):
    """
    This is a helper method that helps you pick out a certain calculation
    from an array of calc_locs.

    There are three modes:
        - If you set target_name to a String, search for most recent calc_loc
            with matching name
        - Otherwise, return most recent calc_loc overall

    Args:
        target_name: (bool or str) If str, will search for calc_loc with
            matching name, else use most recent calc_loc
        calc_locs: (dict) The dictionary of all calc_locs

    Returns:
        (dict) dict with subkeys path, filesystem, and name
    """

    if isinstance(target_name, six.string_types):
        for doc in reversed(calc_locs):
            if doc["name"] == target_name:
                return doc
        raise ValueError("Could not find the target_name: {}".format(target_name))
    else:
        return calc_locs[-1]


def get_mongolike(d, key):
    if "." in key:
        i, j = key.split(".", 1)
        try:
            i = int(i)
        except:
            pass
        return get_mongolike(d[i], j)
    return d[key]


def get_logger(name, level=logging.DEBUG, format='%(asctime)s %(levelname)s %(name)s %(message)s',
               stream=sys.stdout):
    logger = logging.getLogger(name)
    logger.setLevel(level)
    formatter = logging.Formatter(format)
    sh = logging.StreamHandler(stream=stream)
    sh.setFormatter(formatter)
    logger.addHandler(sh)
    return logger


def get_meta_from_structure(structure):
    comp = structure.composition
    elsyms = sorted(set([e.symbol for e in comp.elements]))
    meta = {'nsites': len(structure),
            'elements': elsyms,
            'nelements': len(elsyms),
            'formula': comp.formula,
            'formula_reduced': comp.reduced_formula,
            'formula_reduced_abc': Composition(comp.reduced_formula)
            .alphabetical_formula,
            'formula_anonymous': comp.anonymized_formula,
            'chemsys': '-'.join(elsyms),
            'is_ordered': structure.is_ordered,
            'is_valid': structure.is_valid()}
    return meta


# TODO: move this code as soon as possible.
# It causes errors in particular Python installations - and even though there is a try-catch, it does not recover!!!
#   from matplotlib.backends import _macosx
# RuntimeError: Python is not installed as a framework. The Mac OS X backend will not be able to function correctly if Python is not installed as a framework. See the Python documentation for more information on installing Python as a framework on Mac OS X. Please either reinstall Python as a framework, or try one of the other backends. If you are Working with Matplotlib in a virtual enviroment see 'Working with Matplotlib in Virtual environments' in the Matplotlib FAQ

def plot_wf(wf, depth_factor=1.0, breadth_factor=2.0, labels_on=True, numerical_label=False,
            text_loc_factor=1.0, save_as=None, style='rD--', markersize=10, markerfacecolor='blue',
            fontsize=12,):
    """
    Generate a visual representation of the workflow. Useful for checking whether the firework
    connections are in order before launching the workflow.

    Args:
        wf (Workflow): workflow object.
        depth_factor (float): adjust this to stretch the plot in y direction.
        breadth_factor (float): adjust this to stretch the plot in x direction.
        labels_on (bool): whether to label the nodes or not. The default is to lable the nodes
            using the firework names.
        numerical_label (bool): set this to label the nodes using the firework ids.
        text_loc_factor (float): adjust the label location.
        save_as (str): save the figure to the given name.
        style (str): marker style.
        markersize (int): marker size.
        markerfacecolor (str): marker face color.
        fontsize (int): font size for the node label.
    """
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        plt = None

    keys = sorted(wf.links.keys(), reverse=True)

    # set (x,y) coordinates for each node in the workflow links
    points_map = {}
    points_map.update({keys[0]: (-0.5*breadth_factor, (keys[0]+1)*depth_factor)})
    for k in keys:
        if wf.links[k]:
            for i, j in enumerate(wf.links[k]):
                if not points_map.get(j, None):
                    points_map[j] = ((i-len(wf.links[k])/2.0)*breadth_factor, k*depth_factor)

    # connect the dots
    for k in keys:
        for i in wf.links[k]:
            plt.plot([points_map[k][0], points_map[i][0]], [points_map[k][1], points_map[i][1]],
                     style, markersize=markersize, markerfacecolor=markerfacecolor)
            if labels_on:
                label1 = wf.id_fw[k].name
                label2 = wf.id_fw[i].name
                if numerical_label:
                    label1 = str(k)
                    label2 = str(i)
                plt.text(points_map[k][0] * text_loc_factor, points_map[k][1] * text_loc_factor,
                         label1, fontsize=fontsize)
                plt.text(points_map[i][0] * text_loc_factor, points_map[i][1] * text_loc_factor,
                         label2, fontsize=fontsize)

    plt.axis('scaled')
    plt.axis('off')

    if save_as:
        plt.savefig(save_as)
