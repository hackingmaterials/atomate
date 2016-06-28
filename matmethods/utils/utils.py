# coding: utf-8

from __future__ import division, print_function, unicode_literals, absolute_import

import logging
import sys
import six

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
