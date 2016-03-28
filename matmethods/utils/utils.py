# coding: utf-8

from __future__ import division, print_function, unicode_literals, \
    absolute_import

import logging
import sys

import six

__author__ = 'Anubhav Jain, Kiran Mathew'
__email__ = 'ajain@lbl.gov, kmathew@lbl.gov'


def env_chk(val, fw_spec, strict=True):
    """
    env_chk() is a way to set different values for a property depending
    on the worker machine. For example, you might have slightly different
    executable names or scratch directories on different machines.

    env_chk() works using the principles of the FWorker env in FireWorks.
    For more details, see:
    https://pythonhosted.org/FireWorks/worker_tutorial.html

    This helper method translates string values that look like this:
    ">>ENV_KEY<<"
    to the contents of:
    fw_spec["_fw_env"][ENV_KEY]

    Since the latter can be set differently for each FireWorker, one can
    use this method to translate a single value into multiple possibilities,
    thus achieving different behavior on different machines.

    Args:
        val: any value, with ">><<" notation reserved for special env lookup
            values
        fw_spec: fw_spec where one can find the _fw_env keys
        strict(bool): if True, errors if env value cannot be found
    """

    if isinstance(val, six.string_types) and val.startswith(
            ">>") and val.endswith("<<"):
        if strict:
            return fw_spec['_fw_env'][val[2:-2]]
        return fw_spec.get('_fw_env', {}).get(val[2:-2])
    return val


def get_logger(name, level=logging.DEBUG,
               format='%(asctime)s %(levelname)s %(name)s %(message)s',
               stream=sys.stdout):
    logger = logging.getLogger(name)
    logger.setLevel(level)
    formatter = logging.Formatter(format)
    sh = logging.StreamHandler(stream=stream)
    sh.setFormatter(formatter)
    logger.addHandler(sh)
    return logger
