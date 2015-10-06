__author__ = 'Anubhav Jain <ajain@lbl.gov>'

def env_chk(val, fw_spec, strict=True):
    """
    env_chk() is a way to set different values for a property depending on the worker machine. For example,
    you might have slightly different executable names or scratch directories on different machines.

    env_chk() works using the principles of the FWorker env in FireWorks. For more details, see:
    https://pythonhosted.org/FireWorks/worker_tutorial.html

    This helper method translates string values that look like this:
    ">>ENV_KEY<<"
    to the contents of:
    fw_spec["_fw_env"][ENV_KEY]

    Since the latter can be set differently for each FireWorker, one can use this method to
    translate a single value into multiple possiblities, thus achieving different behavior
    on different machines.

    :param val: any value, with ">><<" notation reserved for special env lookup values
    :param fw_spec: fw_spec where one can find the _fw_env keys
    :param strict: (bool) if True, errors if env value cannot be found
    :return:
    """

    if isinstance(val, basestring) and val.startswith(">>") and val.endswith("<<"):
        if strict:
            return fw_spec['_fw_env'][val[2:-2]]
        return fw_spec.get('_fw_env', {}).get(val[2:-2])
    return val
