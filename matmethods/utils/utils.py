# coding: utf-8

from __future__ import division, print_function, unicode_literals, \
    absolute_import

import six
import logging
import sys
import os
import shutil
import paramiko

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


def get_ssh_connection(username, host, pkey_file):
    """
    Setup ssh connection using paramiko and return the channel
    """
    try:
        privatekeyfile = os.path.expanduser(pkey_file)
    except:
        print("the private ket file {} doesnot exist".format(pkey_file))
        return None
    mykey = paramiko.RSAKey.from_private_key_file(privatekeyfile)
    ssh = paramiko.SSHClient()
    ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    try:
        ssh.connect(host, username=username, pkey=mykey)
    except paramiko.SSHException:
        print("Connection Error: host: {}, username: {}".format(host, username))
        return None
    return ssh


def list_dir(ldir, filesystem=None, pkey_file='~/.ssh/id_rsa'):
    """
    Wrapper of getting the directory listing from either the local or
    remote filesystem.

    Args:
        ldir (string): full path to the directory
        filesystem (string): remote filesystem, e.g. username@remote_host
        pkey_file (string): path to the private key file.
            Note: passwordless ssh login must be setup

    Returns:
        list of full paths to the files
    """
    if filesystem:
        tokens = filesystem.split('@')
        username = tokens[0]
        host = tokens[1]
        ssh = get_ssh_connection(username, host, pkey_file)
        if ssh:
            stdin, stdout, stderr = ssh.exec_command("ls "+ldir)
            return [os.path.join(ldir, l.split('\n')[0]) for l in stdout]
        else:
            print("paramiko connection error. Make sure that passwordless ssh login is setup and "
                  "yor private key is in standard location '~/.ssh/id_rsa'")
            return []
    else:
        return [os.path.join(ldir, f) for f in os.listdir(ldir)]


def copy_file(source, dest, filesystem=None, pkey_file='~/.ssh/id_rsa'):
    """
    Wrapper for copying from source to destination. The source can be
    a remote filesystem

    Args:
        source (string): source full path
        dest (string): destination file full path
        filesystem (string): remote filesystem, e.g. username@remote_host
        pkey_file (string): path to the private key file.
            Note: passwordless ssh login must be setup
    """
    if filesystem:
        tokens = filesystem.split('@')
        username = tokens[0]
        host = tokens[1]
        ssh = get_ssh_connection(username, host, pkey_file)
        if ssh:
            sftp = ssh.open_sftp()
            sftp.get(source, dest)
        else:
            print("paramiko connection error. Make sure that passwordless ssh login is setup and "
                  "yor private key is in standard location '~/.ssh/id_rsa'")
            raise IOError
    else:
        shutil.copy2(source, dest)


