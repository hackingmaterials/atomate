# coding: utf-8

from __future__ import division, print_function, unicode_literals, absolute_import

import glob
import os
import shutil

"""
This module defines the wrapper class for remote file io using paramiko.
"""

__author__ = 'Kiran Mathew'
__email__ = 'kmathew@lbl.gov'


class MMos(object):
    """
    paramiko wrapper
    """

    def __init__(self, filesystem=None, pkey_file="~/.ssh/id_rsa"):
        """
        filesystem (string): remote filesystem, e.g. username@remote_host
        pkey_file (string): path to the private key file.
            Note: passwordless ssh login must be setup
        """
        self.ssh = None
        username = None
        host = None
        if filesystem:
            tokens = filesystem.split('@')
            username = tokens[0]
            host = tokens[1]
        if username and host:
            self.ssh = self._get_ssh_connection(username, host, pkey_file)

    def _get_ssh_connection(self, username, host, pkey_file):
        """
        Setup ssh connection using paramiko and return the channel
        """
        import paramiko
        privatekeyfile = os.path.expanduser(pkey_file)
        if not os.path.exists(privatekeyfile):
            possible_keys = ["~/.ssh/id_rsa", "~/.ssh/id_dsa", "/etc/ssh/id_rsa", "/etc/ssh/id_dsa"]
            for key in possible_keys:
                if os.path.exists(os.path.expanduser(key)):
                    privatekeyfile = os.path.expanduser(key)
                    break
        tokens = privatekeyfile.split("id_")
        try:
            if tokens[1] == "rsa":
                mykey = paramiko.RSAKey.from_private_key_file(privatekeyfile)
            elif tokens[1] == "dsa":
                mykey = paramiko.DSSKey.from_private_key_file(privatekeyfile)
            else:
                print("Unknown private key format. Must be either rsa(preferred) or dsa")
        except:
            print("Found the private key file {}, but not able to load".format(pkey_file))
            return None
        ssh = paramiko.SSHClient()
        ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
        try:
            ssh.connect(host, username=username, pkey=mykey)
        except paramiko.SSHException:
            print("Connection Error: host: {}, username: {}".format(host, username))
            return None
        return ssh

    def listdir(self, ldir):
        """
        Get the directory listing from either the local or remote filesystem.

        Args:
            ldir (string): full path to the directory

        Returns:
            list of filenames
        """
        if self.ssh:
            try:
                #command = ". ./.bashrc; for i in {}/*; do readlink -f $i; done".format(ldir)
                command = ". ./.bashrc; ls {}".format(ldir)
                stdin, stdout, stderr = self.ssh.exec_command(command)
                return [l.split('\n')[0] for l in stdout]
            except:
                print("paramiko connection error. Make sure that passwordless ssh login is setup and "
                      "yor private key is in standard location. e.g. '~/.ssh/id_rsa'")
                return []
        else:
            return [f for f in os.listdir(ldir)]

    def copy(self, source, dest):
        """
        Copy from source to destination.

        Args:
            source (string): source full path
            dest (string): destination file full path

        """
        if self.ssh:
            try:
                command = ". ./.bashrc; readlink -f {}".format(source)
                stdin, stdout, stderr = self.ssh.exec_command(command)
                source_full_path = [l.split('\n')[0] for l in stdout]
                sftp = self.ssh.open_sftp()
                sftp.get(source_full_path[0], dest)
            except:
                print("paramiko connection error. Make sure that passwordless ssh login is setup and "
                      "yor private key is in standard location. e.g. '~/.ssh/id_rsa'")
                raise IOError
        else:
            shutil.copy2(source, dest)

    def abspath(self, path):
        """
        return the absolute path
        """
        if self.ssh:
            try:
                command = ". ./.bashrc; readlink -f {}".format(path)
                stdin, stdout, stderr = self.ssh.exec_command(command)
                full_path = [l.split('\n')[0] for l in stdout]
                return full_path[0]
            except:
                print("paramiko connection error. Make sure that passwordless ssh login is setup and "
                      "yor private key is in standard location. e.g. '~/.ssh/id_rsa'")
                raise IOError
        else:
            return os.path.abspath(path)

    def glob(self, path):
        """
        return the glob
        """
        if self.ssh:
            try:
                command = ". ./.bashrc; for i in $(ls {}); do readlink -f $i; done".format(path)
                stdin, stdout, stderr = self.ssh.exec_command(command)
                return [l.split('\n')[0] for l in stdout]
            except:
                print("paramiko connection error. Make sure that passwordless ssh login is setup and "
                      "yor private key is in standard location. e.g. '~/.ssh/id_rsa'")
                raise IOError
        else:
            return glob.glob(path)