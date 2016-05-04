# coding: utf-8

from __future__ import division, print_function, unicode_literals, absolute_import

import glob
import os
import shutil

"""
This module defines the wrapper class for remote file io using paramiko.
"""

__author__ = 'Kiran Mathew'
__credits__ = 'Anubhav Jain <ajain@lbl.gov>'
__email__ = 'kmathew@lbl.gov'


# TODO: properly document this (i.e., all methods)

class FileClient(object):
    """
    A client that allows performing many file operations while being agnostic of whether
    those operations are happening locally or via SSH
    """

    def __init__(self, filesystem=None, pkey_file="~/.ssh/id_rsa"):
        """
        Args:
            filesystem (str): remote filesystem, e.g. username@remote_host. If None, use local
            pkey_file (str): path to the private key file (for remote connections only)
                Note: passwordless ssh login must be setup
        """
        self.ssh = None
        if filesystem:
            if '@' in filesystem:
                username, host = filesystem.split('@')
            else:
                username = None  # paramiko sets default username
                host = filesystem

            self.ssh = FileClient.get_ssh_connection(username, host, pkey_file)
            self.sftp = self.ssh.open_sftp()

    @staticmethod
    def get_ssh_connection(username, host, pkey_file):
        """
        Args:
            username (str):
            host (str):
            pkey_file (str):  path to private key file

        Returns:
            SSHClient

        """
        import paramiko
        pkey_file = os.path.expanduser(pkey_file)
        if not os.path.exists(pkey_file):
            raise ValueError("Cannot locate private key file: {}".format(pkey_file))

        ssh = paramiko.SSHClient()
        ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
        return ssh.connect(host, username=username, key_filename=pkey_file)

    @staticmethod
    def exists(sftp, path):
        """
        os.path.exists() for paramiko's SCP object

        Args:
            sftp (SFTPClient):
            path (str): path to check existence of
        """
        try:
            sftp.stat(path)
        except IOError as e:
            if e[0] == 2:
                return False
            raise
        else:
            return True

    def listdir(self, ldir):
        """
        Get the directory listing from either the local or remote filesystem.

        Args:
            ldir (str): full path to the directory

        Returns:
            iterator of filenames
        """

        # TODO: this pattern of "if self.ssh: self.ssh.X() else os.X() could be generalized beyond X()=listdir() in much more general code than this

        if not self.ssh:
            return os.listdir(ldir)
        else:
            return self.sftp.listdir()

    def copy(self, src, dest):
        """
        Copy from source to destination.

        Args:
            src (str): source full path
            dest (str): destination file full path

        """
        if not self.ssh:
            shutil.copy2(src, dest)

        else:
            if os.path.isdir(src):
                if not FileClient.exists(self.sftp, dest):
                    self.sftp.mkdir(dest)
                for f in os.listdir(src):
                    if os.path.isfile(os.path.join(src, f)):
                        self.sftp.put(os.path.join(src, f), os.path.join(dest, f))
            else:
                self.sftp.put(src, os.path.join(dest, os.path.basename(src)))

    def abspath(self, path):
        """
        return the absolute path

        Args:
            path (str):
        """
        if not self.ssh:
            return os.path.abspath(path)

        else:
            command = ". ./.bashrc; readlink -f {}".format(path)
            stdin, stdout, stderr = self.ssh.exec_command(command)
            full_path = [l.split('\n')[0] for l in stdout]
            return full_path[0]

    def glob(self, path):
        """
        return the glob

        Args:
            path (str):
        """
        if not self.ssh:
            return glob.glob(path)

        else:
            command = ". ./.bashrc; for i in $(ls {}); do readlink -f $i; done".format(path)
            stdin, stdout, stderr = self.ssh.exec_command(command)
            return [l.split('\n')[0] for l in stdout]
