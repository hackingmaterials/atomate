# coding: utf-8


import glob
import os
import shutil

"""
This module defines the wrapper class for remote file io using paramiko.
"""

__author__ = 'Kiran Mathew'
__credits__ = 'Anubhav Jain <ajain@lbl.gov>'
__email__ = 'kmathew@lbl.gov'


class FileClient(object):
    """
    A client for performing many file operations while being agnostic
    of whether those operations are happening locally or via SSH
    """

    def __init__(self, filesystem=None, private_key="~/.ssh/id_rsa"):
        """
        Args:
            filesystem (str): remote filesystem, e.g. username@remote_host.
                If None, use local
            private_key (str): path to the private key file (for remote
                connections only). Note: passwordless ssh login must be setup
        """
        self.ssh = None

        if filesystem:
            if '@' in filesystem:
                username, host = filesystem.split('@', 1)
            else:
                username = None  # paramiko sets default username
                host = filesystem

            self.ssh = FileClient.get_ssh_connection(username, host, private_key)
            self.sftp = self.ssh.open_sftp()

    @staticmethod
    def get_ssh_connection(username, host, private_key):
        """
        Connect to the remote host via paramiko using the private key.
        If the host key is not present it will be added automatically.

        Args:
            username (str):
            host (str):

            private_key (str):  path to private key file

        Returns:
            SSHClient

        """
        import paramiko
        private_key = os.path.expanduser(private_key)
        if not os.path.exists(private_key):
            raise ValueError("Cannot locate private key file: {}".format(private_key))

        ssh = paramiko.SSHClient()
        ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
        return ssh.connect(host, username=username, key_filename=private_key)

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
            path (str): path to get absolute string of
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
            path (str): path to glob
        """
        if not self.ssh:
            return glob.glob(path)
        else:
            command = ". ./.bashrc; for i in $(ls {}); do readlink -f $i; done".format(path)
            stdin, stdout, stderr = self.ssh.exec_command(command)
            return [l.split('\n')[0] for l in stdout]
