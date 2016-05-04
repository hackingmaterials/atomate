import unittest
import getpass

from matmethods.utils.fileio import FileClient


# TODO: Kiran please either activate this test or remove it ...

"""
class FileIOTests(unittest.TestCase):

    #@unittest.skipIf(not os.path.exists(os.path.expanduser("~/.ssh/id_rsa")) and not
    #os.path.exists(os.path.expanduser("~/.ssh/authorized_keys")),
    #                 "no '~/.ssh/id_rsa' private key file paramiko test skipped")
    @unittest.skip("paramiko test skipped")
    def test_remote_filesystem(self):
        username = getpass.getuser()
        host = "localhost"
        filesystem = "{}@{}".format(username, host)
        fileclient = FileClient(filesystem)
        self.assertIsNotNone(fileclient.ssh)

"""