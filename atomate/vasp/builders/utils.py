"""
This class contains common functions for builders
"""

__author__ = 'Anubhav Jain <ajain@lbl.gov>'


def dbid_to_str(prefix, dbid):
    # converts int dbid to string (adds prefix)
    return "{}-{}".format(prefix, dbid)


def dbid_to_int(dbid):
    # converts string dbid to int (removes prefix)
    return int(dbid.split("-")[1])
