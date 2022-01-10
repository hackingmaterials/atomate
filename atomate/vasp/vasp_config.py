__author__ = "Anubhav Jain <ajain@lbl.gov>"

import warnings

warnings.warn("vasp_config renamed to config.")
from .config import *  # noqa: F401, F403, E402
