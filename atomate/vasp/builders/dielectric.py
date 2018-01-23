from tqdm import tqdm

from atomate.utils.utils import get_logger

import numpy as np

from atomate.utils.utils import get_database

logger = get_logger(__name__)

__author__ = 'Shyue Ping Ong <ongsp@uscd.edu>, Anubhav Jain <ajain@lbl.gov>'


class DielectricBuilder:

    def __init__(self, materials_write):
        """
        Starting with an existing materials collection, adds some averages and 
        eigenvalues for dielectric constants rather than just the tensor
        
        Args:
            materials_write: mongodb collection for materials (write access needed)
        """
        self._materials = materials_write

    def run(self):
        logger.info("EpsilonBuilder starting...")
        q = {"dielectric": {"$exists": True}, "dielectric.eps_ionic_avg": {"$exists": False}}

        for m in tqdm(self._materials.find(q, projection=["material_id", "dielectric"])):
            try:
                eps = m["dielectric"]
                d = {}
                eig_ionic = np.linalg.eig(eps["epsilon_ionic"])[0]
                eig_static = np.linalg.eig(eps["epsilon_static"])[0]

                d["dielectric.epsilon_ionic_avg"] = float(np.average(eig_ionic))
                d["dielectric.epsilon_static_avg"] = float(np.average(eig_static))
                d["dielectric.epsilon_avg"] = d["dielectric.epsilon_ionic_avg"] + \
                                              d["dielectric.epsilon_static_avg"]
                d["dielectric.has_neg_eps"] = bool(np.any(eig_ionic < -0.1) or
                                                   np.any(eig_static < -0.1))

                self._materials.update_one({"material_id": m["material_id"]}, {"$set": d})

            except:
                import traceback
                logger.exception(traceback.format_exc())

        logger.info("EpsilonBuilder finished processing.")

    def reset(self):
        logger.info("Resetting EpsilonBuilder")
        keys = ["dielectric.epsilon_ionic_avg",
                "dielectric.epsilon_static_avg",
                "dielectric.epsilon_avg",
                "dielectric.has_neg_eps"]

        self._materials.update_many({}, {"$unset": {k: "" for k in keys}})
        logger.info("Finished resetting EpsilonBuilder")

    @staticmethod
    def from_file(db_file, m="materials", **kwargs):
        """
        Get a MaterialsEhullBuilder using only a db file.

        Args:
            db_file: (str) path to db file
            m: (str) name of "materials" collection
            **kwargs: other parameters to feed into the builder, e.g. mapi_key

        Returns:
            DielectricBuilder
        """
        db_write = get_database(db_file, admin=True)
        return DielectricBuilder(db_write[m], **kwargs)
