
from tqdm import tqdm
import math
from atomate.utils.utils import get_logger, get_database

logger = get_logger(__name__)

__author__ = 'Anubhav Jain <ajain@lbl.gov>'

"""
Estimates band gap according to equations in:
Ravindra, N. M., Ganapathy, P. & Choi, J. Energy gap-refractive index relations in semiconductors 
- An overview. Infrared Phys. Technol. 50, 21-29 (2007).

Often these estimates are comparable to HSE gap estimates even if based on GGA dielectric constants 
for gaps approximately ~2 eV and higher. Smaller gaps are less accurate. Preliminary testing 
(A. Jain) suggests use of Reddy-Anjaneyulu relation as most reliable.
"""


class BandgapEstimationBuilder:
    def __init__(self, materials_write):
        """
        Starting with an existing materials collection with dielectric constant data, adds
        estimated band gaps that may be more accurate than typical GGA calculations.

        Run the "DielectricBuilder" before running this builder.

        Args:
            materials_write: mongodb collection for materials (write access needed)
        """
        self._materials = materials_write

    def run(self):
        logger.info("{} starting...".format(self.__class__.__name__))
        q = {"dielectric.epsilon_static_avg": {"$gt": 0}, "bandgap_estimation": {"$exists": False}}

        for m in tqdm(self._materials.find(q, projection=["material_id", "dielectric"])):
            try:
                eps = m["dielectric"]["epsilon_static_avg"]  # electronic portion of eps ("eps_static") approximates eps_inf
                n = math.sqrt(eps)  # sqrt(eps_inf) to get refractive index
                d = {}
                d["gap_moss"] = (95 / n**4) if n > 0 else None
                d["gap_gupta-ravindra"] = (4.16-n)/0.85 if n <= 4.16 else None
                d["gap_reddy-anjaneyulu"] = 36.3/math.exp(n)
                d["gap_reddy-ahamed"] = 154/n**4+0.365 if n > 0 else None
                d["gap_herve_vandamme"] = 13.47/math.sqrt(n**2-1)-3.47 if n > 1 else None

                d = {"bandgap_estimation": d}
                self._materials.update_one({"material_id": m["material_id"]}, {"$set": d})

            except:
                import traceback
                logger.exception(traceback.format_exc())

        logger.info("{} finished.".format(self.__class__.__name__))

    def reset(self):
        logger.info("Resetting {} starting!".format(self.__class__.__name__))
        self._materials.update_many({}, {"$unset": {"bandgap_estimation": 1}})
        logger.info("Resetting {} finished!".format(self.__class__.__name__))

    @staticmethod
    def from_file(db_file, m="materials", **kwargs):
        """
        Get builder using only a db file.

        Args:
            db_file: (str) path to db file
            m: (str) name of "materials" collection
            **kwargs: other parameters to feed into the builder, e.g. mapi_key

        Returns:
            BandgapEstimationBuilder
        """
        db_write = get_database(db_file, admin=True)
        return BandgapEstimationBuilder(db_write[m], **kwargs)