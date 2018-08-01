from pymatgen.io.qchem.outputs import QCOutput
from atomate.qchem.drones import QChemDrone

#out = QCOutput("mol.qout.gz")
#print(out.data)


drone = QChemDrone()
doc = drone.assimilate("/scratch2/scratchdirs/sblau/fragment_wf/block_2018-07-04-03-30-43-540608/launcher_2018-07-06-17-03-13-189067", "mol.qin.gz", "mol.qout.gz", False)
print(doc) 
