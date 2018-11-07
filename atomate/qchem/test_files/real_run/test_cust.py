from custodian import Custodian
from custodian.qchem.jobs import QCJob
from custodian.qchem.handlers import QChemErrorHandler

my_input = "mol.qin"
my_output = "mol.qout"

myjob = QCJob(qchem_command="qchem -slurm",multimode="openmp",input_file=my_input,output_file=my_output,max_cores=12)
myhandler = QChemErrorHandler(input_file=my_input,output_file=my_output)

c = Custodian([myhandler],[myjob],max_errors_per_job=10,max_errors=10)

c.run()

