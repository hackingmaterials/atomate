from custodian import Custodian
from custodian.qchem.jobs import QCJob
from custodian.qchem.handlers import QChemErrorHandler

my_input = "test.qin"
my_output = "test.qout"

myjob = QCJob.opt_with_frequency_flattener(qchem_command="qchem -slurm",multimode="openmp",input_file=my_input,output_file=my_output,max_iterations=10,max_molecule_perturb_scale=0.3,max_cores=12)
myhandler = QChemErrorHandler(input_file=my_input,output_file=my_output)

c = Custodian([myhandler],myjob,max_errors_per_job=10,max_errors=10)

c.run()

