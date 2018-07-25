#!/bin/bash -l
#SBATCH -q debug
#SBATCH -N 1
#SBATCH -t 00:30:00  
#SBATCH -J new_def    
#SBATCH -A jcesr
#SBATCH -C haswell 

#set QC run time environment
#export QC=/global/project/projectdirs/qchem/qc51beta
#export QCAUX=$QC/qcaux
#. qcenv.sh

module load qchem

#load intel mpi
module load impi    
#export OMP_PROC_BIND=true     which one is better?
export OMP_PROC_BIND=spread
export OMP_PLACES=threads
export I_MPI_FABRICS=ofi
export I_MPI_OFI_PROVIDER=gni
export I_MPI_OFI_LIBRARY=/usr/common/software/libfabric/1.5.0/gnu/lib/libfabric.so
export I_MPI_PMI_LIBRARY=/usr/lib64/slurmpmi/libpmi.so


#invoke QChem
#export QCSCRATCH=/global/homes/z/zgan/proj/qc50/scratch
export QCSCRATCH=/global/cscratch1/sd/sblau/

qchem -slurm -nt 12 julian.qin
