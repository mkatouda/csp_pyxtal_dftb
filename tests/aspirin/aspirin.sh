#!/bin/bash

LANG=C
NUM_THREADS=1

export MKL_NUM_THREADS=${NUM_THREADS}
export OMP_NUM_THREADS=${NUM_THREADS},1
export OMP_STACKSIZE=4G

. /home/app/intel/intel2020_up1/compilers_and_libraries_2020.1.217/linux/bin/compilervars.sh intel64
#/home/center/opt/x86_64/cores/intel/compilers_and_libraries_2019.5.281/linux/bin/compilervars.sh intel64
. ${HOME}/miniconda3/etc/profile.d/conda.sh

export CSP_PYXTAL_DFTB_PREFIX=${HOME}/csp_pyxtal_dftb/csp_pyxtal_dftb
export ASE_DFTB_COMMAND="${HOME}/bin/dftbplus/21.1/intel2021.1_serial/bin/dftb+ > PREFIX.out"
export DFTB_PREFIX=${HOME}/dftbplus-21.1/external/slakos/origin/3ob-3-1

conda activate py38-csp

in_file=aspirin.yaml
out_file=aspirin.log

time python ${CSP_PYXTAL_DFTB_PREFIX}/csp_pyxtal_dftb.py -c ${in_file} >& ${out_file}

conda deactivate
