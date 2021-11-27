#!/bin/bash
#PJM -L "rscunit=ito-a"
#PJM -L "rscgrp=ito-ss"
#PJM -L "vnode=1"
#PJM -L "vnode-core=36"
#PJM -L "elapse=%TIME%"
#PJM -j
#PJM -S

LANG=C

NUM_NODES=1 #${PJM_VNODES}
NUM_CORES=36
NUM_THREADS=36
NUM_PROCS=`expr ${NUM_NODES} "*" ${NUM_CORES} / ${NUM_THREADS}`
echo NUM_NODES=${NUM_NODES} NUM_CORES=${NUM_CORES} NUM_PROCS=${NUM_PROCS} NUM_THREADS=${NUM_THREADS}

export MKL_NUM_THREADS=${NUM_THREADS}
export OMP_NUM_THREADS=${NUM_THREADS},1
export OMP_STACKSIZE=64G
#ulimit -s unlimited

. /home/app/intel/intel2020_up1/compilers_and_libraries_2020.1.217/linux/bin/compilervars.sh intel64
. ${HOME}/miniconda3/etc/profile.d/conda.sh

export CSP_PYXTAL_DFTB_PREFIX=${HOME}/csp_pyxtal_dftb/csp_pyxtal_dftb
export ASE_DFTB_COMMAND="${HOME}/bin/dftbplus/21.1/intel2021.1_serial/bin/dftb+ > PREFIX.out"
export DFTB_PREFIX=${HOME}/dftbplus-21.1/external/slakos/origin/3ob-3-1

conda activate py38-cryspy

in_file=%IN_FILE%
out_file=%OUT_FILE%

time python ${CSP_PYXTAL_DFTB_PREFIX}/csp_pyxtal_dftb.py -c ${in_file} %RUNOPT% >& ${out_file}

conda deactivate
