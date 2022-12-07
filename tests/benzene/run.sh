#!/bin/bash

in_file=benzene.yml
out_file=benzene.log

LANG=C
NUM_THREADS=4

export MKL_NUM_THREADS=${NUM_THREADS}
export OMP_NUM_THREADS=${NUM_THREADS},1
export OMP_STACKSIZE=4G

csp_pyxtal -i ${in_file} >& ${out_file}
