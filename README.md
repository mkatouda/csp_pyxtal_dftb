# csp_pyxtal_dftb
Crytal structure prediction using PyXtal and DFTB+

##  Requirements 
1.  [Python]             (https://www.anaconda.com/download/)
2.  [PyXtal]             (https://github.com/qzhu2017/PyXtal)
3.  [DFTB+]              (https://dftbplus.org/)
4.  [rdkit]              (https://anaconda.org/rdkit/rdkit)
5.  [ASE]                (https://pypi.org/project/ase/)
6.  [pymatgen]           (https://anaconda.org/conda-forge/pymatgen)
7.  [numpy]              (https://anaconda.org/conda-forge/numpy)
8.  [scipy]              (https://anaconda.org/conda-forge/scipy)
9.  [matplotlib]         (https://anaconda.org/conda-forge/matplotlib)
10. [pandas]             (https://anaconda.org/conda-forge/pandas)
11. [acpype]             (https://anaconda.org/conda-forge/acpype)
12. [zipp]               (https://anaconda.org/conda-forge/zipp)
13. [llvmlite]           (https://anaconda.org/conda-forge/llvmlite)
14. [py3Dmol]            (https://anaconda.org/conda-forge/py3Dmol)
15. [numba]              (https://anaconda.org/conda-forge/numba)
16. [importlib-metadata] (https://anaconda.org/conda-forge/importlib-metadata)

## Optional
1.  [CP2K]               (https://www.cp2k.org/)
2.  [QuantumEspressso]   (https://www.quantum-espresso.org/)
3.  [Critic2]            (https://aoterodelaroza.github.io/critic2/)

## Installation of python envriromntent for csp_pyxtal_dftb
1. Install anaconda3 or miniconda3

2. Creaate python 3.10 virtual environment
'''
   $ conda create -n py310-molcryspred python=3.10
   $ conda activate py310-molcryspred  
   $ conda install -c conda-forge pyyaml openbabel rdkit ambertools lammps dftbplus qe cp2k ase pymatgen sssp  
   $ pip install pyxtal  
   $ pip install git+https://github.com/shirtsgroup/InterMol.git
   $ pip install git+https://github.com/mkatouda/csp_pyxtal_dftb.git
'''


## How to run csp_pyxtal sample job

### Crystal prediction of benzene
1. Move working dirctory
   $ cd ./tests/benzene  

2. Run python script csp_pyxtal_dftb.py with input yaml file (benzene.yml)
   $ bash ./run.sh
