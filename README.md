# csp_pyxtal_dftb
Crytal structure prediction using PyXtal and DFTB+

##  Requirements 
1.  [Python]             (https://www.anaconda.com/download/)
2.  [PyXtal]             (https://github.com/qzhu2017/PyXtal)
3.  [rdkit 2020.09.5]    (https://anaconda.org/rdkit/rdkit)
4.  [ASE]                (https://pypi.org/project/ase/)
5.  [pymatgen]           (https://anaconda.org/conda-forge/pymatgen)
6.  [numpy]              (https://anaconda.org/conda-forge/numpy)
7.  [scipy]              (https://anaconda.org/conda-forge/scipy)
8.  [matplotlib]         (https://anaconda.org/conda-forge/matplotlib)
9.  [pandas]             (https://anaconda.org/conda-forge/pandas)
10. [acpype]             (https://anaconda.org/conda-forge/acpype)
11. [zipp]               (https://anaconda.org/conda-forge/zipp)
12. [llvmlite]           (https://anaconda.org/conda-forge/llvmlite)
13. [py3Dmol]            (https://anaconda.org/conda-forge/py3Dmol)
14. [numba]              (https://anaconda.org/conda-forge/numba)
15. [importlib-metadata] (https://anaconda.org/conda-forge/importlib-metadata)

## Optional
1.  [xTB]                (https://anaconda.org/conda-forge/xtb)
2.  [xtb-python]         (https://anaconda.org/conda-forge/xtb-python)


## Installation of python envriromntent for csp_pyxtal_dftb
1. Install anaconda3 or miniconda3
2. Creaate python 3.8 virtual environment
   $ conda create -n py38-csp python=3.8  
   $ conda activate py38-csp  
   $ conda install numpy scipy matplotlib pandas xtb xtb-python acpype pymatgen zipp llvmlite py3Dmol numba importlib-metadata  
   $ conda install -c conda-forge rdkit==2020.09.5  
   $ pip install --upgrade --user ase  
   $ pip install pyxtal  
3. Install DFTB+ following DFTB+ Recipes> Introduction (https://dftbplus-recipes.readthedocs.io/en/latest/introduction.html#before-you-start)


## How to run LAQA confopt job

### Crystal prediction of benzene
1. Move working dirctory
   $ cd ./tests/benzene  
2. Set environmental variables for DFTB+ code (bash or ksh cases).
   $ export ASE_DFTB_COMMAND="/path/to/dftb+/dftb+ > PREFIX.out"  
   $ export DFTB_PREFIX=/path/to/slakos/origin/3ob-3-1  
   $ export OMP_NUM_THREADS=8  
   $ export OMP_STACKSIZE=2GB  
   
3. Run python script csp_pyxtal_dftb.py with input yaml file
   $ python /path/to/csp_pyxtal_dftb/csp_pyxtal_dftb/csp_pyxtal_dftb.py benzene.yaml  
