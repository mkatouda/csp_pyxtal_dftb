# pull Ubuntu 22.04
FROM ubuntu:22.04

# copy csp_pyxtal_dft and critic2
COPY csp_pyxtal_dftb /csp_pyxtal_dftb
COPY critic2 /critic2

ENV PATH=/opt/miniconda3/bin:/critic2/bin:$PATH

# install miniconda3 and libraries
RUN apt-get update \
    && apt-get upgrade -y \
    && apt-get install -y curl wget git zip \
    && curl -sSL https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -o /tmp/miniconda.sh \
    && bash /tmp/miniconda.sh -bfp /opt/miniconda3 \
    && rm -rf /tmp/miniconda.sh \
    && conda config --add channels conda-forge \
    && conda update conda -y \
    && conda install -c conda-forge mamba -y \
    && mamba install -c conda-forge python=3.10 -y \
    && mamba install -c conda-forge pyyaml openbabel rdkit ambertools lammps dftbplus qe cp2k ase pymatgen sssp -y \
    && pip install pyxtal \
    && pip install git+https://github.com/shirtsgroup/InterMol.git \
    && cd csp_pyxtal_dftb && pip install . && cd ../
