TIES
==============================
[//]: # (Badges)
[![CI](https://github.com/UCL-CCS/TIES20/actions/workflows/CI.yml/badge.svg)]
[![codecov](https://codecov.io/gh/REPLACE_WITH_OWNER_ACCOUNT/ties/branch/main/graph/badge.svg)](https://codecov.io/gh/REPLACE_WITH_OWNER_ACCOUNT/ties/branch/main)

Package for relative binding free energy calculations. See: 

TIES 20: Relative Binding Free Energy with a Flexible Superimposition Algorithm and Partial Ring Morphing
Mateusz K. Bieniek, Agastya P. Bhati, Shunzhou Wan, and Peter V. Coveney
Journal of Chemical Theory and Computation 2021 17 (2), 1250-1265
DOI: [10.1021/acs.jctc.0c01179](https://doi.org/10.1021/acs.jctc.0c01179)


### Install TIES

The easiest way to install TIES is with conda:

`conda install ties`

You can use your own charges with the .mol2 format:

`ties -l l1.mol2 l2.mol2 l3.mol2 --ligand-net-charge -2`

A PDB will automatically assume the need for BCC charges:

`ties -l left_coor.pdb right_coor.pdb -p protein.pdb -nc 0`

For help use `ties -h`

### Install TIES (dev)

Clone the repository:

`git clone https://github.com/UCL-CCS/TIES20.git`

It is simplest to use the same conda environment where `ambertools` was installed. Once in the right environment, go to the cloned directory and type:

```
conda env create -f environment.yml 
conda activate ties
pip install .
```

This should install the dependancies and make 
ties available in the environment. 

Please see the `examples` directory and `test cases` for more.

### Local ambertools and export the AMBERHOME variable

Conda can be used to easily install ambertools:
`conda install ambertools -c conda-forge`

Note that `sqm` in ambertools might need `libgfortran`. This can be installed with `apt-get install libgfortran3`. However, the conda version `conda install libgfortran -c conda-forge` should also work. 

This should set the variable AMBERHOME which will be picked up on the fly. 

### Copyright

Copyright (c) 2025, UCL CCS 


#### Acknowledgements
 
Project based on the 
[Computational Molecular Science Python Cookiecutter](https://github.com/molssi/cookiecutter-cms) version 1.10.