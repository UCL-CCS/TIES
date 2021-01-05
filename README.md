# TIES 20
Topology Superimposition based on joint graph traversal. 

## How to install ambertools and export the AMBERHOME variable

Conda can be used to easily install ambertools:
`conda install ambertools -c conda-forge`

Note that `sqm` in ambertools might need `libgfortran`. This can be installed with `apt-get install libgfortran3`. However, the conda version `conda install libgfortran` should also work. 

This should set the variable AMBERHOME which will be picked up on the fly. 

# Install TIES 20

Clone the repository:

`git clone https://github.com/UCL-CCS/TIES20.git`

It is simplest to use the same conda environment where `ambertools` was installed. Once in the right environment, go to the cloned directory and type:

`pip install .` 

This should install the dependancies and make 
ties available in the environment. For help use: 

`ties -h`

I recommend renaming the molecules first manually:

`ties rename -l left_coor.pdb right_coor.pdb`

Then in the simplest case use the BCC charges. In the following example the net charge is set to -2:

`ties create -l left_coor.pdb right_coor.pdb -p protein.pdb -nc -2`

The protein is not necessary:
`ties create -l left_coor.pdb right_coor.pdb -nc -2`

You can use your own charges with the .mol2 format:
`ties create -l l1.mol2 l2.mol2 l3.mol2 -nc -2`

Please see the `examples` folder with more cases. 
