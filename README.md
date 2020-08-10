# TIES 20
Topology Superimposition based on joint graph traversal. 

Clone the repository:

`git clone https://github.com/UCL-CCS/TIES20.git`

To install ties go to the cloned directory and type:   

`python install .` 

This should install the dependancies and make 
ties available in the environment. For help use: 

`ties -h`

I recommend renaming the molecules first manually:

`ties rename -l left_coor.pdb -r right_coor.pdb`

Then in the simplest case use the BCC charges:

`ties create -l left_coor.pdb -r right_coor.pdb -p protein.pdb -nc -2`


## How to install ambertools and export the AMBERHOME variable

Conda can be used to easily install ambertools:
`conda install ambertools`

Note that `sqm` in ambertools might need `libgfortran`. This can be installed with `apt-get install libgfortran3`. However, the conda version `conda install libgfortran` should also work. 

This should set the variable AMBERHOME which will be picked up on the fly. 
