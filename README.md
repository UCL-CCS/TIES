# BAC2020
Topology Superimposition based on joint graph traversal

Clone the repository:

`git clone https://github.com/UCL-CCS/TIES20.git`

Change the directory. Several packages will be installed as dependancies. 
You can also install them manually. See "setup.py" which contains the list. 
To install ties: 

`python install . # will be tied to the github`

To get help use: 

`ties.py -h`

I recommend renaming the molecules first

`ties.py rename -l left_coor.pdb -r right_coor.pdb`

Then in the simplest case the charge

`ties.py create -l left_coor.pdb -r right_coor.pdb -p protein.pdb -nc -2`