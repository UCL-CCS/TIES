source ~/software/amber18install/amber.sh
# Just the first stage of providing the charges and .mol2 file

# todo - check if these are files with the right extension etc
LEFT=$1
RIGHT=$2

# update to eval rather than ``
CURRENT_DIR=`dirname $0`
cd $CURRENT_DIR
echo 'Current Dir' $CURRENT_DIR

# generate the mol2 topology
antechamber -i $LEFT -fi pdb -o left.mol2 -fo mol2 -c bcc -at gaff2 -nc {net_charge} -dr n  > ante_left.log
antechamber -i $RIGHT -fi pdb -o right.mol2 -fo mol2 -c bcc -at gaff2 -nc {net_charge} -dr n > ante_right.log
