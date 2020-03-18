
# load the tools
source ~/software/amber18install/amber.sh

# update to eval rather than ``
CURRENT_DIR=`dirname $0`
cd $CURRENT_DIR
echo 'Current Dir' $CURRENT_DIR

# generate the param file
tleap -s -f leap_complex.in