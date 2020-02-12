
# load the tools
source /opt/ambertools/amber.sh

# update to eval rather than ``
CURRENT_DIR=`dirname $0`
cd $CURRENT_DIR
echo 'Current Dir' $CURRENT_DIR

# generate the param file
tleap -s -f leap.in