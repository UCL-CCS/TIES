script=$1

for D in */; 
do
    cd $D
    echo "Next Dir: $D"
    echo "Executing $script"
    python $script
    cd ..
done
