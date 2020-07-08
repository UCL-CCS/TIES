script=$1

# execute the script in each dir

for D in */; 
do
    cd $D
    echo "Next Dir: $D"
    echo "Executing $script"
    python $script
    cd ..
done
