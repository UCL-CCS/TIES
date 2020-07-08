#!/bin/bash
## !!!! NO ARGUMENTS REQUIRED !!!!
## $1 and $2 are the first and second command line arguments which are starting and ending replica index.
## $3: number of REST2 replicas

#if [ $# -ne 3 ];then
#    echo "Incorrect no of arguments"
#    exit
#else
#    rep_start=$1
#    rep_end=$2
#    num_rep_rest=$3
#fi

num_rep_rest=13

dwork="$PWD"
mkdir -p $PWD/analysis/t-ser-40ns
dir_out="$PWD/analysis/t-ser-40ns"

#for ((i=rep_start; i<=rep_end; i++)); do
for rep in `ls -d rep{1..10}`; do
    for ((v=0; v<num_rep_rest; v++)); do
        cd $dwork/$rep/simulation/$v/energy_all-40ns
        l1=`grep "PARTITION 1 VDW" sim1-alch.dat | awk '{print $NF}' | head -1`
        l2=`grep "PARTITION 2 VDW" sim1-alch.dat | awk '{print $NF}' | head -1`
        le1=`grep "PARTITION 1 ELEC" sim1-alch.dat | awk '{print $NF}' | head -1`
        le2=`grep "PARTITION 2 ELEC" sim1-alch.dat | awk '{print $NF}' | head -1`
        for f in `ls sim1-alch.dat`; do
            if [[ $l1 == 0 && $l2 == 0 ]]; then
               awk '($1=="TI") {print 0.0}' $f >> $dir_out/$rep-vdw-$v.dat
            elif [[ $l1 == 0 ]]; then
               awk '($1=="TI") {print -$9}' $f >> $dir_out/$rep-vdw-$v.dat
            elif [[ $l2 == 0 ]]; then
               awk '($1=="TI") {print $5}' $f >> $dir_out/$rep-vdw-$v.dat
            else
               awk '($1=="TI") {print $5-$9}' $f >> $dir_out/$rep-vdw-$v.dat
            fi

            if [[ $le1 == 0 && $le2 == 0 ]]; then
               awk '($1=="TI") {print 0.0}' $f >> $dir_out/$rep-ele-$v.dat
            elif [[ $le1 == 0 ]]; then
               awk '($1=="TI") {print -$7}' $f >> $dir_out/$rep-ele-$v.dat
            elif [[ $le2 == 0 ]]; then
               awk '($1=="TI") {print $3}' $f >> $dir_out/$rep-ele-$v.dat
            else
               awk '($1=="TI") {print $3-$7}' $f >> $dir_out/$rep-ele-$v.dat
            fi
        done
    done
done
