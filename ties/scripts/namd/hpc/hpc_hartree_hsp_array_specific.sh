#!/bin/bash
# define an array job
#BSUB -J mcl_spec[1-1]
#BSUB -o std.%J.%I.o
#BSUB -e std.%J.%I.e
#BSUB -R "span[ptile=32]"
#BSUB -n 32
##BSUB -q compbiomed
#BSUB -q scafellpikeSKL # can be used a litte, 10x4, try 65x1 (does not work). Seems to allow only for 10 at a time.
##BSUB -q scafellpikeKNL # can be used a little, 12x4, try 65x1, but the performance is terrible - it does not really work
##BSUB -q universeScafellPike - cannot be used
##BSUB -q scafellpikeI - cannot be used
#BSUB -W 30:00
#BSUB -x

SYSTEM=$HCBASE/resp_ff/mcl1_l12_l35/lig/lambda_0.90/rep1
NP=32

#Load modules
source /etc/profile.d/modules.sh
module load namd-gcc/2.12

cmd_mpinamd="mpiexec -n $NP namd2 +pemap 1-31 +commap 0"

# check the ligand first
# lig min
if test -f min.log && grep -q "WRITING VELOCITIES TO OUTPUT FILE" min.log; then
    echo "LIG MIN ALREADY DONE"
else
    echo "LIG MIN"
	${cmd_mpinamd} min.namd > min.log
fi
# eq1
if test -f eq_step1.log && grep -q "WRITING EXTENDED SYSTEM TO OUTPUT FILE AT STEP 100000" eq_step1.log; then
    echo "LIG EQ1 ALREADY DONE"
else
    echo "LIG EQ1"
	${cmd_mpinamd} eq_step1.namd > eq_step1.log
fi
# eq2
if test -f eq_step2.log && grep -q "WRITING EXTENDED SYSTEM TO OUTPUT FILE AT STEP 100000" eq_step2.log; then
    echo "LIG EQ2 ALREADY DONE"
else
    echo "LIG EQ2"
	${cmd_mpinamd} eq_step2.namd > eq_step2.log
fi
# eq3
if test -f eq_step3.log && grep -q "WRITING EXTENDED SYSTEM TO OUTPUT FILE AT STEP 100000" eq_step3.log; then
    echo "LIG EQ3 ALREADY DONE"
else
    echo "LIG EQ3"
	${cmd_mpinamd} eq_step3.namd > eq_step3.log
fi
# eq4
if test -f eq_step4.log && grep -q "WRITING EXTENDED SYSTEM TO OUTPUT FILE AT STEP 100000" eq_step4.log; then
    echo "LIG EQ4 ALREADY DONE"
else
    echo "LIG EQ4"
	${cmd_mpinamd} eq_step4.namd > eq_step4.log
fi

# prod
if test -f prod.log && grep -q "WRITING VELOCITIES TO OUTPUT FILE AT STEP 3000000" prod.log; then
    echo "LIG PROD ALREADY DONE"
else
    echo "LIG PROD"
	${cmd_mpinamd} prod.namd > prod.log
fi