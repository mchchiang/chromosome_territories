#!/bin/bash

ehh_start=$1
ehh_end=$2
ehh_inc=$3
ehl_start=$4
ehl_end=$5
ehl_inc=$6
run_start=$7
run_end=$8
run_inc=$9

ehh=$(python -c "print '%.1f' % ($ehh_start)")
ehl=$(python -c "print '%.1f' % ($ehl_start)")

chr=20
L=40
N=6303

combine_exe="../analysis/bin/exe/CombinePosXYZ"

while (( $(bc <<< "$ehh < $ehh_end") ))
do
    ehl=$(python -c "print '%.1f' % ($ehl_start)")
    while (( $(bc <<< "$ehl < $ehl_end") ))
    do
	for (( run=$run_start; $run<=$run_end; run+=$run_inc ))
	do
	    dir="sene_chr_${chr}_L_${L}_HH_${ehh}_HL_${ehl}_run_${run}"
	    run_file="${dir}/run_${dir}.lammpstrj"
	    pos_file="${dir}/pos_${dir}.dat"
	    echo "Generating position file for ehh ${ehh} ehl ${ehl} run ${run} ... "
	    $combine_exe $N $N $run_file $pos_file
	done
	ehl=$(python -c "print '%.1f' % ($ehl + $ehl_inc)")
    done
    ehh=$(python -c "print '%.1f' % ($ehh + $ehh_inc)")
done
