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
dir=${10}

ehh=$(python -c "print '%.2f' % ($ehh_start)")
ehl=$(python -c "print '%.2f' % ($ehl_start)")
run=$run_start

while (( $(bc <<< "$ehh<=$ehh_end") ))
do
    ehl=$(python -c "print '%.2f' % ($ehl_start)")
    while (( $(bc <<< "$ehl<=$ehl_end") ))
    do
	for (( run=$run_start; $run<=$run_end; run+=$run_inc ))
	do
	    echo "Creating files for HH = $ehh HL = $ehl run = $run"
	    bash senescence.sh $ehh $ehl $run $dir
	done
	ehl=$(python -c "print '%.2f' % ($ehl + $ehl_inc)")
    done
    ehh=$(python -c "print '%.2f' % ($ehh + $ehh_inc)")
done
