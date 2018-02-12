#!/bin/bash

ehh_start=$1
ehh_end=$2
ehh_inc=$3
ehl_start=$4
ehl_end=$5
ehl_inc=$6
dir=$7

ehh=$(python -c "print '%.1f' % ($ehh_start)")
ehl=$(python -c "print '%.1f' % ($ehl_start)")

average_py="./GetAverage.py"

L=40
chr=20

while (( $(bc <<< "$ehh < $ehh_end") ))
do
    ehl=$(python -c "print '%.1f' % ($ehl_start)")
    while (( $(bc <<< "$ehl < $ehl_end") ))
    do
	echo "Averaging loop prob for HH = ${ehh} HL = ${ehl}"

	name="sene_chr_${chr}_L_${L}_HH_${ehh}_HL_${ehl}"
	file="${dir}/loop-prob_${name}"
	
	python $average_py 0 1 -1 -1 ${file}_avg.dat ${file}_run*.dat
	
	ehl=$(python -c "print '%.1f' % ($ehl + $ehl_inc)")
    done
    ehh=$(python -c "print '%.1f' % ($ehh + $ehh_inc)")
done

