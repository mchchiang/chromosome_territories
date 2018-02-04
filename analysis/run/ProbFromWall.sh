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

prob_py="./ProbFromWall.py"

L=40
chr=20
min_val=0.0
max_val=40.0
binsize=1.0
wall_pos=20.0

while (( $(bc <<< "$ehh < $ehh_end") ))
do
    ehl=$(python -c "print '%.1f' % ($ehl_start)")
    while (( $(bc <<< "$ehl < $ehl_end") ))
    do
	echo "Calculating prob. for HH = ${ehh} HL = ${ehl}"
	name="cluster_chr_${chr}_L_${L}_HH_${ehh}_HL_${ehl}"
	pos_file="${dir}/selected_pos_${name}"
	out_file="${dir}/wall-prob_${name}.dat"
	python $prob_py $wall_pos $min_val $max_val $binsize $out_file ${pos_file}*.dat
	ehl=$(python -c "print '%.1f' % ($ehl + $ehl_inc)")
    done
    ehh=$(python -c "print '%.1f' % ($ehh + $ehh_inc)")
done

