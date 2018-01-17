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

ehh=$(python -c "print '%.1f' % ($ehh_start)")
ehl=$(python -c "print '%.1f' % ($ehl_start)")

time_avg_py="./GetTimeAverage.py"
avg_py="./GetAverage.py"

L=40
chr=20
tstart=150000

avg_file="${dir}/wall-frac_avg.dat"
> $avg_file

while (( $(bc <<< "$ehh < $ehh_end") ))
do
    ehl=$(python -c "print '%.1f' % ($ehl_start)")
    while (( $(bc <<< "$ehl < $ehl_end") ))
    do
	for (( run=$run_start; $run<=$run_end; run+=$run_inc ))
	do
	    echo "Averaging contact fraction for HH = ${ehh} HL = ${ehl} run ${run}"
	    
	    name="cluster_chr_${chr}_L_${L}_HH_${ehh}_HL_${ehl}"
	    file="${dir}/wall-frac_${name}"
	    python $time_avg_py $tstart 0 1 "${file}_run_${run}.dat" "${file}_run_${run}_avg.dat"
	done

	python $avg_py -1 0 -1 -1 ${file}_avg.dat ${file}_run*_avg.dat
	rm ${file}_run*_avg.dat
	
	data=$(cat ${file}_avg.dat) 
	echo "${ehh} ${ehl} ${data}" >> $avg_file
	rm ${file}_avg.dat

	ehl=$(python -c "print '%.1f' % ($ehl + $ehl_inc)")
    done

    echo >> $avg_file

    ehh=$(python -c "print '%.1f' % ($ehh + $ehh_inc)")
done

