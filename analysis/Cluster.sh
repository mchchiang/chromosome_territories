#!/bin/bash

rc=$1
ehh_start=$2
ehh_end=$3
ehh_inc=$4
ehl_start=$5
ehl_end=$6
ehl_inc=$7
run_start=$8
run_end=$9
run_inc=${10}
dir=${11}

ehh=$(python -c "print '%.1f' % ($ehh_start)")
ehl=$(python -c "print '%.1f' % ($ehl_start)")

cluster_exe="./Cluster"
average_py="./GetAverage.py"
num_avg_file="${dir}/boxwide-cluster-num_avg_rc_${rc}.dat"
> $num_avg_file
size_avg_file="${dir}/boxwide-cluster-size_avg_rc_${rc}.dat"
> $size_avg_file

L=40
chr=20

while (( $(bc <<< "$ehh < $ehh_end") ))
do
    ehl=$(python -c "print '%.1f' % ($ehl_start)")
    while (( $(bc <<< "$ehl < $ehl_end") ))
    do
	name="cluster_chr_${chr}_L_${L}_HH_${ehh}_HL_${ehl}"
	for (( run=$run_start; $run<=$run_end; run+=$run_inc ))
	do
	    echo "Selecting beads for HH = ${ehh} HL = ${ehl} run = ${run}"
	    pos_file="${dir}/selected_pos_${name}_run_${run}.dat"
	    out_file="${dir}/cluster-stats_${name}_run_${run}.dat"

	    if [ -e $pos_file ]; then
		$cluster_exe $rc $pos_file $out_file
	    fi
	done
	file="${dir}/cluster-stats_${name}"
	python $average_py -1 0 -1 -1 ${file}_num-avg.dat ${file}*.dat
	data=$(cat ${file}_num-avg.dat)
	echo "${ehh} ${ehl} ${data}" >> $num_avg_file
	python $average_py -1 1 -1 -1 ${file}_size-avg.dat ${file}*.dat
	data=$(cat ${file}_size-avg.dat)
	echo "${ehh} ${ehl} ${data}" >> $size_avg_file
	ehl=$(python -c "print '%.1f' % ($ehl + $ehl_inc)")
    done
    echo >> $num_avg_file
    echo >> $size_avg_file
    ehh=$(python -c "print '%.1f' % ($ehh + $ehh_inc)")
done

