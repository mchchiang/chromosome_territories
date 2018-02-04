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

locality_exe="./Locality"
average_py="./GetAverage.py"
avg_file="${dir}/locality_avg_ld_${rc}.dat"
> $avg_file

N=6303
L=40
chr=20
rc=5
ld=200
t_start=150000
t_end=200000
t_inc=1000

while (( $(bc <<< "$ehh < $ehh_end") ))
do
    ehl=$(python -c "print '%.1f' % ($ehl_start)")
    while (( $(bc <<< "$ehl < $ehl_end") ))
    do
	name="cluster_chr_${chr}_L_${L}_HH_${ehh}_HL_${ehl}"
	for (( run=$run_start; $run<=$run_end; run+=$run_inc ))
	do
	    echo "Calculating locality for HH = ${ehh} HL = ${ehl} run = ${run}"
	    pos_file="${dir}/pos_${name}_run_${run}.dat"
	    out_file="${dir}/locality_${name}_run_${run}.dat"

	    if [ -e $pos_file ]; then
		$locality_exe $N $L $L $L $rc $ld $t_start $t_end $t_inc $pos_file $out_file
	    fi
	done
	file="${dir}/locality_${name}"
	python $average_py -1 2 -1 -1 ${file}_avg.dat ${file}*.dat
	data=$(cat ${file}_avg.dat)
	echo "${ehh} ${ehl} ${data}" >> $avg_file
	ehl=$(python -c "print '%.1f' % ($ehl + $ehl_inc)")
    done
    echo >> $avg_file
    ehh=$(python -c "print '%.1f' % ($ehh + $ehh_inc)")
done

