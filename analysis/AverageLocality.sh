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

N=6303
L=40
chr=20
rc=7
ld=2000
t_start=190000
t_end=200000
t_inc=1000

ehh=$(python -c "print '%.1f' % ($ehh_start)")
ehl=$(python -c "print '%.1f' % ($ehl_start)")

locality_exe="./Locality"
average_py="./GetAverage.py"
avg_file="${dir}/locality_avg_rc_${rc}_ld_${ld}.dat"
> $avg_file

max_jobs=8
cmd=()
jobid=0

# Average the data
echo "Averaging data"

ehh=$(python -c "print '%.1f' % ($ehh_start)")
ehl=$(python -c "print '%.1f' % ($ehl_start)")

while (( $(bc <<< "$ehh < $ehh_end") ))
do
    ehl=$(python -c "print '%.1f' % ($ehl_start)")
    while (( $(bc <<< "$ehl < $ehl_end") ))
    do
	echo "Averaging HH = ${ehh} HL = ${ehl}"
	file="${dir}/locality_cluster_chr_${chr}_L_${L}_HH_${ehh}_HL_${ehl}"
	python $average_py -1 2 -1 -1 ${file}_avg.dat ${file}*.dat
	data=$(cat ${file}_avg.dat)
	echo "${ehh} ${ehl} ${data}" >> $avg_file
	ehl=$(python -c "print '%.1f' % ($ehl + $ehl_inc)")
    done
    echo >> $avg_file
    ehh=$(python -c "print '%.1f' % ($ehh + $ehh_inc)")
done


