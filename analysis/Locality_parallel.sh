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
ld=20
t_start=150000
t_end=200000
t_inc=1000

ehh=$(python -c "print '%.1f' % ($ehh_start)")
ehl=$(python -c "print '%.1f' % ($ehl_start)")

locality_exe="./Locality"
average_py="./GetAverage.py"
local_avg_file="${dir}/locality-local_avg_rc_${rc}_ld_${ld}.dat"
distal_avg_file="${dir}/locality-distal_avg_rc_${rc}_ld_${ld}.dat"
ratio_avg_file="${dir}/locality-ratio_avg_rc_${rc}_ld_${ld}.dat"
> $local_avg_file
> $distal_avg_file
> $ratio_avg_file

max_jobs=8
cmd=()
jobid=0

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
		cmd[$jobid]="$locality_exe $N $L $L $L $rc $ld $t_start $t_end $t_inc $pos_file $out_file"
		jobid=$(bc <<< "$jobid + 1")
	    fi
	done
	ehl=$(python -c "print '%.1f' % ($ehl + $ehl_inc)")
    done
    ehh=$(python -c "print '%.1f' % ($ehh + $ehh_inc)")
done

# Parallel runs

total_jobs=$jobid
jobid=0

while (( $(bc <<< "$jobid < $total_jobs") ))
do
    for (( i=0; i<$max_jobs && $jobid < $total_jobs; i++))
    do
	echo "${cmd[jobid]} &"
	${cmd[jobid]} &
	jobid=$(bc <<< "$jobid + 1")
    done
    wait
done

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
	python $average_py -1 0 -1 -1 ${file}_avg.dat ${file}*.dat
	data=$(cat ${file}_avg.dat)
	echo "${ehh} ${ehl} ${data}" >> $local_avg_file
	rm ${file}_avg.dat

	python $average_py -1 1 -1 -1 ${file}_avg.dat ${file}*.dat
	data=$(cat ${file}_avg.dat)
	echo "${ehh} ${ehl} ${data}" >> $distal_avg_file
	rm ${file}_avg.dat

	python $average_py -1 2 -1 -1 ${file}_avg.dat ${file}*.dat
	data=$(cat ${file}_avg.dat)
	echo "${ehh} ${ehl} ${data}" >> $ratio_avg_file
	rm ${file}_avg.dat

	ehl=$(python -c "print '%.1f' % ($ehl + $ehl_inc)")
    done
    echo >> $local_avg_file
    echo >> $distal_avg_file
    echo >> $ratio_avg_file
    ehh=$(python -c "print '%.1f' % ($ehh + $ehh_inc)")
done


