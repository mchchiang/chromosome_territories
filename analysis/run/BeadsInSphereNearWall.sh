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
rc=5
ld=200
t_start=150000
t_end=200000
t_inc=1000
z=3

ehh=$(python -c "print '%.1f' % ($ehh_start)")
ehl=$(python -c "print '%.1f' % ($ehl_start)")

source 'runconfig.cfg'
bis_exe="${bin}/contact/exe/BeadsInSphereNearWall"

max_jobs=8
cmd=()
jobid=0

while (( $(bc <<< "$ehh < $ehh_end") ))
do
    ehl=$(python -c "print '%.1f' % ($ehl_start)")
    while (( $(bc <<< "$ehl < $ehl_end") ))
    do
	name="wall_L_${L}_HH_${ehh}_HL_${ehl}"
	for (( run=$run_start; $run<=$run_end; run+=$run_inc ))
	do
	    echo "Calculating beads in sphere for HH = ${ehh} HL = ${ehl} run = ${run}"
	    pos_file="${dir}/pos_${name}_run_${run}.dat"
	    out_file="${dir}/bis_${name}_run_${run}.dat"
	    wall_file="${dir}/z_${name}_run_${run}.dat"
	    if [ -e $pos_file ]; then
		cmd[$jobid]="$bis_exe $N $L $L $L $rc $z $ld $t_start $t_end $t_inc $pos_file $out_file"
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
