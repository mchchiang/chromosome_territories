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
L=32
N=6303

combine_exe="../analysis/bin/exe/CombinePosXYZ"

max_jobs=10
cmd=()
jobid=0

while (( $(bc <<< "$ehh < $ehh_end") ))
do
    ehl=$(python -c "print '%.1f' % ($ehl_start)")
    while (( $(bc <<< "$ehl < $ehl_end") ))
    do
	for (( run=$run_start; $run<=$run_end; run+=$run_inc ))
	do
	    name="sene_chr_${chr}_L_${L}_HH_${ehh}_HL_${ehl}_run_${run}"
	    prep1_traj="${name}/prep1_${name}.lammpstrj"
	    prep2_traj="${name}/prep2_${name}.lammpstrj"
	    prep3_traj="${name}/prep3_${name}.lammpstrj"
	    run_traj="${name}/run_${name}.lammpstrj"
	    prep1_file="${name}/prep1_${name}.dat"
	    prep2_file="${name}/prep2_${name}.dat"
	    prep3_file="${name}/prep3_${name}.dat"
	    pos_file="${name}/pos_${name}.dat"
	    cmd[$jobid]="$combine_exe $N $N $prep1_traj $prep1_file"
	    jobid=$(bc <<< "$jobid + 1")
	    cmd[$jobid]="$combine_exe $N $N $prep2_traj $prep2_file"
	    jobid=$(bc <<< "$jobid + 1")
	    cmd[$jobid]="$combine_exe $N $N $prep3_traj $prep3_file"
	    jobid=$(bc <<< "$jobid + 1")
	    cmd[$jobid]="$combine_exe $N $N $run_traj $pos_file"
	    jobid=$(bc <<< "$jobid + 1")
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
