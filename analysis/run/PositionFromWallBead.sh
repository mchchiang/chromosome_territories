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
in_dir=${10}
out_dir=${11}

ehh=$(python -c "print '%.2f' % ($ehh_start)")
ehl=$(python -c "print '%.2f' % ($ehl_start)")

source 'runconfig.cfg'
pos_exe="${bin}/contact/exe/PositionFromWallBead"

# Selection arguments
#t_start=150000
#t_end=200000
t_start=2000
t_end=2000
t_inc=10

N=6303
L=35
chr=20

max_jobs=10
cmd=()
jobid=0

while (( $(bc <<< "$ehh < $ehh_end") ))
do
    ehh=$(python -c "print '%.2f' % (0.01 if abs($ehh)<0.0001 else $ehh)")
    ehl=$(python -c "print '%.2f' % ($ehl_start)")
    while (( $(bc <<< "$ehl < $ehl_end") ))
    do
	ehl=$(python -c "print '%.2f' % (0.01 if abs($ehl)<0.0001 else $ehl)")
	for (( run=$run_start; $run<=$run_end; run+=$run_inc ))
	do
	    name="sene_chr_${chr}_L_${L}_HH_${ehh}_HL_${ehl}_run_${run}"
	    pos_file="${in_dir}/pos_${name}.dat"
	    #frac_file="${out_dir}/zpos-bead_${name}.dat"
	    frac_file="${out_dir}/zpos-bead_${name}_t_${t_start}.dat"
	    cmd[$jobid]="$pos_exe $N $L $L $L $t_start $t_end $t_inc $frac_file $pos_file"
	    jobid=$(bc <<< "$jobid + 1")

	done
	ehl=$(python -c "print '%.2f' % (0.0 if abs($ehl-0.01)<0.0001 else $ehl)")	
	ehl=$(python -c "print '%.2f' % ($ehl + $ehl_inc)")
    done
    ehh=$(python -c "print '%.2f' % (0.0 if abs($ehh-0.01)<0.0001 else $ehh)")
    ehh=$(python -c "print '%.2f' % ($ehh + $ehh_inc)")
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
