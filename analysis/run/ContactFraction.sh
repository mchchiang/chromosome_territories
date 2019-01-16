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

ehh=$(python -c "print '%.1f' % ($ehh_start)")
ehl=$(python -c "print '%.1f' % ($ehl_start)")

source 'runconfig.cfg'
frac_exe="${bin}/contact/exe/ContactFraction"

# Selection arguments
wall_dist=3.0
bead_type="all"
t_start=0
t_end=2200000
t_inc=1000

N=6303
Nhet=3565 # 3079
Neu=2738 # 2923
Ncent=0 # 301
L=40
chr=20

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
	    pos_file="${in_dir}/pos_${name}.dat"
	    frac_file="${out_dir}/wall-frac_${name}.dat"
	    cmd[$jobid]="$frac_exe $N $Neu $Nhet $L $L $L $wall_dist $bead_type $t_start $t_end $t_inc $pos_file $frac_file"
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
