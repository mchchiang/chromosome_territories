#!/bin/bash

ehh_start=$1
ehh_end=$2
ehh_inc=$3
ehl_start=$4
ehl_end=$5
ehl_inc=$6
in_dir=$7
out_dir=$8

zero=0.0
near_zero=0.01
diff_zero=0.00001

N=6303
Nhet=3513 # 3565 # 3079
Neu=2790 # 2738 # 2923
Ncent=0 # 301
bead_type="all"

L=35
chr=20
t_start=150000
t_end=200000
t_inc=1000
min=0.0
max=35.0
binsize=1.0

source 'runconfig.cfg'
prob_exe="${bin}/contact/exe/ProbFromWall"

ehh=$(python -c "print '%.2f' % ($ehh_start)")
ehl=$(python -c "print '%.2f' % ($ehl_start)")

max_jobs=8
cmd=()
jobid=0

while (( $(bc <<< "$ehh < $ehh_end") ))
do
    ehh=$(python -c "print '%.2f' % ($near_zero if abs($ehh)<$diff_zero else $ehh)")
    ehl=$(python -c "print '%.2f' % ($ehl_start)")
    while (( $(bc <<< "$ehl < $ehl_end") ))
    do
	ehl=$(python -c "print '%.2f' % ($near_zero if abs($ehl)<$diff_zero else $ehl)")
	name="sene_chr_${chr}_L_${L}_HH_${ehh}_HL_${ehl}"
	pos_file="${in_dir}/pos_${name}"*.dat
	out_file="${out_dir}/zpos-prob_${name}.dat"
	cmd[$jobid]="$prob_exe $N $L $L $L $min $max $binsize $bead_type $t_start $t_end $t_inc $out_file ${pos_file}"
	jobid=$(bc <<< "$jobid + 1")
	
	ehl=$(python -c "print '%.2f' % ($zero if abs($ehl-$near_zero)<$diff_zero else $ehl)")	
	ehl=$(python -c "print '%.2f' % ($ehl + $ehl_inc)")
    done
    ehh=$(python -c "print '%.2f' % ($zero if abs($ehh-$near_zero)<$diff_zero else $ehh)")
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
