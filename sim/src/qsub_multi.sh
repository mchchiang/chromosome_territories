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

chr=20
L=32

qsub_sh="./qsub_script.2.sh"

while (( $(bc <<< "$ehh < $ehh_end") ))
do
    ehl=$(python -c "print '%.1f' % ($ehl_start)")
    while (( $(bc <<< "$ehl < $ehl_end") ))
    do
	for (( run=$run_start; $run<=$run_end; run+=$run_inc ))
	do
	    sim_name="sene_chr_${chr}_L_${L}_HH_${ehh}_HL_${ehl}_run_${run}"
	    sim_in_dir="${in_dir}/${sim_name}"
	    sim_out_dir="${out_dir}/${sim_name}"
	    if [ ! -d $sim_out_dir ]; then
		mkdir -p $sim_out_dir
	    fi
	    qsub -o $sim_out_dir/stdout.log -e $sim_out_dir/stderr.log $qsub_sh $sim_name $sim_in_dir $sim_out_dir
	done
	ehl=$(python -c "print '%.1f' % ($ehl + $ehl_inc)")
    done
    ehh=$(python -c "print '%.1f' % ($ehh + $ehh_inc)")
done
