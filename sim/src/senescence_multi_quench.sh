#!/bin/bash

ehh1=1.0
ehl1=1.2
ehh2=1.0
ehl2=0.2
run_start=$1
run_end=$2
run_inc=$3
dir=$4

ehh1=$(python -c "print '%.2f' % ($ehh1)")
ehl1=$(python -c "print '%.2f' % ($ehl1)")
ehh2=$(python -c "print '%.2f' % ($ehh2)")
ehl2=$(python -c "print '%.2f' % ($ehl2)")
run=$run_start

for (( run=$run_start; $run<=$run_end; run+=$run_inc ))
do
    echo "Creating files for HH = $ehh HL = $ehl run = $run"
    bash senescence_quench.sh $ehh1 $ehl1 $ehh2 $ehl2 $run $dir
done
