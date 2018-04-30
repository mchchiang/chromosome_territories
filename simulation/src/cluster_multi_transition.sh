#!/bin/bash

run_start=$1
run_end=$2
run_inc=$3
init_mode=$4
dir=$5

ehh=1.4
ehl=1.8
ehh2=1.4
ehl2=0.2
run=$run_start


for (( run=$run_start; $run<=$run_end; run+=$run_inc ))
do
    echo "Creating files for HH = $ehh HL = $ehl run = $run"
    bash cluster_transition.sh $ehh $ehl $ehh2 $ehl2 $run $init_mode $dir
done
