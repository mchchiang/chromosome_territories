#!/bin/bash

run_start=$1
run_end=$2
run_inc=$3
init_mode=$4
dir=$5

ehh=1.4
ehl_start=1.2
ehl_end=0.6
run=$run_start


for (( run=$run_start; $run<=$run_end; run+=$run_inc ))
do
    echo "Creating files for HH = $ehh HL = $ehl run = $run"
    bash cluster_hysteresis.sh $ehh $ehl_start $ehl_end $run $init_mode $dir
done
