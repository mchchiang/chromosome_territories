#!/bin/bash

run_start=$1
run_end=$2
run_inc=$3
dir=$4

ehh=1.0
ehl_start=1.2
ehl_end=0.2
run=$run_start


for (( run=$run_start; $run<=$run_end; run+=$run_inc ))
do
    echo "Creating files for HH = $ehh HL = $ehl_start run = $run"
    bash senescence_hysteresis.sh $ehh $ehl_start $ehl_end $run $dir
done
