#!/bin/bash

run_start=$1
run_end=$2
run_inc=$3
run_dir=$4

for (( i=$run_start;i<=$run_end;i+=$run_inc ))
do
    bash init.sh $i $run_dir
done
