#!/bin/bash

ice_file=$1 # Iced matrix
norm_file=$2 # Normalised matrix

# Renormalise the contact matrix such that each column/row sums (close) to 1    
# Find the average column sum                                                   
echo "Renormalising contact matrix ..."
avg=$(awk '{if (NF==3) {if ($1 > max) {max = $1}; arr[$1]+=$3}} END {for (i = 0; i <= max; i++) {sum += arr[i]; if (arr[i]!=0) count++} print sum/count}' ${ice_file})

# Normalise the map by this average value                                       
awk -v norm="$avg" '{if (NF==3){printf("%d %d %.10g\n", $1, $2, $3/norm)} else {print}}' ${ice_file} > ${norm_file}

