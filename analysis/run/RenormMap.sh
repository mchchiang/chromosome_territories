#!/bin/bash

ice_file=$1 # Iced matrix
norm_file=$2 # Normalised matrix

# Renormalise the contact matrix such that each column/row sums (close) to 1    
# Find the average column sum                                                   
echo "Renormalising contact matrix ..."
#max=$(awk '{if (NF==3) {if ($1 > max) {max = $1}; arr[$1]+=$3}} END {for (i = 0; i <= max; i++) {if (arr[i] > maxCount) {maxCount = arr[i]}} print maxCount}' ${ice_file})

#avg=$(awk '{if (NF==3) {if ($1 > max) {max = $1}; arr[$1]+=$3}} END {for (i = 0; i <= max; i++) {sum += arr[i]; if (arr[i]!=0) count++} print sum/count}' ${ice_file})

diagavg=$(awk '{if(NF==3){if($1>max){max=$1};if(arr[$1]<$3){arr[$1]=$3}}}END{for(i=0;i<=max;i++){if(arr[i]>0){sum+=arr[i];count++}};print sum/count}' ${ice_file})

# Normalise the map by this average value
awk -v norm="$diagavg" '{if (NF==3){printf("%d %d %.10g\n", $1, $2, $3/norm)} else {print}}' ${ice_file} > ${norm_file}
