#!/bin/bash

contact_file=$1
dir=$(realpath "${1}")
dir="${dir%/*}"
name=$(basename "$contact_file")
name="${name%.*}"
name="${name#*_}"

nocent_file="${dir}/contact_${name}_nocent.dat"
norm_file="${dir}/contact_${name}_nocent_norm.dat"
eigen_file="${dir}/eigen_${name}_nocent_norm.dat"
corr_file="${dir}/corr_${name}_nocent_norm.dat"
localdistal_file="${dir}/local-distal_${name}_nocent_norm.dat"

source 'runconfig.cfg'
norm_exe="${bin}/contact/exe/VanillaNorm"
corr_exe="${bin}/contact/exe/CorrelationMatrix"
eigen_exe="${bin}/contact/exe/EigenAnalysis"
localdistal_exe="${bin}/contact/exe/LocalDistalScore"

# Selection arguments
N=313
L=40
chr=20
rc=5
ld=10
matrix_mode="full"

echo "Remove centromeric region ..."
# 200 kb (131-146 for progeria)
# 50 kb (528-588)
# 10 kb (2640-2940)
awk '{if (($1>=131&&$1<=146)||($2>=131&&$2<=136)) {$3=0.0; print} else {print} }' $contact_file > $nocent_file

echo "Normalising ..."
$norm_exe $N $matrix_mode $nocent_file $norm_file

echo "Calculating correlation ..."
$corr_exe $N $matrix_mode $norm_file $corr_file

echo "Doing principal component analysis ..."
$eigen_exe $N $matrix_mode $norm_file $eigen_file

echo "Calculating local/distal scores ..."
$localdistal_exe $N $ld $matrix_mode $norm_file $localdistal_file
