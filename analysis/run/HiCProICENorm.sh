#!/bin/bash

ehh=$1
ehl=$2
config_file=$3
in_dir=$4
out_dir=$5

ehh=$(python -c "print '%.1f' % ($ehh)")
ehl=$(python -c "print '%.1f' % ($ehl)")

source 'runconfig.cfg'
convert_exe="${bin}/contact/exe/ConvertMatrix"

# Selection arguments
N=1261
L=40
chr=20
reads=10000

name="contact_sene_chr_${chr}_L_${L}_HH_${ehh}_HL_${ehl}"
contact_file="${in_dir}/${name}_avg_nocent.dat"

# Create the directory for storing the result
echo "Creating the directory for storing the result ..."
mkdir ${out_dir}/${name}
out_dir=${out_dir}/${name}

raw_dir=${out_dir}/raw
ice_dir=${out_dir}/ice

mkdir ${raw_dir}
mkdir ${raw_dir}/${name}
mkdir ${ice_dir}

cp $contact_file ${raw_dir}/${name}/${name}.dat

# Amplify the probability by a constant to mimic reads from HiC
echo "Amplifying the contact probability matrix ..."
awk -v n="$reads" '{if(NF==3){print $1, $2, $3*n}}' ${raw_dir}/${name}/${name}.dat > ${raw_dir}/${name}/${name}_50000.matrix
rm ${raw_dir}/${name}/${name}.dat 

echo "Doing ICE normalistion with HiC-Pro ..."
HiC-Pro -i ${raw_dir} -o ${ice_dir} -c $config_file -s ice_norm
	
# Move results to the top directory within ice
mv ${ice_dir}/hic_results/matrix/${name}/iced/50000/* ${out_dir}

# Convert the normalised map to the format read by gnuplot
echo "Converting result to full matrix ..."
$convert_exe $N upper full ${out_dir}/${name}_50000_iced.matrix ${out_dir}/${name}_avg_nocent_ice.dat

# Renormalise the contact matrix such that each column/row sums (close) to 1
# Find the average column sum
echo "Renormalising contact matrix ..."
avg=$(awk '{if (NF==3) {if ($1 > max) {max = $1}; arr[$1]+=$3}} END {for (i = 0; i <= max; i++) {sum += arr[i]; if (arr[i]!=0) count++} print sum/count}' ${out_dir}/${name}_avg_nocent_ice.dat)

# Normalise the map by this average value
awk -v norm="$avg" '{if (NF==3){printf("%d %d %.10g\n", $1, $2, $3/norm)} else {print}}' ${out_dir}/${name}_avg_nocent_ice.dat > ${out_dir}/${name}_avg_nocent_ice_norm.dat

# Clean up
echo "Clean up ..."

rm -r ${ice_dir}
rm -r ${raw_dir}

mv ${out_dir}/${name}_50000_iced.matrix.biases ${out_dir}/${name}_avg_nocent_ice_bias.dat
rm ${out_dir}/${name}_50000_iced.matrix
