#!/bin/bash

ehh_start=$1
ehh_end=$2
ehh_inc=$3
ehl_start=$4
ehl_end=$5
ehl_inc=$6
dir=$7

ehh=$(python -c "print '%.1f' % ($ehh_start)")
ehl=$(python -c "print '%.1f' % ($ehl_start)")

source 'runconfig.cfg'
norm_exe="${bin}/contact/exe/VanillaNorm"
corr_exe="${bin}/contact/exe/CorrelationMatrix"
eigen_exe="${bin}/contact/exe/EigenAnalysis"
localdistal_exe="${bin}/contact/exe/LocalDistalScore"

# Selection arguments
N=1261
L=40
chr=20
rc=5
ld=40
matrix_mode="full"

while (( $(bc <<< "$ehh < $ehh_end") ))
do
    ehl=$(python -c "print '%.1f' % ($ehl_start)")
    while (( $(bc <<< "$ehl < $ehl_end") ))
    do
	echo "Processing contact for ehh = ${ehh} ehl = ${ehl}:"
	name="sene_chr_${chr}_L_${L}_HH_${ehh}_HL_${ehl}"
	contact_file="${dir}/contact_${name}_avg.dat"
	nocent_file="${dir}/contact_${name}_avg_nocent.dat"
	norm_file="${dir}/contact_${name}_avg_nocent_ice_norm.dat"
	eigen_file="${dir}/eigen_${name}_avg_nocent_ice_norm.dat"
	corr_file="${dir}/corr_${name}_avg_nocent_ice_norm.dat"
	localdistal_file="${dir}/local-distal_${name}_avg_nocent_ice_norm.dat"
	
#	echo "Remove centromeric region ..."
	# This assumes the bead index starts from 0 (not 1)
#	awk '{if (($1>=515&&$1<=595)||($2>=515&&$2<=595)) {$3=0.0; print} else {print} }' $contact_file > $nocent_file
	
#	echo "Normalising ..."
#	$norm_exe $N $matrix_mode $nocent_file $norm_file

	echo "Calculating correlation ..."
	$corr_exe $N $matrix_mode $norm_file $corr_file
	
	echo "Doing principal component analysis ..."
	$eigen_exe $N $matrix_mode $norm_file $eigen_file

	echo "Calculating local/distal scores ..."
	$localdistal_exe $N $ld $matrix_mode $norm_file $localdistal_file
	
	ehl=$(python -c "print '%.1f' % ($ehl + $ehl_inc)")
    done
    ehh=$(python -c "print '%.1f' % ($ehh + $ehh_inc)")
done

