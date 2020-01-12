#!/bin/bash

in_dir=$1
out_dir=$2

ehh=1.0
ehl=1.6

ehh=$(python -c "print '%.2f' % ($ehh)")
ehl=$(python -c "print '%.2f' % ($ehl)")

source 'runconfig.cfg'
locality_exe="${bin}/contact/exe/LocalDistalScore" # Program that computes OCI
average_py="../src/AverageOneFile.py"

N=1261 # For senescence experiments
#N=316 # For progeria experiments
Nm1=$(python -c "print $N-1")
L=35
chr=20
rc=3
ld_start=1
ld_end=100 # For senescence experiments
#ld_end=25 # For progeria experiments (lower resolutions 200kb)
ld_inc=1
ld=$ld_start
mode="full"

# For simulations
name="sene_chr_${chr}_L_${L}_HH_${ehh}_HL_${ehl}"
contact_file="${in_dir}/contact_${name}_avg_nocent.dat"

# For experiments
#name="hic_growing_50000_iced_chr_20_norm"
#name="hic_senescence_50000_iced_chr_20_norm"
#name="hic_age-control_iced_chr_20_norm"
#name="hic_hgps-p19_iced_chr_20_norm"
#contact_file="${in_dir}/contact_${name}.dat"

#name="hic_growing"
#name="hic_senescence"
#name="hic_age-control"
#name="hic_hgps-p19"

vary_file="${out_dir}/vary-locality_${name}.dat"

> $vary_file
for (( ld=$ld_start; $ld<=$ld_end; ld+=$ld_inc ))
do
    echo "Doing HH = ${ehh} HL = ${ehl} ld = ${ld}"
    oci_file="${out_dir}/locality_${name}_ld_${ld}.dat"
    oci_avg_file="${out_dir}/locality_${name}_ld_${ld}_avg.dat"
    $locality_exe $N $ld $mode $contact_file $oci_file
    
    # Average the OCI values over all beads
    python $average_py 3 -1 0 $Nm1 $oci_file $oci_avg_file
    data=$(cat $oci_avg_file)
    echo "${ld} ${data}" >> $vary_file
    rm $oci_file $oci_avg_file
done

