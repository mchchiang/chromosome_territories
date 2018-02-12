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

L=40
chr=20

fit_file="${dir}/loop-prob_exponent.dat"
> $fit_file

while (( $(bc <<< "$ehh < $ehh_end") ))
do
    ehl=$(python -c "print '%.1f' % ($ehl_start)")
    while (( $(bc <<< "$ehl < $ehl_end") ))
    do
	echo "Fitting loop prob exponent for HH = ${ehh} HL = ${ehl}"

	name="sene_chr_${chr}_L_${L}_HH_${ehh}_HL_${ehl}"
	file="${dir}/loop-prob_${name}_avg.dat"

	fitparameters=$( echo "
f(x)=a*x+b
set fit errorvariables
fit [1.65:2.65] f(x) '$file' u (log10(\$1)):(log10(\$2)) via a,b
set print '-'
print a,a_err,b,b_err
" | gnuplot 2> /dev/null )
	
	echo "$ehh $ehl $fitparameters" >> $fit_file

	ehl=$(python -c "print '%.1f' % ($ehl + $ehl_inc)")
    done

    echo "" >> $fit_file

    ehh=$(python -c "print '%.1f' % ($ehh + $ehh_inc)")
done
