#!/bin/bash

ehh_start=$1
ehh_end=$2
ehh_inc=$3
ehl_start=$4
ehl_end=$5
ehl_inc=$6
dir=$7

zero=0.0
near_zero=0.01
diff_zero=0.00001

N=6303
L=35
chr=20
rc=5
ld=200
t_start=150000
t_end=200000
t_inc=1000

ehh=$(python -c "print '%.1f' % ($ehh_start)")
ehl=$(python -c "print '%.1f' % ($ehl_start)")

source 'runconfig.cfg'
average_py="../src/AverageMultiFiles.py"

total_file="${dir}/bis_sene_rc_${rc}_ld_${ld}.dat"
> $total_file


# Average the data
echo "Averaging data"

ehh=$(python -c "print '%.2f' % ($ehh_start)")
ehl=$(python -c "print '%.2f' % ($ehl_start)")

while (( $(bc <<< "$ehh < $ehh_end") ))
do
    ehh=$(python -c "print '%.2f' % ($near_zero if abs($ehh)<$diff_zero else $ehh)")
    ehl=$(python -c "print '%.2f' % ($ehl_start)")
    while (( $(bc <<< "$ehl < $ehl_end") ))
    do
	ehl=$(python -c "print '%.2f' % ($near_zero if abs($ehl)<$diff_zero else $ehl)")
	echo "Averaging HH = ${ehh} HL = ${ehl}"
	file="${dir}/bis_sene_chr_${chr}_L_${L}_HH_${ehh}_HL_${ehl}"
	local_bis_avg_file=${dir}/local-bis_avg.dat
	distal_bis_avg_file=${dir}/distal-bis_avg.dat
	local_density_avg_file=${dir}/local-density_avg.dat
	distal_density_avg_file=${dir}/distal-density_avg.dat
	bis_avg_file=${dir}/bis_avg.dat
	density_avg_file=${dir}/density_avg.dat

	# Average all columns
	# Local beads in sphere
	python $average_py -1 0 -1 -1 $local_bis_avg_file ${file}*.dat
	python $average_py -1 1 -1 -1 $distal_bis_avg_file ${file}*.dat
	python $average_py -1 2 -1 -1 $bis_avg_file ${file}*.dat
	python $average_py -1 3 -1 -1 $local_density_avg_file ${file}*.dat
	python $average_py -1 4 -1 -1 $distal_density_avg_file ${file}*.dat
	python $average_py -1 5 -1 -1 $density_avg_file ${file}*.dat

	data=$(paste -d" " $local_bis_avg_file $distal_bis_avg_file $bis_avg_file $local_density_avg_file $distal_density_avg_file $density_avg_file)
	echo $ehh $ehl $data >> $total_file
	rm $local_bis_avg_file
	rm $distal_bis_avg_file
	rm $bis_avg_file
	rm $local_density_avg_file
	rm $distal_density_avg_file
	rm $density_avg_file
	
	ehl=$(python -c "print '%.2f' % ($zero if abs($ehl-$near_zero)<$diff_zero else $ehl)")	
	ehl=$(python -c "print '%.2f' % ($ehl + $ehl_inc)")
    done
    echo >> $total_file
    ehh=$(python -c "print '%.2f' % ($zero if abs($ehh-$near_zero)<$diff_zero else $ehh)")
    ehh=$(python -c "print '%.2f' % ($ehh + $ehh_inc)")
done

