#!/bin/bash

file_type=$1
ehh_start=$2
ehh_end=$3
ehh_inc=$4
ehl_start=$5
ehl_end=$6
ehl_inc=$7
run_start=$8
run_end=$9
run_inc=${10}
dir=${11}

ehh=$(python -c "print '%.1f' % ($ehh_start)")
ehl=$(python -c "print '%.1f' % ($ehl_start)")

select_lammps_exe="./GetBeadPositionLAMMPS"
select_pos_exe="./GetBeadPositionXYZ"

polymer_num=1
polymers="20"
type_num=1
types="2"

select_mode="cubic"
select_args="-20.0 20.0 -20.0 20.0 17.0 20.0"
L=40
chr=20

while (( $(bc <<< "$ehh < $ehh_end") ))
do
    ehl=$(python -c "print '%.1f' % ($ehl_start)")
    while (( $(bc <<< "$ehl < $ehl_end") ))
    do
	for (( run=$run_start; $run<=$run_end; run+=$run_inc ))
	do
	    echo "Selecting beads for HH = ${ehh} HL = ${ehl} run = ${run}"

	    name="sene_chr_${chr}_L_${L}_HH_${ehh}_HL_${ehl}_run_${run}"
	    lammps_file="${dir}/end_${name}.out"
	    map_file="${dir}/${name}.lammpsmap"
	    out_file="${dir}/selected_pos_${name}.dat"

	    echo $lammps_file
	    if [ -e $lammps_file ]; then
		temp="${dir}/selection.in~"
		echo > $temp
		echo $lammps_file >> $temp
		echo $map_file >> $temp
		echo $out_file >> $temp
		echo $polymer_num >> $temp
		echo $polymers >> $temp
		echo $type_num >> $temp
		echo $types >> $temp
		echo $select_mode >> $temp
		echo $select_args >> $temp
		
		if [ $file_type="lammps" ]; then
		    $select_lammps_exe $temp
		elif [ $file_type="pos" ]; then
		    $select_pos_exe $temp
		fi
		rm $temp
	    else
		echo "Cannot find position file!"
	    fi
	done
	ehl=$(python -c "print '%.1f' % ($ehl + $ehl_inc)")
    done
    ehh=$(python -c "print '%.1f' % ($ehh + $ehh_inc)")
done

