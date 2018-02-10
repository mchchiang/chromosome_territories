#!/bin/bash

# An input script for setting parameters in the LAMMPS config file

# Read in and set parameters
run=$1            # trial number
src_dir=$2        # source directory - where the sphere files are
run_dir=$3        # run directory

gen_chromo_exe="../exe/Gen_46Chr_FromBalls"
chromo_file="../data/chromo_length.dat"
prom_file="../data/H3K4me3.Pk.genome.full.dat"
lam_file="../data/LAD.Pk.genome.full.dat"
pcg_file="../data/PCG.genome.full.dat"

box_size=200
lo=$(bc <<< "-$box_size/2.0")
hi=$(bc <<< "$box_size/2.0")

restart_freq=100

prep1_printfreq=1
prep1_seed=23763
prep1_time=10

prep2_printfreq=1
prep2_seed=9832
prep2_time=100

prep3_seed=11873
prep3_time=10

# Dimensions of the laminar shell
out_radius=71
in_radius=70.5

lam_atoms=5000
lam_seed=7364

delta_t=0.01       # time step size in Brownian time units

# Source file names
sphere_file="${src_dir}/sample_$(basename ${src_dir})_${run}.dat"
sphere_map_file="${src_dir}/$(basename ${src_dir}).lammpsmap"

# Set output file names
sim_name="46Chr_FromBall_NLAM_${lam_atoms}"
init_file="chromo_${sim_name}.in"
restart_file="restart_${sim_name}"
prep1_outfile="prep1_${sim_name}.lammpstrj"
prep2_outfile="prep2_${sim_name}.lammpstrj"
pos_file="pos_${sim_name}.dat"
precompress_file="precompress_${sim_name}.dat"
compress_file="compress_${sim_name}.dat"
warmedup_file="warmedup_${sim_name}.dat"
map_file="run_${sim_name}.lammpsmap"

# Convert all time values to simulation time (i.e. rescale by delta t)
restart_freq=$(bc <<< "$restart_freq/$delta_t")
prep1_time=$(bc <<< "$prep1_time/$delta_t")
prep1_printfreq=$(bc <<< "$prep1_printfreq/$delta_t")
prep2_time=$(bc <<< "$prep2_time/$delta_t")
prep2_printfreq=$(bc <<< "$prep2_printfreq/$delta_t")

# Make execution directory
run_dir="${run_dir}/${sim_name}_run_${run}"
mkdir $run_dir

# Create the lammps command file based on template
lammps_file="${sim_name}.lam"
file="${run_dir}/${lammps_file}"
cp Chromo_with_Laminar_WarmUp.lam $file

# Replace macros in template with input values
sed -i -- "s/INIT_FILE/${init_file}/g" $file
sed -i -- "s/RESTART_FILE/${restart_file}/g" $file

sed -i -- "s/XLO/${lo}/g" $file
sed -i -- "s/XHI/${hi}/g" $file
sed -i -- "s/YLO/${lo}/g" $file
sed -i -- "s/YHI/${hi}/g" $file
sed -i -- "s/ZLO/${lo}/g" $file
sed -i -- "s/ZHI/${hi}/g" $file

sed -i -- "s/RESTART_FREQ/${restart_freq}/g" $file

sed -i -- "s/PREP1_PRINTFREQ/${prep1_printfreq}/g" $file
sed -i -- "s/PREP1_OUTFILE/${prep1_outfile}/g" $file
sed -i -- "s/PREP1_SEED/${prep1_seed}/g" $file
sed -i -- "s/PREP1_TIME/${prep1_time}/g" $file

sed -i -- "s/PREP2_PRINTFREQ/${prep2_printfreq}/g" $file
sed -i -- "s/PREP2_OUTFILE/${prep2_outfile}/g" $file
sed -i -- "s/PREP2_SEED/${prep2_seed}/g" $file
sed -i -- "s/PREP2_TIME/${prep2_time}/g" $file

sed -i -- "s/PREP3_SEED/${prep3_seed}/g" $file
sed -i -- "s/PREP3_TIME/${prep3_time}/g" $file

sed -i -- "s/OUT_RADIUS/${out_radius}/g" $file
sed -i -- "s/IN_RADIUS/${in_radius}/g" $file
sed -i -- "s/LAM_ATOMS/${lam_atoms}/g" $file
sed -i -- "s/LAM_SEED/${lam_seed}/g" $file

sed -i -- "s/PRECOMPRESS_FILE/${precompress_file}/g" $file
sed -i -- "s/COMPRESS_FILE/${compress_file}/g" $file
sed -i -- "s/WARMEDUP_FILE/${warmedup_file}/g" $file

sed -i -- "s/DELTA_T/${delta_t}/g" $file

# Generate polymers

${gen_chromo_exe} $chromo_file $prom_file $lam_file $pcg_file $sphere_file $sphere_map_file "${run_dir}/${init_file}" "${run_dir}/${map_file}"
