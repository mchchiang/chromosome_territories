#!/bin/bash

# An input script for setting parameters in the LAMMPS config file

# Read in and set parameters
run=$1            # trial number
run_dir=$2        # run directory

gen_sphere_exe="../exe/Gen_46Chr_1BallEach"
chromo_file="../data/chromo_length.dat"
lam_file="../data/LAD.Pk.genome.full.dat"

box_size=200
lo=$(bc <<< "-$box_size/2.0")
hi=$(bc <<< "$box_size/2.0")

restart_freq=100

prep1_printfreq=1
prep1_seed=72394
prep1_time=10

prep2_printfreq=1
prep2_seed=53642
prep2_time=10

run_printfreq=1000
run_seed=184536
run_time=1000
max_samples=10

# Dimensions of the laminar shell
out_radius=71
in_radius=70.5

lam_atoms=5000
lam_seed=56213

delta_t=0.01       # time step size in Brownian time units

# Interaction energy parameters
# LJ potentials
sigma=13.0
low_e=0.5
low_cutoff=$(python -c "print '%.13f' % (2.0**(1.0/6.0)*$sigma)") # repulsive
mid_e=0.5
mid_cutoff=$(python -c "print '%.13f' % (1.8*$sigma)")
high_e=1.0
high_cutoff=$(python -c "print '%.13f' % (1.8*$sigma)")

# Normalisation (ensure minimum of potential is actually epsilon)
mid_norm=$(python -c "print '%.13f' % (1.0 + 4.0*(($sigma/$mid_cutoff)**12-($sigma/$mid_cutoff)**6))")
high_norm=$(python -c "print '%.13f' % (1.0 + 4.0*(($sigma/$high_cutoff)**12-($sigma/$high_cutoff)**6))")

echo "Mid norm = $mid_norm"
echo "High norm = $high_norm"

mid_e=$(python -c "print '%.13f' % ($mid_e/$mid_norm)")
high_e=$(python -c "print '%.13f' % ($high_e/$high_norm)")

echo "Mid E = $mid_e"
echo "High E = $high_e"

# Soft potentials
soft_a=10.0
soft_cutoff=50.0

# Set output file names
sim_name="46Chr_1BallEach_NLAM_${lam_atoms}_run_${run}"
init_file="sphere_${sim_name}.in"
restart_file="restart_${sim_name}"
prep1_outfile="prep1_${sim_name}.lammpstrj"
prep2_outfile="prep2_${sim_name}.lammpstrj"
run_outfile="run_${sim_name}.lammpstrj"
pos_file="pos_${sim_name}.dat"
sample_file="sample_${sim_name}"
map_file="run_${sim_name}.lammpsmap"

# Convert all time values to simulation time (i.e. rescale by delta t)
restart_freq=$(bc <<< "$restart_freq/$delta_t")
prep1_time=$(bc <<< "$prep1_time/$delta_t")
prep1_printfreq=$(bc <<< "$prep1_printfreq/$delta_t")
prep2_time=$(bc <<< "$prep2_time/$delta_t")
prep2_printfreq=$(bc <<< "$prep2_printfreq/$delta_t")
run_time=$(bc <<< "$run_time/$delta_t")
run_printfreq=$(bc <<< "$run_printfreq/$delta_t")

# Make execution directory
run_dir="${run_dir}/${sim_name}_run_${run}"
mkdir $run_dir

# Create the lammps command file based on template
lammps_file="${sim_name}.lam"
file="${run_dir}/${lammps_file}"
cp Sphere_with_Laminar.lam $file

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

sed -i -- "s/RUN_PRINTFREQ/${run_printfreq}/g" $file
sed -i -- "s/RUN_OUTFILE/${run_outfile}/g" $file
sed -i -- "s/RUN_SEED/${run_seed}/g" $file
sed -i -- "s/RUN_TIME/${run_time}/g" $file

sed -i -- "s/OUT_RADIUS/${out_radius}/g" $file
sed -i -- "s/IN_RADIUS/${in_radius}/g" $file
sed -i -- "s/LAM_ATOMS/${lam_atoms}/g" $file
sed -i -- "s/LAM_SEED/${lam_seed}/g" $file

sed -i -- "s/MAX_SAMPLES/${max_samples}/g" $file
sed -i -- "s/SAMPLE_FILE/${sample_file}/g" $file
sed -i -- "s/POS_FILE/${pos_file}/g" $file

sed -i -- "s/DELTA_T/${delta_t}/g" $file

sed -i -- "s/SIGMA/${sigma}/g" $file

sed -i -- "s/LOW_E/${low_e}/g" $file
sed -i -- "s/MID_E/${mid_e}/g" $file
sed -i -- "s/HIGH_E/${high_e}/g" $file
sed -i -- "s/LOW_CUTOFF/${low_cutoff}/g" $file
sed -i -- "s/MID_CUTOFF/${mid_cutoff}/g" $file
sed -i -- "s/HIGH_CUTOFF/${high_cutoff}/g" $file

sed -i -- "s/SOFT_A/${soft_a}/g" $file
sed -i -- "s/SOFT_CUTOFF/${soft_cutoff}/g" $file

# Generate spheres

${gen_sphere_exe} $chromo_file $lam_file "${run_dir}/${init_file}" "${run_dir}/${map_file}"
