#!/bin/bash

# An input script for setting parameters in the LAMMPS config file

# Read in and set parameters
run=$1            # trial number
run_dir=$2        # run directory

gen_chromo_exe="../exe/Gen_Chr_Het-LAD"
chromo_file="../data/chromo_length.dat"
lam_file="../data/LAD.Pk.genome.full.dat"
het_file="../data/GM12878.H3K9me3.Pk.full.dat"
chr_num=20

box_size=50
buffer=5
lo=$(bc <<< "-$box_size/2.0")
hi=$(bc <<< "$box_size/2.0")

max_seed=100000

restart_freq=100

prep1_printfreq=10
prep1_seed=$(python GetRandom.py $max_seed)
prep1_time=1000

prep2_printfreq=10
prep2_seed=$(python GetRandom.py $max_seed)
prep2_time=1000

run_printfreq=1000
run1_seed=$(python GetRandom.py $max_seed)
run1_time=10000
run2_seed=$(python GetRandom.py $max_seed)
run2_time=250000
run3_seed=$(python GetRandom.py $max_seed)
run3_time=250000

lam_atoms=2000
lam_seed=$(python GetRandom.py $max_seed)

delta_t=0.01       # time step size in Brownian time units

# Interaction energy parameters
# LJ potentials
sigma=1.0
cutoff=$(python -c "print '%.13f' % (1.8*$sigma)")
e_hethet=1.0
e_hetlam=1.0
e_hetlad=1.0
e_ladlad=1.0
e_ladlam=2.0

# Normalisation (ensure minimum of potential is actually epsilon)
norm=$(python -c "print '%.13f' % (1.0 + 4.0*(($sigma/$cutoff)**12-($sigma/$cutoff)**6))")

e_hethet_norm=$(python -c "print '%.13f' % ($e_hethet/$norm)")
e_hetlam_norm=$(python -c "print '%.13f' % ($e_hetlam/$norm)")
e_hetlad_norm=$(python -c "print '%.13f' % ($e_hetlad/$norm)")
e_ladlad_norm=$(python -c "print '%.13f' % ($e_ladlad/$norm)")
e_ladlam_norm=$(python -c "print '%.13f' % ($e_ladlam/$norm)")

# Set output file names
sim_name="sene_NLAM_${lam_atoms}_run_${run}"
init_file="init_${sim_name}.in"
restart_file="restart_${sim_name}"
prep1_outfile="prep1_${sim_name}.lammpstrj"
prep2_outfile="prep2_${sim_name}.lammpstrj"
run_outfile="run_${sim_name}.lammpstrj"
pos_file="pos_${sim_name}.dat"
map_file="${sim_name}.lammpsmap"
equil_simfile="equil_${sim_name}.out"
run1_simfile="preInterxn_${sim_name}.out"
run2_simfile="postInterxn_${sim_name}.out"
run3_simfile="end_${sim_name}.out"

# Convert all time values to simulation time (i.e. rescale by delta t)
restart_freq=$(bc <<< "$restart_freq/$delta_t")
prep1_time=$(bc <<< "$prep1_time/$delta_t")
prep1_printfreq=$(bc <<< "$prep1_printfreq/$delta_t")
prep2_time=$(bc <<< "$prep2_time/$delta_t")
prep2_printfreq=$(bc <<< "$prep2_printfreq/$delta_t")
run1_time=$(bc <<< "$run1_time/$delta_t")
run2_time=$(bc <<< "$run2_time/$delta_t")
run3_time=$(bc <<< "$run3_time/$delta_t")
run_printfreq=$(bc <<< "$run_printfreq/$delta_t")

# Make execution directory
run_dir="${run_dir}/${sim_name}"
mkdir $run_dir

# Create the lammps command file based on template
lammps_file="${sim_name}.lam"
file="${run_dir}/${lammps_file}"
cp Senescence.lam $file

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

sed -i -- "s/POS_FILE/${pos_file}/g" $file

sed -i -- "s/RUN_PRINTFREQ/${run_printfreq}/g" $file
sed -i -- "s/RUN_OUTFILE/${run_outfile}/g" $file

sed -i -- "s/RUN1_SEED/${run1_seed}/g" $file
sed -i -- "s/RUN1_TIME/${run1_time}/g" $file

sed -i -- "s/RUN2_SEED/${run2_seed}/g" $file
sed -i -- "s/RUN2_TIME/${run2_time}/g" $file

sed -i -- "s/RUN3_SEED/${run3_seed}/g" $file
sed -i -- "s/RUN3_TIME/${run3_time}/g" $file

sed -i -- "s/LAM_ATOMS/${lam_atoms}/g" $file
sed -i -- "s/LAM_SEED/${lam_seed}/g" $file

sed -i -- "s/DELTA_T/${delta_t}/g" $file

sed -i -- "s/SIGMA/${sigma}/g" $file

sed -i -- "s/EHETHET/${e_hethet_norm}/g" $file
sed -i -- "s/EHETLAM/${e_hetlam_norm}/g" $file
sed -i -- "s/EHETLAD/${e_hetlad_norm}/g" $file
sed -i -- "s/ELADLAD/${e_ladlad_norm}/g" $file
sed -i -- "s/ELADLAM/${e_ladlam_norm}/g" $file

sed -i -- "s/CUTOFF/${cutoff}/g" $file

sed -i -- "s/EQUIL_SIMFILE/${equil_simfile}/g" $file
sed -i -- "s/RUN1_SIMFILE/${run1_simfile}/g" $file
sed -i -- "s/RUN2_SIMFILE/${run2_simfile}/g" $file
sed -i -- "s/RUN3_SIMFILE/${run3_simfile}/g" $file

# Generate chromatin

${gen_chromo_exe} $chromo_file $lam_file $het_file $chr_num $box_size $box_size $box_size $buffer "${run_dir}/${init_file}" "${run_dir}/${map_file}"
