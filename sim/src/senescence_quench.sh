#!/bin/bash

# An input script for setting parameters in the LAMMPS config file

if [ "$#" != 6 ]; then
    echo "Usage: senescence_quench.sh [e_hethet] [e_hetlam] [e_hethet2] [e_hetlam2] [run] [run_dir]"
    exit 1
fi

# Read in and set parameters
e_hethet=$1       # HH interaction in growing
e_hetlam=$2       # HL interaction in growing
e_hethet2=$3      # HH interaction in senescence
e_hetlam2=$4      # HL interaction in senescence
run=$5            # trial number
run_dir=$6        # run directory

# Interaction energy for euchromatin and centromere region
e_eueu=0.2
e_cent=1.0

# Init config (default is random walk)
gen_chromo_exe="../bin/exe/Gen_Chr_Het"
lammps_script="senescence_quench.lam"

chromo_file="../../data/chromo_length.dat"
lam_file="../../data/tamir_data/LaminB1.no-4OHT.Rep1.domains.hg19.nopatch.bed"
het_file="../../data/tamir_data/H3K9me3.Pk.full.hg19.nopatch.bed"
chr_num=20

init_box_size=100
box_size=35
buffer=5 

# A function for generating random numbers
function get_rand(){
    rand=$(python -c "import random, sys; print random.randint(0,$max_seed)")
    echo $rand
}

# Initial box size
ilo=$(python -c "print -$init_box_size/2.0")
ihi=$(python -c "print $init_box_size/2.0")

# Final box size
lo=$(python -c "print -$box_size/2.0")
hi=$(python -c "print $box_size/2.0")

max_seed=100000

restart_freq=1000

# Setting up the chromatin fibre
prep1_printfreq=1000
prep1_seed=$(get_rand)
prep1_time_1=4000 # 4000
prep1_time_2=2000 # 2000
prep1_time_3=4000 # 4000

# Reducing the box size
prep2_printfreq=1000
prep2_seed=$(get_rand)
prep2_time=5000 # 5000

# Equilibrating with the lamina beads
prep3_printfreq=1000
prep3_seed=$(get_rand)
prep3_time=5000 # 5000

# Main simulation
run_printfreq=1000
run_seed=$(get_rand)
run_time=200000

# Timesteps for the quench
run2_time_1=1000 # 1000
run2_time_2=30000 # 30000
run2_printfreq=10 # 10

lam_atoms=2000
lam_seed=$(get_rand)

delta_t=0.01       # time step size in Brownian time units

# Interaction energy parameters
# Harmonic potential
e_harm=100.0

# Soft potential
e_soft=100.0

# LJ potentials
sigma=1.0
cutoff=$(python -c "print '%.13f' % (1.8*$sigma)")
# See input arguments for HH and HL

# Normalisation (ensure minimum of potential is actually epsilon)
norm=$(python -c "print '%.13f' % (1.0 + 4.0*(($sigma/$cutoff)**12-($sigma/$cutoff)**6))")

e_hethet_norm=$(python -c "print '%.13f' % ($e_hethet/$norm)")
e_hetlam_norm=$(python -c "print '%.13f' % ($e_hetlam/$norm)")
e_eueu_norm=$(python -c "print '%.13f' % ($e_eueu/$norm)")
e_cent_norm=$(python -c "print '%.13f' % ($e_cent/$norm)")

e_hethet2_norm=$(python -c "print '%.13f' % ($e_hethet2/$norm)")
e_hetlam2_norm=$(python -c "print '%.13f' % ($e_hetlam2/$norm)")

# Set output file names
sim_name="sene_chr_${chr_num}_L_${box_size}_HH_${e_hethet}_HL_${e_hetlam}_run_${run}"
init_file="init_${sim_name}.in"
restart_file="restart_${sim_name}"
prep1_outfile="prep1_${sim_name}.lammpstrj"
prep2_outfile="prep2_${sim_name}.lammpstrj"
prep3_outfile="prep3_${sim_name}.lammpstrj"
run_outfile="run_${sim_name}.lammpstrj"
run2_outfile="quench_${sim_name}.lammpstrj"
pos_file="pos_${sim_name}.dat"
map_file="${sim_name}.lammpsmap"
equil_simfile="equil_${sim_name}.out"
run_simfile="run_${sim_name}.out"
run2_simfile="end_${sim_name}.out"

# Convert all time values to simulation time (i.e. rescale by delta t)
restart_freq=$(bc <<< "$restart_freq/$delta_t")
prep1_time_1=$(bc <<< "$prep1_time_1/$delta_t")
prep1_time_2=$(bc <<< "$prep1_time_2/$delta_t")
prep1_time_3=$(bc <<< "$prep1_time_3/$delta_t")
prep1_printfreq=$(bc <<< "$prep1_printfreq/$delta_t")
prep2_time=$(bc <<< "$prep2_time/$delta_t")
prep2_printfreq=$(bc <<< "$prep2_printfreq/$delta_t")
prep3_time=$(bc <<< "$prep3_time/$delta_t")
prep3_printfreq=$(bc <<< "$prep3_printfreq/$delta_t")
run_time=$(bc <<< "$run_time/$delta_t")
run_printfreq=$(bc <<< "$run_printfreq/$delta_t")
run2_time_1=$(bc <<< "$run2_time_1/$delta_t")
run2_time_2=$(bc <<< "$run2_time_2/$delta_t")
run2_printfreq=$(bc <<< "$run2_printfreq/$delta_t")

# Make execution directory
run_dir="${run_dir}/${sim_name}"

if [ ! -d $run_dir ]; then
mkdir -p $run_dir
fi

# Create the lammps command file based on template
lammps_file="${sim_name}.lam"
file="${run_dir}/${lammps_file}"

# Copy template lammps script
echo "Copying senescence.lam lammps script"
cp $lammps_script $file

# Replace macros in template with input values
sed -i -- "s/INIT_FILE/${init_file}/g" $file
sed -i -- "s/RESTART_FILE/${restart_file}/g" $file

sed -i -- "s/IXLO/${ilo}/g" $file
sed -i -- "s/IXHI/${ihi}/g" $file
sed -i -- "s/IYLO/${ilo}/g" $file
sed -i -- "s/IYHI/${ihi}/g" $file
sed -i -- "s/IZLO/${ilo}/g" $file
sed -i -- "s/IZHI/${ihi}/g" $file

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
sed -i -- "s/PREP1_TIME_1/${prep1_time_1}/g" $file
sed -i -- "s/PREP1_TIME_2/${prep1_time_2}/g" $file
sed -i -- "s/PREP1_TIME_3/${prep1_time_3}/g" $file

sed -i -- "s/PREP2_PRINTFREQ/${prep2_printfreq}/g" $file
sed -i -- "s/PREP2_OUTFILE/${prep2_outfile}/g" $file
sed -i -- "s/PREP2_SEED/${prep2_seed}/g" $file
sed -i -- "s/PREP2_TIME/${prep2_time}/g" $file

sed -i -- "s/PREP3_PRINTFREQ/${prep3_printfreq}/g" $file
sed -i -- "s/PREP3_OUTFILE/${prep3_outfile}/g" $file
sed -i -- "s/PREP3_SEED/${prep3_seed}/g" $file
sed -i -- "s/PREP3_TIME/${prep3_time}/g" $file

sed -i -- "s/POS_FILE/${pos_file}/g" $file

sed -i -- "s/RUN_PRINTFREQ/${run_printfreq}/g" $file
sed -i -- "s/RUN_OUTFILE/${run_outfile}/g" $file

sed -i -- "s/RUN2_PRINTFREQ/${run2_printfreq}/g" $file
sed -i -- "s/RUN2_OUTFILE/${run2_outfile}/g" $file

sed -i -- "s/RUN_SEED/${run_seed}/g" $file
sed -i -- "s/RUN_TIME/${run_time}/g" $file

sed -i -- "s/RUN2_TIME_1/${run2_time_1}/g" $file
sed -i -- "s/RUN2_TIME_2/${run2_time_2}/g" $file

sed -i -- "s/LAM_ATOMS/${lam_atoms}/g" $file
sed -i -- "s/LAM_SEED/${lam_seed}/g" $file

sed -i -- "s/DELTA_T/${delta_t}/g" $file

sed -i -- "s/EHARM/${e_harm}/g" $file
sed -i -- "s/ESOFT/${e_soft}/g" $file

sed -i -- "s/EHETHET2/${e_hethet2_norm}/g" $file
sed -i -- "s/EHETLAM2/${e_hetlam2_norm}/g" $file

sed -i -- "s/EHETHET/${e_hethet_norm}/g" $file
sed -i -- "s/EHETLAM/${e_hetlam_norm}/g" $file
sed -i -- "s/EEUEU/${e_eueu_norm}/g" $file
sed -i -- "s/ECENT/${e_cent_norm}/g" $file

sed -i -- "s/SIGMA/${sigma}/g" $file

sed -i -- "s/CUTOFF/${cutoff}/g" $file

sed -i -- "s/EQUIL_SIMFILE/${equil_simfile}/g" $file
sed -i -- "s/RUN_SIMFILE/${run_simfile}/g" $file
sed -i -- "s/RUN2_SIMFILE/${run2_simfile}/g" $file


# Generate chromatin

${gen_chromo_exe} $chromo_file $lam_file $het_file $chr_num $init_box_size $init_box_size $init_box_size $buffer "${run_dir}/${init_file}" "${run_dir}/${map_file}"

# Relabel centromere region (as heterochromatin)
awk '{if (NF==9&&$1>=2640&&$1<=2940) {$3=2; print} else {print}}' ${run_dir}/${init_file} > ${run_dir}/temp
mv ${run_dir}/temp ${run_dir}/${init_file}
