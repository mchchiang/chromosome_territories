#!/bin/bash
#
# SGE (Sun Grid Engine) submission script for an array of runs
#
# Usage: qsub -e [error_log] -o [output_log] ./qsub_script.sh [sim_name] [job_in_dir] [job_out_dir]
#
#$ -cwd              # run from current directory
#$ -N sene           # job name
#$ -V                # use all shell environment variables
#
# Choose a queue:
#$ -q cdt.7.day      # cdt.7.day, cm.7.day, sopa.1.day 
#
# Set job runtime
#$ -l h_rt=00:15:00  # time limit = 7days  (in sopa.1.day the time limit = 24h)
#
# Specify the standard output and error log file using -e and -o in qsub
#

sim_name=$1
job_in_dir=$2
job_out_dir=$3

local_dir=$PWD
node_job_dir=/scratch/s1309877_sim/$sim_name/

# Find the remote scratch directory when the job is running
remote_scratch=/Disk/$(hostname -s)_staging
# Replace 'scratch' with the correct remote scratch directory
remote_job_dir=${node_job_dir/\/scratch/$remote_scratch}


# Stagein data
stagein() {
    echo "Staging data from $job_in_dir to $remote_job_dir"
    mkdir -p $remote_job_dir \
	&& rsync -ravg $job_in_dir/* $remote_job_dir
    echo
}

# Stageout data
stageout() {
    echo "Staging data from $remote_job_dir to $job_out_dir"
    mkdir -p $job_out_dir \
	&& rsync -ravg $remote_job_dir/* $job_out_dir \
	&& rm -rf $remote_job_dir
    echo
}

# Actual computation
compute() {
    echo "Running computation"
    cd $node_job_dir
    lmp_serial -screen none -in $node_job_dir/$sim_name.lam -log $node_job_dir/$sim_name.log
    cd $local_dir
    echo
}

###############################################################################

# Run the stagein, compute, and stageout steps

stagein || {
    echo "ERROR! Initial data stagein failed!"
    echo "Aborting job."
    exit 1
}
compute || {
    echo "ERROR! Compute job failed!"
    echo "Job's staging area has been left alone."
    echo "You can access this at: $remote_job_dir"
    exit 1
}
stageout || {
    echo "ERROR! Final data stageout failed!"
    echo "You may need to manually rescue data from job's staging area."
    echo "You can access this at: $remote_job_dir"
    exit 1
}
