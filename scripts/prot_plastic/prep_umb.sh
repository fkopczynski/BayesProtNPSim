#!/bin/bash

# Script preparing the system for umbrella sampling and running it
# Runs in the optimiser
# Assumes that umbrella sampling is intended to last 50 ns
# Requires one command line parameter which defines the value of COM separation to sample, taken care of by the optimiser
# Important! Change the path to the cloned repo

# path to cloned repo
proj_path="$HOME/project"

# pulling velocity
vel=0.35

# timestep in ns
timestep_ns=0.000002

# step at which equilibration is completed and restarted for umbrella sampling
restart_step=50000000

# copy the tpr and cpt file
cp ../md/md.tpr .
cp ../md/md.cpt .

# load the gromacs executable
gmx_mpi="/home/spack-user/spack/opt/spack/linux-centos7-zen3/aocc-3.1.0/gromacs-2020.4-z7lmmyeup2uhxfy2mr3bwi2dt6k4grzy/bin/gmx_mpi"

# prepare the tpr file, extending the simulation by 50 ns
srun --mpi=pmix $gmx_mpi convert-tpr -s md.tpr -extend 50000 -o 50ns.tpr

# prepare the results folder, copy required files
mkdir $1
cd $1
cp ../50ns.tpr .
cp ../md.cpt . 
cp $proj_path/plumed_files/* .

# set plumed.dat
# save calculated params for debugging
dist=$(echo "$1 - 1.6" | bc)
echo "distance to pull: ${dist}" > params.log
time_ns=$(echo "scale=3; ${dist}/${vel}" | bc)
echo "time in ns: ${time_ns}" >> params.log
step=$(echo "${restart_step} + ${time_ns}/${timestep_ns}" | bc)
echo "step: ${step}" >> params.log
sed -i "15s/.*/   STEP1=${step}  AT1=$1/" plumed.dat

# perform the biased simulation
srun --mpi=pmix $gmx_mpi mdrun -s 50ns.tpr -cpi md.cpt -plumed plumed.dat -noappend >& plumed.out

# clean up
rm 50ns.tpr md.cpt md.tpr
