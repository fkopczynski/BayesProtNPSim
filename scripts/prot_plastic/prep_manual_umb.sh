#!/bin/bash

# This script prepares the equilibrated protein+plastic system for a series of 12x 75 ns umbrella sampling runs to generate a full free-energy profile for a given system
# It is not included in the optimiser and should be used manually.
# Important! Remember to change the project path, otherwise required wiles will not be copied
# Time-related parameters will be saved in params.log file

# path to cloned repo
proj_path="$HOME/project"

# pulling velocity
vel=0.35

# timestep in ns
timestep_ns=0.000002

# step at which equilibration is completed and restarted for umbrella sampling
restart_step=50000000

# copy the tpr file
cp ../md/md.tpr .
cp ../md/md.cpt .

# prepare the tpr file
gmx convert-tpr -s md.tpr -extend 75000 -o 75ns.tpr

# loop over umbrellas and prepare each system
for i in $(seq 1.6 0.5 6.6); do
mkdir $i
cd $i
cp ../75ns.tpr .
cp ../md.cpt . 
cp $proj_path/plumed_files/* .

# set plumed.dat
# save calculated params for debugging
dist=$(echo "${i} - 1.6" | bc)
echo "distance to pull: ${dist}" > params.log
time_ns=$(echo "scale=3; ${dist}/${vel}" | bc)
echo "time in ns: ${time_ns}" >> params.log
step=$(echo "${restart_step} + ${time_ns}/${timestep_ns}" | bc)
echo "step: ${step}" >> params.log
sed -i "15s/.*/   STEP1=${step}  AT1=${i}/" plumed.dat

# set job names
sed -i "2s/.*/#SBATCH --job-name=u${i}/" *gpu_umbrella
cd ../
done

# clean up
rm 75ns.tpr md.cpt md.tpr
