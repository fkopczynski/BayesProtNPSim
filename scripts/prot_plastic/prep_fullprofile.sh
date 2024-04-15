#!/bin/bash

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
cp ~/project/plumed_files/* .

# set plumed.dat, velocity of pulling set to 0.35 nm/ns
# save calculated params for debugging
vel=0.35
timestep_ns=0.000002
restart_step=50000000
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
