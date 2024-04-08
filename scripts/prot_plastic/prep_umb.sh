#!/bin/bash

# copy the tpr file
cp ../md/md.tpr .
cp ../md/md.cpt .

# prepare the tpr file
gmx convert-tpr -s md.tpr -extend 75000 -o 75ns.tpr

for i in $(seq 1.6 0.5 6.6); do
mkdir $i
cd $i
cp ../75ns.tpr .
cp ../md.cpt .
cp ~/project/plumed_files/$i/plumed.dat . 
cp ~/project/plumed_files/*gpu_umbrella .
sed -i "2s/.*/#SBATCH --job-name=u${i}/" 16gpu_umbrella
sed -i "2s/.*/#SBATCH --job-name=u${i}/" 20gpu_umbrella
sed -i "2s/.*/#SBATCH --job-name=u${i}/" 24gpu_umbrella
sed -i "2s/.*/#SBATCH --job-name=u${i}/" 32gpu_umbrella
cd ../
done
rm 75ns.tpr md.cpt
