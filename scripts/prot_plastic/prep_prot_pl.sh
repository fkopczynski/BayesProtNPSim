#!/bin/bash

# Script that prepares and runs protein equilibration with the plastic.
# Runs in the optimiser, can also be used for manual simulations (as an alternative to simulate.sh)
# Remember to change the path to your cloned repo!
# Parameters:
# -p: plastic name, e.g. ps [REQUIRED]
# -l: plastic length [REQUIRED]
# -d: distance between the molecule and the edge of the box, 1.5 nm by default
# -c: concentration of NaCl, by default 150 mM
# -r: molecules which are replaced with ions during GROMACS pre-processing, by default SOL
# this sctipt generated a prep.log file which contains the outputs and potential error messages of GROMACS commands

# path to the cloned repo
# change according to your placement of the folder!
proj_path="$HOME/project"

# set up gromacs executable
gmx_mpi="/home/spack-user/spack/opt/spack/linux-centos7-zen3/aocc-3.1.0/gromacs-2020.4-z7lmmyeup2uhxfy2mr3bwi2dt6k4grzy/bin/gmx_mpi"

# some defaults
ion_conc=0.15
ion_repl='SOL'
plast=""
len=0

# get script parameters
while getopts "p:l:c:r:" opt; do
case $opt in
	p) plast="$OPTARG" ;; 
	l) len=$OPTARG ;;
	c) ion_conc=$OPTARG ;;
	r) ion_repl=$OPTARG ;;
	?) echo "Invalid option: -$OPTARG" ;;
	:) echo "Arguments -p and -l must be given"; exit 1 ;;
esac
done


# check if plastic and its length are given, otherwise exit
if [ -z $plast ] || [ $len -eq 0 ]; then
echo "Arguments -p and -l must be given"
exit 1
fi

# set up the structure
cp $proj_path/mdps/prot_plastic/* .

# Copy a correct set of input files and the modified forcefield
cp ../plastic/plastic.itp .
cp -r $proj_path/charmm27.ff .

# change the topol.top file so that the plastic molecule is added
# the plastic needs to be added to the [molecules] section to the last line
echo "$plast$len	1" >> topol.top

# put in a box, solvate
srun --mpi=pmix $gmx_mpi editconf -f conf.pdb -bt cubic -d 2.0 -o box.gro &> prep.log
srun --mpi=pmix $gmx_mpi solvate -cs -cp box.gro -p topol.top -o solve.gro &>> prep.log

# energy minimisation (emw)
srun --mpi=pmix $gmx_mpi grompp -f emw.mdp -c solve.gro -p topol.top -o emw.tpr -maxwarn 4 &>> prep.log
srun --mpi=pmix $gmx_mpi mdrun -deffnm emw >& emw.out
 
# add ions
srun --mpi=pmix $gmx_mpi grompp -f emw.mdp -c emw.gro -p topol.top -o ion.tpr -maxwarn 3 &>> prep.log
echo $ion_repl  | srun --mpi=pmix $gmx_mpi genion -s ion.tpr -p topol.top -neutral -conc $ion_conc -o ion.gro &>> prep.log

# em with ions
srun --mpi=pmix $gmx_mpi grompp -f em.mdp -c ion.gro -p topol.top -o em.tpr -maxwarn 2 &>> prep.log
srun --mpi=pmix $gmx_mpi mdrun -deffnm em >& em.out

# water equilibration (position restraint run)
srun --mpi=pmix $gmx_mpi grompp -f posre.mdp -c em.gro -r em.gro -p topol.top -o posre.tpr -maxwarn 5 &>> prep.log
srun --mpi=pmix $gmx_mpi mdrun -deffnm posre >& posre.out

# final posre
srun --mpi=pmix $gmx_mpi grompp -f md.mdp -c posre.gro -p topol.top -o md/md.tpr &>> prep.log
cd md

# production run
srun --mpi=pmix $gmx_mpi mdrun -deffnm md >& md.out
