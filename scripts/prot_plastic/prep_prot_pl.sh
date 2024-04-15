#!/bin/bash

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
cp ~/project/mdps/prot_plastic/* .
cp ~/project/scripts/prot_plastic/*gpu_submit md/.
cp ~/project/scripts/prot_plastic/simulate.sh md/. 

# Copy a correct set of input files and the modified forcefield
cp ~/results/plastic/${plast}/${plast}${len}/plastic/plastic.itp .
cp -r ~/project/charmm27.ff .

# change the topol.top file so that the plastic molecule is added
# the plastic needs to be added to the [molecules] section to the last line
echo "$plast$len	1" >> topol.top

# put in a box, solvate
gmx editconf -f conf.pdb -bt cubic -d 2.0 -o box.gro &> prep.log
gmx solvate -cs -cp box.gro -p topol.top -o solve.gro &>> prep.log

# energy minimisation (emw)
gmx grompp -f emw.mdp -c solve.gro -p topol.top -o emw.tpr -maxwarn 4 &>> prep.log
mpirun -np $SLURM_NTASKS gmx_mpi mdrun -deffnm emw -ntomp $SLURM_CPUS_PER_TASK >& emw.out
 
# add ions
gmx grompp -f emw.mdp -c emw.gro -p topol.top -o ion.tpr -maxwarn 3 &>> prep.log
echo $ion_repl  | gmx genion -s ion.tpr -p topol.top -neutral -conc $ion_conc -o ion.gro &>> prep.log

# em with ions
gmx grompp -f em.mdp -c ion.gro -p topol.top -o em.tpr -maxwarn 2 &>> prep.log
mpirun -np $SLURM_NTASKS gmx_mpi mdrun -deffnm em -ntomp $SLURM_CPUS_PER_TASK >& em.out

# water equilibration (position restraint run)
gmx grompp -f posre.mdp -c em.gro -r em.gro -p topol.top -o posre.tpr -maxwarn 5 &>> prep.log
mpirun -np $SLURM_NTASKS gmx_mpi mdrun -deffnm posre -ntomp $SLURM_CPUS_PER_TASK >& posre.out

# final equilibration
gmx grompp -f md.mdp -c posre.gro -p topol.top -o md/md.tpr &>> prep.log
cd md

sed -i "2s/.*/#SBATCH --job-name=${plast}${len}prot/" *gpu_submit
sbatch 16gpu_submit
