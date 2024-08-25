#!/bin/bash

# Script to run for manual calculations without the optimiser (e.g. when full free-energy profiles are determined for a plastic length), assumes that the simulation will run for 75 ns and that the appropriate .tpr file with the simulation extension has been generated and names 75ns.tpr
# Yields the trajectory without PBC
# Also gives a post_proc.log file where any potential errors are gathered

# run the simulation
srun --mpi=pmix gmx_mpi mdrun -s 75ns.tpr -cpi md.cpt -plumed plumed.dat -noappend -ntomp $SLURM_CPUS_PER_TASK >& md.out

# convert trajectory to remove PBC, rotations and translations
echo -e '"Protein" | 13 | "Other" \nq ' | srun --mpi=pmix gmx_mpi make_ndx -f 75ns.tpr -o dry.ndx &> post_proc.log
echo "Protein_CAL_Other Protein System" | srun --mpi=pmix gmx_mpi trjconv -f traj_comp.part0002.xtc -s 75ns.tpr -dt 100 -pbc cluster -center -o cluster.xtc -n dry.ndx &>> post_proc.log
echo "Protein Protein_CAL_Other" | srun --mpi=pmix gmx_mpi trjconv -f cluster.xtc -n dry.ndx -s 75ns.tpr -o dry.xtc -fit rot+trans &>> post_proc.log

# create first snapshot for VMD
echo "Protein_CAL_Other " | srun --mpi=pmix gmx_mpi trjconv -f cluster.xtc -s 75ns.tpr -o dry.gro -dump 0 -n dry.ndx &>> post_proc.log

# clean up a bit
rm cluster.xtc
