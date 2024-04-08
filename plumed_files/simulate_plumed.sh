#!/bin/bash

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
