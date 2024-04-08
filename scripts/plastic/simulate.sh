#!/bin/bash
srun --mpi=pmix gmx_mpi mdrun -deffnm md -v -ntomp $SLURM_CPUS_PER_TASK >& md.out

# extract trajectories
echo -e '"Other"|"Ion" \nq' | srun --mpi=pmix gmx_mpi make_ndx -f md.tpr -o pl_ion.ndx &> post_proc.out
echo -e "Other_Ion \nSystem \nSystem" | srun --mpi=pmix gmx_mpi trjconv -f md.xtc -s md.tpr -center -pbc cluster -n pl_ion.ndx -dt 1000 -o dry_cluster.xtc &>> post_proc.out
echo -e "Other \nOther \nOther_Ion" | srun --mpi=pmix gmx_mpi trjconv -f dry_cluster -s md.tpr -fit rot+trans -center -o dry.xtc -n pl_ion.ndx &>> post_proc.out
echo "Other_Ion" | srun --mpi=pmix gmx_mpi trjconv -f dry_cluster.xtc -s md.tpr -n pl_ion.ndx -dump 0 -o dry.gro &>> post_proc.out

# extract structural parameters
echo 'Other Other' | srun --mpi=pmix gmx_mpi rms -f dry_cluster.xtc -s md.tpr -o rmsd.xvg -xvg none &>> post_proc.out
echo 'Other' | srun --mpi=pmix gmx_mpi gyrate -f dry_cluster.xtc -s md.tpr -o gyr.xvg -xvg none &>> post_proc.out
rm dry_cluster.xtc

# energy-related params
echo 'Potential ' | srun --mpi=pmix gmx_mpi energy -f md.edr -o pe.xvg -xvg none &>> post_proc.out
echo 'Kinetic ' | srun --mpi=pmix gmx_mpi energy -f md.edr -o ke.xvg -xvg none &>> post_proc.out
echo 'Temperature ' | srun --mpi=pmix gmx_mpi energy -f md.edr -o temp.xvg -xvg none &>> post_proc.out
echo 'Pressure ' | srun --mpi=pmix gmx_mpi energy -f md.edr -o pres.xvg -xvg none &>> post_proc.out

# create a final plastic molecule for the next stage
echo 'Other' | srun --mpi=pmix gmx_mpi editconf -f md.gro -n pl_ion.ndx -o plastic.pdb &>> post_proc.out
