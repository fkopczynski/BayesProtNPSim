#!/bin/bash

# run the simulation
srun --mpi=pmix gmx_mpi mdrun -deffnm md >& md.out

# convert trajectory to remove PBC, rotations and translations
echo -e '"Protein" | 13 | "Other" \nq ' | srun --mpi=pmix gmx_mpi make_ndx -f md.tpr -o dry.ndx &> post_proc.log
echo "Protein_CAL_Other Protein System" | srun --mpi=pmix gmx_mpi trjconv -f md.xtc -s md.tpr -dt 100 -pbc cluster -center -o cluster.xtc -n dry.ndx &>> post_proc.log
echo "Protein Protein_CAL_Other" | srun --mpi=pmix gmx_mpi trjconv -f cluster.xtc -n dry.ndx -s md.tpr -o dry.xtc -fit rot+trans &>> post_proc.log

# create first snapshot for VMD
echo "Protein_CAL_Other " | srun --mpi=pmix gmx_mpi trjconv -f cluster.xtc -s md.tpr -o dry.gro -dump 0 -n dry.ndx &>> post_proc.log

# calculate RMSD and RMSF
echo -e '3\n3' | srun --mpi=pmix gmx_mpi rms -f cluster.xtc -s md.tpr -o rmsd_prot.xvg -xvg none &>> post_proc.log
echo -e 'Other\nOther' | srun --mpi=pmix gmx_mpi rms -f cluster.xtc -s md.tpr -o rmsd_pl.xvg -xvg none &>> post_proc.log
echo 'Protein' | srun --mpi=pmix gmx_mpi rmsf -f cluster.xtc -s md.tpr -o rmsf_prot.xvg -xvg none &>> post_proc.log
echo 'Other' | srun --mpi=pmix gmx_mpi rmsf -f md.xtc -s md.tpr -o rmsf_pl.xvg -xvg none &>> post_proc.log

# calculate radius of gyration
echo -e 'Protein\n' | srun --mpi=pmix gmx_mpi gyrate -f cluster.xtc -s md.tpr -o gyr.xvg -xvg none &>> post_proc.log

# calculate p, T, pe, ke
echo 'Pressure' | srun --mpi=pmix gmx_mpi energy -f md.edr -o pres.xvg -xvg none &>> post_proc.log
echo 'Temperature' | srun --mpi=pmix gmx_mpi energy -f md.edr -o temp.xvg -xvg none &>> post_proc.log
echo 'Potential' | srun --mpi=pmix gmx_mpi energy -f md.edr -o pe.xvg -xvg none &>> post_proc.log
echo 'Kinetic-En.' | srun --mpi=pmix gmx_mpi energy -f md.edr -o ke.xvg -xvg none &>> post_proc.log

rm cluster.xtc
