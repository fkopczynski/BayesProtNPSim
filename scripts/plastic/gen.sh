#!/bin/bash

# some defaults
dist=1.5
ion_conc=0.15
ion_repl='SOL'
plast=""
len=0

# get script parameters
while getopts "p:l:d:c:r:" opt; do
case $opt in
	p) plast="$OPTARG" ;;
	l) len=$OPTARG ;;
	d) dist=$OPTARG ;;
	c) ion_conc=$OPTARG ;;
	r) ion_repl=$OPTARG ;;
	?) echo "Invalid option: -$OPTARG" ;;
	:) echo "Arguments -p and -l must be given"; exit 1 ;;
esac
done

echo ${plast}

# check if plastic and its length are given because they are essential, otherwise exit
if [ -z $plast ] || [ $len -eq 0 ]; then
echo "Arguments -p and -l must be given"
exit 1
elif [[ $plast != "ps" && $plast != "idk" && $plast != "pm"]]; then 
echo "Please choose -p between ps, ?, pm"
exit 1
fi

# define standard residues
struct1=$plast
head1="h${plast}"
tail1="t${plast}"
norm1=$plast

# define carboxlated residues
struct2="${plast}_carb"
head2="h${plast:1:1}c"
tail2="t${plast:1:1}c"
norm2="${plast}c"

# activate amber env
source ~/miniconda3/etc/profile.d/conda.sh
conda activate amber

# prepare standard residues
antechamber -fi mol2 -fo ac -i ${struct1}.mol2 -o ${struct1}.ac -c bcc -pf y
prepgen -i ${struct1}.ac -o ${norm1}.prepi -f prepi -m mainchain.${norm1} -rn ${norm1^^} -rf ${norm1}.res
prepgen -i ${struct1}.ac -o ${head1}.prepi -f prepi -m mainchain.${head1} -rn ${head1^^} -rf ${head1}.res
prepgen -i ${struct1}.ac -o ${tail1}.prepi -f prepi -m mainchain.${tail1} -rn ${tail1^^} -rf ${tail1}.res
rm ${struct1}.ac

# prepare carboxylated residues
antechamber -fi mol2 -fo ac -i ${struct2}.mol2 -o ${struct2}.ac -c bcc -pf y
prepgen -i ${struct2}.ac -o ${norm2}.prepi -f prepi -m mainchain.${norm2} -rn ${norm2^^} -rf ${norm2}.res
prepgen -i ${struct2}.ac -o ${head2}.prepi -f prepi -m mainchain.${head2} -rn ${head2^^} -rf ${head2}.res
prepgen -i ${struct2}.ac -o ${tail2}.prepi -f prepi -m mainchain.${tail2} -rn ${tail2^^} -rf ${tail2}.res
rm ${struct2}.ac

# generate the sequence
# stard with the normal head
echo "{ ${head1^^}" > poly

# middle residues
(( mid=(${len}-2) ))
echo ${mid}
for ((res=1; res<=${mid}; res++))
do
if ((${res} % 2 != 0))
then 
echo "${norm2^^}" >> poly
else
echo "${norm1^^}" >> poly
fi
done

# last residue
# include elif for the case then len=2
if [[ $(tail poly -n1) == "${norm1^^}" ]]
then
echo "${tail2^^} }" >> poly
elif [[ $(tail poly -n1) == *"${head1^^}" ]]
then
echo "${tail2^^} }" >> poly
else
echo "${tail1^^} }" >> poly
fi

# make the sequence one line
cat poly | tr "\n" " " > poly.seq
rm poly

cat poly.seq

# generate the tleap.in file
echo "source leaprc.gaff" > tleap.in
for file in ./*.prepi
do
echo "loadamberprep ${file}" >> tleap.in
done
echo "mol = sequence $(cat poly.seq)" >> tleap.in
echo "savepdb mol poly.pdb" >> tleap.in
echo "saveamberparm mol poly.top poly.crd" >> tleap.in
echo "quit" >> tleap.in

# polymerize anionic
tleap -s -f tleap.in > tleap.out

acpype -p poly.top -x poly.crd -b ${struct1}${len}
mv ${struct1}${len}.amb2gmx ${struct1}${len}
cp plastic_prep.sh ${struct1}${len}/.

# deactivate python environment
conda deactivate
cd ${struct1}${len}

# GROMACS PREPARATION PART
# create the directories for future computations
mkdir plastic
mv * plastic/
mkdir prot_pl
cd prot_pl
mkdir md
mkdir md/umbrella

# put all other required files here
#ADD ALSO SOME ANALYSIS SCRIPTS!!!
cp ~/project/protein/* .
cd ../plastic/
cp ~/project/mdps/plastic/* .
cp ~/project/scripts/plastic/* .
cp -r ~/project/charmm27.ff .
mkdir md
mv md.mdp *gpu_submit md/.

# create the topol.top file
acp_top=*GMX.top
echo "; plastic topology" > topol.top
echo '#include "charmm27.ff/forcefield.itp"' >> topol.top
echo "; include plastic parameters" >> topol.top
echo '#include "plastic.itp"' >> topol.top
echo '#include "charmm27.ff/tip3p.itp"' >> topol.top
echo '#include "charmm27.ff/ions.itp"' >> topol.top
tail -n7 $acp_top >> topol.top

# create the plastic.itp file
# find the line number where the section starts
startline=$(awk '/moleculetype/ {print NR}' $acp_top) 
stopline=$(awk '/system/ {print NR}' $acp_top)
awk "NR >= $startline && NR < $stopline" $acp_top > plastic.itp

# put the system in a box
ac_gro=*.gro
gmx editconf -f $ac_gro -bt cubic -d $dist -o box.gro &> prep.log

# solvate
gmx solvate -cs -cp box.gro -p topol.top -o solve.gro &>> prep.log

# energy minimisation (emw)
gmx grompp -f emw.mdp -c solve.gro -p topol.top -o emw.tpr -maxwarn 2 &>> prep.log
mpirun -np $SLURM_NTASKS gmx_mpi mdrun -deffnm emw -v -ntomp $SLURM_CPUS_PER_TASK >& emw.out

#put ion
gmx grompp -f emw.mdp -c emw.gro -p topol.top -o ion.tpr -maxwarn 1 &>> prep.log
echo $ion_repl | gmx genion -s ion.tpr -p topol.top -neutral -conc $ion_conc -o ion.gro &>> prep.log

# em with ions
gmx grompp -f em.mdp -c ion.gro -p topol.top -o em.tpr &>> prep.log
mpirun -np $SLURM_NTASKS gmx_mpi mdrun -deffnm em -v -ntomp $SLURM_CPUS_PER_TASK >& em.out

# water equilibration (position restraint run)
gmx grompp -f posre.mdp -c em.gro -r em.gro -p topol.top -o posre.tpr -maxwarn 3 &>> prep.log
mpirun -np $SLURM_NTASKS gmx_mpi mdrun -deffnm posre -v -ntomp $SLURM_CPUS_PER_TASK >& posre.out

# production run
cd md
gmx grompp -f md.mdp -c ../posre.gro -p ../topol.top -o md.tpr &>> prep.log

# for generating a consistent job name
sed -i "2s/.*/#SBATCH --job-name=${norm1}${len}eq/" *gpu_submit
#sbatch 32gpu_submit

