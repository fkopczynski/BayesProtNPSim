#!/usr/bin/env python3

# The main optimiser script. It carries out the whole process automatically, i.e.:
# - selects the next point to sample, based on currenlty sampled points
# - generates the plastic structure of a desired length (or skips this step, if it has already been generated)
# - equilibrates the plastic
# - puts equilibrated plastic together with the protein
# - equilibrates the protein with the plastic
# - carries out umbrella sampling at the CV selected by the Bayesian optimiser
# - extracts the force from the simulation, using block analysis
# - returns the force value to the Bayesian optimiser
# - the cycle is configured to run 3 times - after it is done, you can restart the cycle to get additional 3 points.
# The results will be in a folder dedicated to a certain plastic length, they are divided into plastic-only and plastic+protein directory.
# At the end of each datapoint collection, the optimiser state is saved as optstate_x, where x is the number of colleted points in total.
# Put this script in your results directory and specify the path to your cloned repo (both here, as well as in the individual scripts inside the cloned repo - proj_path variable, available to you at the top of all scripts that need it).
# Please remember to also put the input files with the plastic structure in the same directory (see: gen.sh script for more info and examples).

# path to cloned repo
proj_path='$HOME/project'

# import essential libraries
import numpy as np
import matplotlib.pyplot as plt
import os
from skopt.space import Real, Integer
from skopt import Optimizer
from skopt.utils import dump, load

# define the main foler where the main script is executed
main_folder = os.getcwd()

# cast all the current optimizers in a list, sort versions numerically to determine the newest one to load
current_optimizers = [int(opt.split('_')[1]) for opt in os.listdir() if opt.startswith('optstate')]
sorted_optimizers = sorted(current_optimizers)

try:
    # try to load the most recent optimizer
    ver = sorted_optimizers[-1]    
    optimizer = load('optstate_'+str(ver))
    print(f"Loading version {ver}...")

except IndexError:
    # if it occurs, it means that it is the initial run and the optimizer must be defined
    print('There were no optimizers before, initialising...')
    ver = 0
    # define the dimensions
    dim1 = Integer(name='Polymer length', low=0, high=40)
    dim2 = Real(name='CV', low=2.1, high=5.6)

    # define the optimizer
    optimizer = Optimizer(dimensions=[dim1, dim2],
        base_estimator="gp",
        n_random_starts=0,
        n_initial_points=0,
        n_jobs=1,
        acq_func="EI",
        acq_optimizer="lbfgs",
        random_state=1999,
        acq_func_kwargs={"xi": 10000})

    # load data onto the optimizer
    ini_points = [[0,2.1],
                [0,2.6],
                [0,3.1],
                [0,3.6],
                [0,4.1],
                [0,4.6],
                [0,5.1],
                [0,5.6],
                [10,2.1],
                [10,2.6],
                [10,3.1],
                [10,3.6],
                [10,4.1],
                [10,4.6],
                [10,5.1],
                [10,5.6],
                [20,2.1],
                [20,2.6],
                [20,3.1],
                [20,3.6],
                [20,4.1],
                [20,4.6],
                [20,5.1],
                [20,5.6],
                [40,2.1],
                [40,2.6],
                [40,3.1],
                [40,3.6],
                [40,4.1],
                [40,4.6],
                [40,5.1],
                [40,5.6]]

    ini_data=[14.579715052631622, 
            19.190460644736888, 
            10.33424946052636, 
            14.856456881578996, 
            21.62012953947351, 
            13.628588289473507, 
            0.012775631578769259, 
            0.3406779078945591, 
            16.960575907894782, 
            21.57553236842109, 
            48.089532302631625, 
            -3.7394148421052193, 
            0.7589639605261377, 
            -0.356555881579125, 
            0.9089245263156126, 
            0.22910882894719065, 
            15.544790342105308, 
            19.071600368421098, 
            11.658016315789519, 
            6.284280171052677, 
            11.020756381578769, 
            -1.3047369078949158, 
            25.45458977631561, 
            0.0019442236840332768, 
            11.282320578947411, 
            18.108763013157937, 
            6.10140898684215, 
            13.511121144736885, 
            7.001024342105086, 
            1.3792392631577166, 
            0.01062426315771658, 
            2.9354461578945594]

    # tell the optimiser the initial points
    optimizer.tell(x=ini_points,y=ini_data)

    # if this is the initial run, gen.sh script needs to be copied to the same directory
    gen_path = proj_path+'/scripts/plastic/gen.sh'
    cmd = f"cp ${gen_path} ."
    os.system(cmd)


# start the datapoints collection
for i in range(3):
    # print out already collected points
    checked_lengths = set([x[0] for x in optimizer.Xi])
    p = optimizer.Xi
    cmd = f'echo "Checked points are: {p}" > POINTS'
    os.system(cmd)

    # ask the optimiser for the next point
    x = optimizer.ask()
    cmd = f'echo "THE FOLLOWING POINT HAS BEEN SELECTED: {x}" >> POINTS'
    os.system(cmd)

    # check what plastic length needs to be generated and equilibrated
    # check if plastic was already generated; if no then do the simulation
    if x[0] in checked_lengths:
        print(f"length {x[0]} already checked before, proceeding to umbrella...")
        os.chdir(main_folder)

    # if not true, simulation needs to be carried out: generation and equilibration taken care of by gen.sh script
    # then, the protein and plastic are put together and equilibrated
    else:
        print(f"point {x[0]} will be checked!")
        cmd = f"./gen.sh -p ps -l {x[0]} &> gen{x[0]}.log"
        print("Generating plastic structure...")
        os.system(cmd)
        print("Plastic structure generated and equilibrated!")
        cmd = f"cp ps{x[0]}/plastic/md/plastic.pdb ps{x[0]}/prot_pl/"
        os.system(cmd)
        print("Plastic structure copied")
        print(os.getcwd())
        # go to the prot_pl folder
        prot_pl = f"ps{x[0]}/prot_pl/"
        os.chdir(prot_pl)
        os.system('./put_together.py')
        print('Structures combined!')
        print(f'Preparation + simulation of ps{x[0]}...')
        cmd = f"./prep_prot_pl.sh -p ps -l {x[0]}"
        os.system(cmd)
        print('Equilibration complete!')
        os.chdir(main_folder)

    # start biased sampling run in plumed folder
    fold = f"ps{x[0]}/prot_pl/plumed/"
    os.chdir(fold)

    # approximate COM separation to 2 decimal places
    com = round(x[1],2)

    # prepare the system and run the simulation
    cmd = f"./prep_umb.sh {com}"
    print('Preparing for umbrella...')
    os.system(cmd)

    # extract the force from the biased sampling run
    path = f"{com}/"
    os.chdir(path)
    os.system('cp ../get_force.py .')

    cmd = f"./get_force.py {com} > FINAL_FORCE"
    print("Extracting the force...")
    y = os.system(cmd)

    with open('FINAL_FORCE') as f:
        y = float(f.readline().strip())

    print('THE FINAL VALUE OF THE FORCE:')
    print(y)
    optimizer.tell(x, y)

    # tell the solution and save the updated optimizer
    os.chdir(main_folder)
    ver = ver+1
    print(f"Saving version{ver}...")
    dump(res=optimizer, filename=f"optstate_{ver}")
