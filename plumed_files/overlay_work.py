#!/usr/bin/env python3

# This script generates overlayed work and bias plot for each umbrella
# Useful in cases where full free-energy profiles are determined (not relevant for optimiser-determined points)

import matplotlib.pyplot as plt
import numpy as np
import os

os.chdir('./umbrella')

pl, axs = plt.subplots(1,1, figsize=(12,12))

# get the CV folders 
folders = os.listdir()
lst = sorted(folders, reverse=True)

# iterate through different points on the CV, extracting forces and work from COLVAR with respect to time
for i,folder in enumerate(lst):
    
    os.chdir(os.getcwd()+'/'+folder)
    print(os.getcwd())

    # if there is no COLVAR file, skip
    if 'COLVAR' not in os.listdir():
        os.chdir('../')
        continue

    path = os.getcwd()+'/COLVAR'

    axs.set_title(folder)
    axs.set_xlabel('t / ps')

    # open COLVAR
    with open(path, 'r') as f:
        lines = f.readlines()[1:]
        t = np.zeros(len(lines))
        cv = np.zeros(len(lines))

        # initiate the arrays
        steer_cntr = np.zeros(len(lines))
        steer_bias = np.zeros(len(lines))
        steer_work = np.zeros(len(lines))

        # get the values for each line in COLVAR
        for j,line in enumerate(lines):
            values = line.strip().split()
            t[j] = float(values[0])
            cv[j] = float(values[1])
            steer_cntr[j] = float(values[4])
            steer_bias[j] = float(values[2])
            steer_work[j] = float(values[5])

        # plot the time evolution of the parameters
        axs.set_ylabel('COM separation / nm')
        axs.plot(t,cv, label=cv)
        axs.plot(t,steer_cntr, label=folder)

        os.chdir('../')

pl.savefig('umb_bias.png')


