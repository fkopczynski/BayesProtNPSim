#!/usr/bin/env python3

# this script generates overlayed work plot for each umbrella

import mdtraj as md
import matplotlib.pyplot as plt
import numpy as np
import os

os.chdir('./umbrella')

pl, axs = plt.subplots(1,1, figsize=(12,12))
#axs = ax.flatten()

folders = os.listdir()
lst = sorted(folders, reverse=True)
print(type(folders))

for i,folder in enumerate(lst):
    
    os.chdir(os.getcwd()+'/'+folder)
    print(os.getcwd())
    if 'COLVAR' not in os.listdir():
        os.chdir('../')
        continue

    path = os.getcwd()+'/COLVAR'

    axs.set_title(folder)
    axs.set_xlabel('t / ps')
    with open(path, 'r') as f:
        lines = f.readlines()[1:]
        t = np.zeros(len(lines))
        cv = np.zeros(len(lines))
        steer_cntr = np.zeros(len(lines))
        steer_bias = np.zeros(len(lines))
        steer_work = np.zeros(len(lines))
        for j,line in enumerate(lines):
            values = line.strip().split()
            t[j] = float(values[0])
            cv[j] = float(values[1])
            steer_cntr[j] = float(values[4])
            steer_bias[j] = float(values[2])
            steer_work[j] = float(values[5])
        axs.set_ylabel('COM separation / nm')
        axs.plot(t,cv, label=cv)
        axs.plot(t,steer_cntr, label=folder)
        #axs[3*i+1].set_ylabel('steer_bias')
        #axs[3*i+1].plot(t,steer_bias, label=steer_bias)
        #axs[3*i+2].set_ylabel('steer.d1_work / kcal')
        #axs[3*i+2].plot(t,steer_work, label=steer_work)

        os.chdir('../')
#pl.legend()
pl.show()
pl.savefig('umb_bias.png')


