#!/usr/bin/env python3

# This script extracts the average force and standard deviation from an umbrella sampling simulation. By default, it starts collecting forces for block analysis at 35 ns and ends at 50 ns. Block size is by default 1 ns.
# Important! This script requires a positional argument, which is the sampled COM separation (and the results folder) name. This is taken care of by the optimiser but if executed manually, this needs to be adjusted. 

import numpy as np
import os
import sys

def extract_force(umb, start=35000, end=50000, size=1000):
    '''
    Calculates the forces for a simulated system. Requirements:
    - folders in which the results are located must be called the same as the CV value at the umbrella, e.g. 2.1: that defines the umbrellas; these folders must contain COLVAR files from PLUMED
    - time resolutino in COLVAR file must be 1 ps so that indexing works properly
    '''
    
    # define the bins / blocks from input parameters, as well as the umbrella centers
    bins = np.arange(start, end, size)

    # initialise numpy arrays that will store n rows of m averaged binned forces (n=number of blocks; m=number of umbrellas to integrate)
    # each row will be an input for integration to get a separate free energy profile
    # std_err is the standard deviation of obtained free energy paths, it will contain n elements
    n_bins = len(bins)-1
    binned_forces = np.zeros(n_bins)

    # open COLVAR file in the umbrella folder
    # extract info about relevant parameters from COLVAR

    with open('COLVAR') as f:
        lines = f.readlines()[1:]

        # initialise the list storing forces from COLVAR
        force_raw = []

        # define the centre of the umbrella and the force constant from the last recorded value
        ref = float(lines[-1].strip().split()[4])
        kappa = float(lines[-1].strip().split()[6])

        # get the remaining params
        for l in lines:
            params = l.strip().split()

            # calculate the force as -k(com-com_ref)
            com_raw = float(params[1])
            force_raw.append(-kappa*(com_raw - ref))
            
        # calculate the average force per bin
        for i in range(n_bins):
            bin_start, bin_end = bins[i], bins[i+1]
            f_bin = np.mean(force_raw[bin_start:bin_end])

            # replace the zeros in the force matrix with the averaged values
            binned_forces[i] = f_bin
        
    os.chdir('../')
    #print(binned_forces)
    # averaged free energy profile
    final_force = np.average(binned_forces)
    final_std = np.std(binned_forces)
    return final_force, final_std

umb = sys.argv[1]
force, std = extract_force(umb=umb)
print(force)
print(std)
