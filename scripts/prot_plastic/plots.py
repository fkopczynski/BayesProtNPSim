#!/usr/bin/env python3

# libraries
import matplotlib.pyplot as plt
import numpy as np
import mdtraj as md

# time definition, second one for RMSD, gyr and COM because based on xtc
t = np.arange(0, 200002, 2)
t2 = np.arange(0, 200002, 500)

# PE
xvg_pe = open('pe.xvg', 'r')
pe = np.empty(0)

for line in xvg_pe:
    datapoints = line.strip().split()
    pe = np.append(pe, float(datapoints[1]))

xvg_pe.close()

# KE
xvg_ke = open('ke.xvg', 'r')
ke = np.empty(0)

for line in xvg_ke:
    datapoints = line.strip().split()
    ke = np.append(ke, float(datapoints[1]))

xvg_ke.close()

# RMSD protein
xvg_rmsd_prot = open('rmsd_prot.xvg', 'r')
rmsd_prot = np.empty(0)

for line in xvg_rmsd_prot:
    datapoints = line.strip().split()
    rmsd_prot = np.append(rmsd_prot, float(datapoints[1]))

xvg_rmsd_prot.close()

# RMSD plastic
xvg_rmsd_pl = open('rmsd_pl.xvg', 'r')
rmsd_pl = np.empty(0)

for line in xvg_rmsd_pl:
    datapoints = line.strip().split()
    rmsd_pl = np.append(rmsd_pl, float(datapoints[1]))

xvg_rmsd_pl.close()

# radius of gyration
xvg_gyr = open('gyr.xvg', 'r')
gyr = np.empty(0)

for line in xvg_gyr:
    datapoints = line.strip().split()
    gyr = np.append(gyr, float(datapoints[1]))

xvg_gyr.close()

# COM separation
traj = md.load('dry.xtc',top='dry.gro')
ch1 = 'index 1 to 1689'
ch2 = 'index 1690 to 3378'
com1 = md.compute_center_of_mass(traj, select=ch1)
com2 = md.compute_center_of_mass(traj, select=ch2)
com_sep = np.sqrt(np.sum((com2-com1)**2, axis=1))

# Temperature
xvg_temp = open('temp.xvg', 'r')
temp = np.empty(0)

for line in xvg_temp:
    datapoints = line.strip().split()
    temp = np.append(temp, float(datapoints[1]))

xvg_temp.close()

# pressure
xvg_pres = open('pres.xvg', 'r')
pres = np.empty(0)

for line in xvg_pres:
    datapoints = line.strip().split()
    pres = np.append(pres, float(datapoints[1]))

xvg_pres.close()



# plotting everything

fig, axs = plt.subplots(4, 2, figsize=(14,28))
ax = axs.flat

ax[0].set_xlabel('Simulation time / ps')
ax[0].set_ylabel('Potential energy / kJ mol-1')
ax[0].plot(t, pe, color='black', linewidth=0.4) 

ax[1].set_xlabel('Simulation time / ps')
ax[1].set_ylabel('Kinetic energy / kJ mol-1')
ax[1].plot(t, ke, color='black', linewidth=0.4) 

ax[2].set_xlabel('Simulation time / ps')
ax[2].set_ylabel('Temperature / K')
ax[2].plot(t, temp, color='black', linewidth=0.4)

ax[3].set_xlabel('Simulation time / ps')
ax[3].set_ylabel('Pressure / bar')
ax[3].plot(t, pres, color='black', linewidth=0.4)

ax[4].set_xlabel('Simulation time / ps')
ax[4].set_ylabel("Protein C_alpha RMSD / nm")
ax[4].plot(t2, rmsd_prot, color='black', linewidth=0.4)

ax[5].set_xlabel('Simulation time / ps')
ax[5].set_ylabel("Plastic RMSD / nm")
ax[5].plot(t2, rmsd_pl, color='black', linewidth=0.4)

ax[6].set_xlabel('Simulation time / ps')
ax[6].set_ylabel('Radius of gyration / nm')
ax[6].plot(t2, gyr, color='black', linewidth=0.4)

ax[7].set_xlabel('Simulation time / ps')
ax[7].set_ylabel('COM separation / nm')
ax[7].plot(t2, com_sep, color='black', linewidth=0.4)

plt.savefig('param_check.png')
