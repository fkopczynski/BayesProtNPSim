#!/usr/bin/env python3

# This script puts the plastic and the protein together in an optimal non-overlapping position, applying translations and random rotations of the plastic.

# import all the libraries
import mdtraj as md
import numpy as np
import os

def move_plastic(prot_struct="protein.pdb", pl_struct="plastic.pdb", output="moved_plastic.pdb", length=3):

    ''' 
    Translates the plastic so that it is positioned near the protein and the binding spot within a specified distance from the protein. 
    - prot_struct, pl_struct: names of the protein and plastic coordinate files (allowed format: .pdb / .gro), respectively. They must be placed in the same directory.
    - output: the name of the moved output structure (.pdb).
    - length: distance (in nm) from the geometric centres of proteins.
    '''
    
    # load the structures
    protein = md.load(prot_struct)
    plastic = md.load(pl_struct)
    top_protein = protein.top
    chainA = top_protein.select("chainid == 0")
    chainB = top_protein.select("chainid == 1")

    # compute the geometric centres
    center_plastic = md.compute_center_of_geometry(plastic)
    centerA = md.compute_center_of_geometry(protein.atom_slice(chainA))
    centerB = md.compute_center_of_geometry(protein.atom_slice(chainB))
    
    # vectors defining the plane above which the plastic should be placed
    end1 = top_protein.select("resid 17 and name CA")
    coord1 = protein.xyz[-1:, end1, :]
    end2 = top_protein.select("resid 4 and name CA")
    coord2 = protein.xyz[-1:, end2, :]
    prot_vec1 = centerA - centerB
    prot_vec2 = coord2 - coord1

    # cross product with z unit vector to define a perpendicular vector, then make it an appropriate length
    perp = np.cross(prot_vec1, prot_vec2)
    perp_scaled = (perp * length) / np.linalg.norm(perp)

    # define the coordinates of the new center of geometry of the plastic 
    new_center = centerB + perp_scaled

    # calculate by how much the plastic center of geometry needs to move to get to the new one defined by the cross product
    move_vec = new_center - center_plastic

    # move plastic
    plastic.xyz += move_vec
    plastic.save_pdb(output)

    return perp, perp_scaled, move_vec

def rotate_plastic(pl_struct="moved_plastic.pdb", axis="x", angle=90, output="rotated.pdb"):
    '''
    Rotates the plastic along a specified axis (x, y, z) by a specified angle in degrees.
    '''
    # load everything
    plastic = md.load(pl_struct)
    angle_rad = np.radians(angle)

    # bring the geometrical center of the plastic molecule to the origin so that the rotations don't happen together with large translations
    center_plastic = md.compute_center_of_geometry(plastic)
    plastic.xyz -= center_plastic
    
    # define rotation matrices
    if axis == "x":
        rot_mat = np.array([[1, 0, 0], [0, np.cos(angle_rad), -np.sin(angle_rad)], [0, np.sin(angle_rad), np.cos(angle_rad)]])
    elif axis == "y":
        rot_mat = np.array([[np.cos(angle_rad), 0, np.sin(angle_rad)], [0, 1, 0], [-np.sin(angle_rad), 0, np.cos(angle_rad)]])
    elif axis == "z":
        rot_mat = np.array([[np.cos(angle_rad), -np.sin(angle_rad), 0], [np.sin(angle_rad), np.cos(angle_rad), 0], [0, 0, 1]])
    else:
        print("Wrong axis specified: x, y or z.")
        exit(1)

    # iterate over all coordinates and multiply them by an appropriate rotation matrix
    rot_coord = np.zeros(np.shape(plastic.xyz))
    for i, coord in enumerate(plastic.xyz[0]):
        rot_coord[0][i] = np.matmul(rot_mat, coord)

    # replace the coordinates with the rotated ones and move the geometrical centre to the same point as before
    plastic.xyz = rot_coord
    plastic.xyz += center_plastic
    plastic.save_pdb(output)

def combine_pdbs(prot="protein.pdb", plast="rot_xyz.pdb", out="conf.pdb"):
    '''Combines the protein and plastic pdbs. Puts the protein coordinates first in the output pdb, then the plastic.'''
    
    with open(prot) as protein:
        lines_protein = protein.readlines()

    with open(plast) as plastic:
        lines_plastic = plastic.readlines()

    new_pdb = open(out, "w")

    for line in lines_protein[:-1]:
        new_pdb.write(line)
    protein.close()

    for line in lines_plastic[3:]:
        new_pdb.write(line)

    plastic.close()
    new_pdb.close()
        
def check_overlap(struct="conf.pdb", cutoff_d=0.5):
    '''
    This function checks if the structures overlap, it helps to arrange the molecules in an iterative way. If molecules overlap, returns True; otherwise False.
    '''
    # load everything
    comb_struct = md.load(struct)
    top = comb_struct.top

    # select protein and plastic
    prot = top.select('protein and symbol !="H"')
    plast = top.select('not protein')

    # calculate neighbours
    neigh = md.compute_neighbors(comb_struct, query_indices=plast, haystack_indices=prot, cutoff=cutoff_d)
    if neigh[0].size == 0:
        return False
    else:
        return True
    
# copy the plastic pdb from the equilibration
os.system('cp ../plastic/md/plastic.pdb .')

# iterate through different distances of the plastic geometrical centre and COM of the protein, starting from 2.5 nm and slowly increasing by 0.001 nm
length = 2.5
while True:
    move_plastic(length=length)
    rand_ang = np.random.randint(0,360)
    rotate_plastic(pl_struct="moved_plastic.pdb",axis="x",angle=rand_ang,output="rot_x.pdb")
    rand_ang = np.random.randint(0,360)
    rotate_plastic(pl_struct="rot_x.pdb",axis="y",angle=rand_ang,output="rot_xy.pdb")
    rand_ang = np.random.randint(0,360)
    rotate_plastic(pl_struct="rot_xy.pdb",axis="z",angle=rand_ang,output="rot_xyz.pdb")
    combine_pdbs()

    # check for overlap
    if check_overlap(cutoff_d=0.5):

        # try again, as for more linear molecules we can end up with unfortunate rotations for an optimal distance
        rand_ang = np.random.randint(0,360)
        rotate_plastic(pl_struct="moved_plastic.pdb",axis="x",angle=rand_ang,output="rot_x.pdb")
        rand_ang = np.random.randint(0,360)
        rotate_plastic(pl_struct="rot_x.pdb",axis="y",angle=rand_ang,output="rot_xy.pdb")
        rand_ang = np.random.randint(0,360)
        rotate_plastic(pl_struct="rot_xy.pdb",axis="z",angle=rand_ang,output="rot_xyz.pdb")
        combine_pdbs()
        if check_overlap(cutoff_d=0.5):
            length += 0.001
        else:
            break

    else:
        break
