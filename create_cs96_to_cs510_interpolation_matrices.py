#!/usr/bin/env python3
''' 
Function to get the nearest neighbour weights to 
interpolate between low and high resolution mitgcm grids
'''
__author__ = "Hamish Hay"

import numpy as np 
from scipy.spatial import cKDTree
from scipy.sparse import coo_matrix
from scipy.sparse import save_npz
from scipy.sparse import identity
from MITgcmutils import mds
from time import time

grid_cs510 = "/home/hay/Research/EuropaOceanDynamics/MITgcm_grids/cs510/"
grid_cs96 = "/home/hay/Research/EuropaOceanDynamics/MITgcm_grids/cs96/"
grid_cs510 = grid_cs96

def sph2xyz(r, theta, phi):
    x = r*np.cos(phi)*np.sin(theta)
    y = r*np.sin(phi)*np.sin(theta)
    z = r*np.cos(theta)

    return x, y, z

radius = 1550e3
d2r = np.pi/180.

def get_nearest_interp_matrix(xName, yName, zName='RC', savename="interp_mat.npz", is_2D=False, p=1.4, neighbour_num=6):
    Ncs96 = 6*96*96
    Ncs510 = 6*510*510
    Ncs510 = Ncs96
    
    X = mds.rdmds(grid_cs96 + xName)
    Y = mds.rdmds(grid_cs96 + yName)
    Z = mds.rdmds(grid_cs96 + zName)
    Z = (Z + radius) / radius

    if is_2D:
        LATr = d2r*(Y+90)
        LONr = d2r*(X+180)
        
        R = 1.0
        

    else:
        R, LONr = np.meshgrid(Z, d2r*(X+180))
        R, LATr = np.meshgrid(Z, d2r*(Y+90))

        R = R.ravel()
        
    LATr = LATr.ravel()
    LONr = LONr.ravel()

    x, y, z = sph2xyz(R, LATr, LONr)

    xyz = np.column_stack((x, y, z))

    tree = cKDTree( xyz, leafsize=100)

    print("GOT TREE!")

    # Load coords of higher res grid
    X = mds.rdmds(grid_cs510 + xName)
    Y = mds.rdmds(grid_cs510 + yName)
    Z = mds.rdmds(grid_cs510 + zName)
    Z = (Z + radius) / radius

    if is_2D:
        LATr = d2r*(Y+90)
        LONr = d2r*(X+180)
        
        R = 1.0

        Nr = 1
        
    else:
        # R, LONr = np.meshgrid(Z, d2r*(X+180))
        # R, LATr = np.meshgrid(Z, d2r*(Y+90))

        LATr = d2r*(Y+90)
        LONr = d2r*(X+180)

        # R = R.ravel()

        R = 1.0

        Nr = 1

    LATr = LATr.ravel()
    LONr = LONr.ravel()

    x2, y2, z2 = sph2xyz(R, LATr, LONr)

    print(x2.size)

    xyz2 = np.column_stack((x2, y2, z2))


    print("QUERYING TREE...")
    s = time()
    distances, ix = tree.query( xyz2, k=neighbour_num, n_jobs=-1)
    e = time()

    print("GOT TREE QUERY AFTER {:1.2f} SECONDS!".format(e-s))
    
    print(distances.shape[0])
    
    print("COMPUTING WEIGHTS")
    w = 1/distances**p
    w /= np.sum(w, axis=-1)[:, np.newaxis]

    print(w.shape, x.size)
    cnt = 0
    for i in range(w.shape[0]):
        if 0.0 in distances[i]:
            # print("Found zero")
            at_p = np.zeros(neighbour_num, dtype=np.float32)
            at_p[0] = 1.0
            w[i] = at_p

            # ix[i, 0] = i
    
    print("CREATING SPARSE INTERPOLATION MATRIX")
    rows = np.zeros(w.size, dtype=np.int)
    cols = np.zeros(w.size, dtype=np.int)
    count = 0
    for ii in range(x2.size):
        for j in range(neighbour_num):
            rows[count] = ii
            cols[count] = ix[ii, j]

            count += 1

    weights_mat = coo_matrix( (w.ravel(), (rows, cols)), shape=(Ncs510*Nr, Ncs96*Nr), dtype=np.float32)
    weights_mat.eliminate_zeros()
    weights_mat = weights_mat.tocsr()


    print("SAVING MATRIX")
    save_npz(savename, weights_mat)


if __name__=='__main__':

    # cell center points at the surface
    get_nearest_interp_matrix('XC', 'YC',       is_2D=True, savename="cs96_to_cs96_eta.npz")

    get_nearest_interp_matrix('XG', 'YC', is_2D=True, savename="cs96_to_cs96_uvel.npz")

    get_nearest_interp_matrix('XC', 'YG', is_2D=True, savename="cs96_to_cs96_vvel.npz")