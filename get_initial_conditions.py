#!/usr/bin/env python3
''' 
Function to get the MITgcm output at a specific time,
interpolate it from cs96 to cs510, and save the 
corresponding initial condition binary file.
'''
__author__ = "Hamish Hay"

import numpy as np 
from scipy.sparse import coo_matrix
from scipy.sparse import load_npz
from MITgcmutils import mds
from time import time
from project_velocity import project_velocity
import sys 

Nr = 58
Nx = 96
print(Nr*Nx*Nx*6)
testing = False

if testing:
        import matplotlib.pyplot as plt
        import matplotlib.tri as tri
        fig, ( (axu1, axv1, axe1), (axu2, axv2, axe2) ) = plt.subplots(ncols=3,nrows=2, figsize=(15,8))

grid_cs510 = "/home/hay/Research/EuropaOceanDynamics/MITgcm_grids/cs510/"
grid_cs96 = "/home/hay/Research/EuropaOceanDynamics/MITgcm_grids/cs96/"
grid_cs510 = grid_cs96

dir_iterp_matrices = "/home/hay/Research/EuropaOceanDynamics/MITgcm_grids/"
dir_data = "./"

iteration = int(sys.argv[1])

interp_u = load_npz(dir_iterp_matrices+"cs96_to_cs96_uvel.npz")
interp_v = load_npz(dir_iterp_matrices+"cs96_to_cs96_vvel.npz")
interp_eta = load_npz(dir_iterp_matrices+"cs96_to_cs96_eta.npz")

if testing:
    CS = mds.rdmds(grid_cs96 + "AngleCS")
    SN = mds.rdmds(grid_cs96 + "AngleSN")

    X1 = mds.rdmds(grid_cs96 + 'XC').ravel()
    Y1 = mds.rdmds(grid_cs96 + 'YC').ravel()
    X2 = mds.rdmds(grid_cs510 + 'XC').ravel()
    Y2 = mds.rdmds(grid_cs510 + 'YC').ravel()

    triang1 = tri.Triangulation(X1, Y1)
    triang2 = tri.Triangulation(X2, Y2)

data_u_all = mds.rdmds(dir_data + "U", itrs=iteration)
data_v_all = mds.rdmds(dir_data + "V", itrs=iteration)
data_eta = mds.rdmds(dir_data + "Eta", itrs=iteration)

u_start = np.zeros( (Nr, Nx, 6*Nx) )
v_start = np.zeros( (Nr, Nx, 6*Nx) )
eta_start = np.zeros( (Nr, Nx, 6*Nx) )

eta_start = interp_eta.dot( data_eta.ravel() ).reshape( (Nx, 6*Nx) )

if testing:
    axe1.tricontourf(triang1, data_eta.ravel())
    axe2.tricontourf(triang2, eta_start.ravel())

for r in range(Nr):
    
    if testing:
        data_u, data_v = project_velocity(data_u_all[r], data_v_all[r], CS, SN)
    else:
        data_u, data_v = (data_u_all[r], data_v_all[r])

    u_start[r] = interp_u.dot( data_u.ravel() ).reshape( ( Nx, 6*Nx) )
    v_start[r] = interp_v.dot( data_v.ravel() ).reshape( ( Nx, 6*Nx) )

    if testing:
        axu1.tricontourf(triang1, data_u.ravel())
        axu2.tricontourf(triang2, u_start[r].ravel())
  
        axv1.tricontourf(triang1, data_v.ravel())
        axv2.tricontourf(triang2, v_start[r].ravel())

        fig.savefig("/home/hay/Research/SyncFolder/init_plot{:02d}.png".format(r),dpi=500)  
        print("saved... ", r)

    # interp_eta = load_npz(dir_iterp_matrices+"cs96_to_cs96_eta.npz")
    # data = mds.rdmds(dir_data + "Eta", itrs=iteration)

    print(r)

eta_start.astype('>f4').tofile( 'etaInit{:d}.bin'.format(iteration) )
u_start.astype('>f4').tofile(   'uInit{:d}.bin'.format(iteration) )
v_start.astype('>f4').tofile(   'vInit{:d}.bin'.format(iteration) )