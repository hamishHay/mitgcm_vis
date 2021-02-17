import numpy as np 
import matplotlib.pyplot as plt 
import argparse
import os
import h5py
import cartopy.crs as ccrs
from plotting_routines import get_slice_axes
import time

plt.rcParams.update({'axes.titlesize': 10})


parser = argparse.ArgumentParser()

parser.add_argument("file", help="File name to read data from.")
parser.add_argument("-f", "--field", action='append', help="Field(s) to read from the data field. Multiple fields can be summed with '+'.")
args = parser.parse_args()

num_fields = len(args.field)
fields = args.field

# Plotting rules!
# Each field will be given one row.
# Each row will be split into four columns:
# (1) Section 
# (2) Same as (1) but not projected
# (3) Projected map view at the surface 
# (4) Projected map view at the middle of the ocean

data_file = h5py.File(args.file, 'r')

fig = plt.figure(figsize=(12, 6*num_fields))
gs = fig.add_gridspec(ncols=3, nrows=2*num_fields, width_ratios=[3.5,4,4])

R = 1550e3
Rb = R - 100e3
for i in range(num_fields):

    if "+" in fields[i]:
        multi_fields = fields[i].split('+')
        data = np.zeros( data_file[multi_fields[0]].shape, dtype=np.float32 )
        for j in range(len(multi_fields)):
            data += data_file[ multi_fields[j] ]
    else:
        data = data_file[ fields[i] ]
    
    x, y = data_file['lon'], data_file['lat']
    RC = data_file['RC'][:].flatten() 

    ax1s = get_slice_axes(fig, 1, R, Rb, gs=gs[2*i:2*i+2,0])
    ax1s[0].set_aspect('equal')

    r_max = np.amax(RC)
    r_min = np.amin(RC)

    RC_norm = (R/Rb - 1.0) * ((RC - r_min) / (r_max - r_min)) + 1.0
    ll2, rr2 = np.meshgrid(np.radians(y), RC_norm)
    ax1s[1].pcolormesh(ll2, rr2, np.mean(data,axis=-1) )

    # (3) Projected map view at the surface 
    ax2 = fig.add_subplot(gs[2*i,1], projection="mollweide")#, projection=ccrs.Mollweide())
    ax2.set_title(fields[i] + ": surface")


    ax2.pcolormesh(np.radians(x), np.radians(y), data[1])#, transform=ccrs.PlateCarree(), fast=True)
    ax2.set_yticklabels([])
    ax2.set_xticklabels([])

    ax3 = fig.add_subplot(gs[2*i,2], projection="mollweide")
    ax3.set_title(fields[i] + ": mid-depth")

    ax3.pcolormesh(np.radians(x), np.radians(y), data[29])#, transform=ccrs.PlateCarree(), fast=True)
    ax3.set_yticklabels([])
    ax3.set_xticklabels([])


    ax4 = fig.add_subplot(gs[2*i+1,1])
    ax4.set_title(fields[i] + ": zonal avg")
    ax4.set_xlabel("Latitude")
    ax4.set_ylabel("Depth")


    ax4.contourf(y, RC/1e3, np.mean(data,axis=-1), levels=11)

    # (4) Projected map view at depth 
    ax5 = fig.add_subplot(gs[2*i+1,2], projection="mollweide")
    ax5.set_title(fields[i] + ": bottom")

    ax5.pcolormesh(np.radians(x), np.radians(y), data[-1])#, transform=ccrs.PlateCarree(), fast=True)
    ax5.set_yticklabels([])
    ax5.set_xticklabels([])




fig.savefig("field_plot.png",dpi=500, bbox_inches='tight')
os.system("echo 'plot' | mailx -s 'plot from plot_field.py' -A field_plot.png hamish.hay@jpl.nasa.gov")