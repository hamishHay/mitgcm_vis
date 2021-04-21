import numpy as np 
import matplotlib.pyplot as plt 
import argparse
import os
import h5py
import cartopy.crs as ccrs
from plotting_routines import get_slice_axes
import time
import seaborn as sns
# plt.rcParams['image.cmap']=sns.color_palette("vlag", as_cmap=True)

sns.set_palette("vlag")

cmap = sns.color_palette("vlag", as_cmap=True)

plt.rcParams.update({'axes.titlesize': 10})


parser = argparse.ArgumentParser()

parser.add_argument("file", help="File name to read data from.")
parser.add_argument("-f", "--field", action='append', help="Field(s) to read from the data field. Multiple fields can be summed with '+'.")
parser.add_argument("-p", "--proj", default='None', help='Projection type to use for map plots')
args = parser.parse_args()

print(args.file)
print(args.field)

if args.proj is 'None':
    proj = None
else:
    proj = args.proj

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

fig = plt.figure(figsize=(5*num_fields, 15))
gs = fig.add_gridspec(ncols=num_fields, nrows=3)#, width_ratios=[4,4,4])

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

    # ax1s = get_slice_axes(fig, 1, R, Rb, gs=gs[2*i:2*i+2,0])
    # ax1s[0].set_aspect('equal')

    x = x[:]
    y = y[:]
    
    xl = -180#-10
    xh = 180#10
    yl = -90#-10
    yh = 90#10
    rl = 0

    xb = (np.abs(x - xl)).argmin()
    xt = (np.abs(x - xh)).argmin()
    yb = (np.abs(y - yl)).argmin()
    yt = (np.abs(y - yh)).argmin()
    rl = (np.abs(RC - rl)).argmin()

    y_eq = (np.abs(y - 0.0)).argmin()

    data = data[:, yb:yt+1, xb:xt+1]
    x = x[xb:xt+1]
    y = y[yb:yt+1]

    cval = max( abs(data.min()), abs(data.max()) )

    # r_max = np.amax(RC)
    # r_min = np.amin(RC)

    # RC_norm = (R/Rb - 1.0) * ((RC - r_min) / (r_max - r_min)) + 1.0
    # ll2, rr2 = np.meshgrid(np.radians(y), RC_norm)
    # ax1s[1].pcolormesh(ll2, rr2, np.mean(data,axis=-1) )


    # map at const R level 
    ax1 = fig.add_subplot(gs[0,i], projection=ccrs.Mollweide())
    ax1.set_title(fields[i] + ": depth = " + str(RC[rl]/1e3) + 'km')
    ax1.pcolormesh(x, y, data[rl], cmap=cmap, vmin=-cval, vmax=cval, transform=ccrs.PlateCarree(), shading="flat")#, fast=True)
    ax1.set_aspect('equal')
    
    
    # map at const lon 
    ax2 = fig.add_subplot(gs[1,i])
    ax2.set_title(fields[i] + ": const lat")
    ax2.set_xlabel("Longitude")
    ax2.set_ylabel("Depth [km]")

    mean_val = abs(np.mean(data))
    c2 = ax2.pcolormesh(x, RC/1e3, data[:,y_eq,:], cmap=cmap, vmin=-cval, vmax=cval)#, vmin=-1, vmax=1)

    print(data[5,y_eq,:])

    # ax2.set_xlim([40,50])
    # ax2.set_ylim([-30,0])
    # ax2.set_aspect(0.3)

    

    # map at const lat
    ax3 = fig.add_subplot(gs[2,i])
    ax3.set_title(fields[i] + ": const lon")
    ax3.set_xlabel("Latitude")
    ax3.set_ylabel("Depth [km]")

    mean_val = abs(np.mean(data))
    c3 = ax3.pcolormesh(y, RC/1e3, data[:,:,0], cmap=cmap, vmin=-cval, vmax=cval)

    # ax3.set_xlim([40,50])
    # ax3.set_ylim([-30,0])
    # ax3.set_aspect(0.3)

    plt.colorbar(c3, orientation='horizontal', label='Velocity [m/s]')

    # (4) Projected map view at depth 
    # ax5 = fig.add_subplot(gs[2*i+1,2], projection=proj)
    # ax5.set_title(fields[i] + ": bottom")

    # ax5.pcolormesh(np.radians(x), np.radians(y), data[-1], cmap=cmap)#, transform=ccrs.PlateCarree(), fast=True)
    # ax5.set_yticklabels([])
    # ax5.set_xticklabels([])


# import pickle
# with open('myplot.pkl','wb') as fid:
#     pickle.dump(fig, fid)
# os.system("echo 'plot' | mailx -s 'plot from plot_field.py' -A myplot.pkl hamish.hay@jpl.nasa.gov")
fig.savefig("/home/hay/Research/SyncFolder/field_plot.png",dpi=500)#, bbox_inches='tight')
# os.system("echo 'plot' | mailx -s 'plot from plot_field.py' -A field_plot.png hamish.hay@jpl.nasa.gov")