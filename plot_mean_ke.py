import mpl_toolkits.axisartist.floating_axes as floating_axes
from matplotlib.projections import PolarAxes
from mpl_toolkits.axisartist.grid_finder import FixedLocator, \
     MaxNLocator, DictFormatter
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.transforms import Affine2D
import matplotlib
from MITgcmutils import mds
from xmitgcm import open_mdsdataset
import pyresample
import sys 
sys.path.append("/home/hay/Research/python/plotting/mitgcm_vis")
from project_velocity import project_velocity, interpolate_data
import os

matplotlib.rcParams['axes.unicode_minus'] = False
matplotlib.rcParams['axes.linewidth'] = 0.8

times = np.arange(100, 3001, 100, dtype=np.int)


var_list = ['U', 'V', 'W']
ds2 = open_mdsdataset('.', prefix=['CN', 'SN', 'XC', 'YC'], iters=100, geometry='curvilinear', ignore_unknown_vars=True)
CS = np.array(ds2.CS)
SN = np.array(ds2.SN)
XC = np.array(ds2.XC)
YC = np.array(ds2.YC)

hres = 4.0                                                                              
lon = np.arange(-180, 180.01, hres) #+ hres/2 
lat = np.arange(-90, 90.01, hres) #+ hres/2
lon_out, lat_out = np.meshgrid(lon, lat)


EK = np.zeros(len(times))

for i in range(len(times)):
  time = times[i]
  print("TIME: ", time/times[-1])
  ds = open_mdsdataset('.', prefix=var_list, iters=time, geometry='curvilinear', ignore_unknown_vars=True, default_dtype=np.float)

  U = np.array(ds.U)[0,:,:]
  V = np.array(ds.V)[0,:,:]
  W = np.array(ds.W)[0,:,:]

  #U, V = project_velocity(U, V, CS, SN)

  #U = interpolate_data(U, XC, YC, lon_out, lat_out)
  #W = interpolate_data(W, XC, YC, lon_out, lat_out)
  #V = interpolate_data(V, XC, YC, lon_out, lat_out)

  U2 = np.sqrt(U*U + V*V + W*W)

  EK[i] = np.sum(U2)

#EK = np.sum(EK, axis=(1, 2, 3))

fig, ax = plt.subplots()

ax.plot(EK)
ax.set_ylabel("Total Kinetic Energy")
ax.set_xlabel("Time [Orbits]")

fig.savefig("avg_EK.pdf", bbox_inches='tight')


os.system("echo 'plot' | mailx -s 'plot from plot_mean_ke.py' -A avg_EK.pdf hamish.hay@jpl.nasa.gov")
