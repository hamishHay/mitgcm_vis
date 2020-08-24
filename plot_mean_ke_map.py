import mpl_toolkits.axisartist.floating_axes as floating_axes
from matplotlib.projections import PolarAxes
from mpl_toolkits.axisartist.grid_finder import FixedLocator, \
     MaxNLocator, DictFormatter
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.transforms import Affine2D
import matplotlib
import matplotlib.tri as tri
from MITgcmutils import mds
from xmitgcm import open_mdsdataset
import pyresample
import sys 
sys.path.append("/home/hay/Research/python/plotting/mitgcm_vis")
from project_velocity import project_velocity, interpolate_data
from calc_tavg import calc_mean_EK
import os
from matplotlib.ticker import FormatStrFormatter


matplotlib.rcParams['axes.unicode_minus'] = False
matplotlib.rcParams['axes.linewidth'] = 0.5
matplotlib.rc('image', cmap='coolwarm')

iter = int(sys.argv[1])

var_list = ['XC', 'YC', 'RC']

ds = open_mdsdataset('.', prefix=var_list, iters=iter, geometry='curvilinear', ignore_unknown_vars=True, default_dtype=np.float)

XC = np.array(ds.XC)
YC = np.array(ds.YC)
RC = -mds.rdmds('RC').flatten()

triang = tri.Triangulation(XC.flatten(), YC.flatten())

hres = 0.1                                                                              
lon = np.arange(-180, 180.01, hres) #+ hres/2 
lat = np.arange(-90, 90.01, hres) #+ hres/2
lon_out, lat_out = np.meshgrid(lon, lat)

EK = calc_mean_EK(start=16000, end=17000)

EK = interpolate_data(EK, XC, YC, lon_out, lat_out)

r_slice = 0
lon_slice = 90
lat_slice = 90




MAX_R = 40000.
radius = np.random.rand(100)*MAX_R
radius = radius/np.max(radius) + 1.

# initialize figure:
fig = plt.figure(figsize=(12,4))
R = 1550e3
Rb = 1450e3
Rb= 0.8*R

def get_slice_axes(fig, num):

  tr_rotate = Affine2D().translate(0, 90)
  # set up polar axis
  tr = PolarAxes.PolarTransform() #+ tr_rotate

  angle_ticks = [(np.radians(i), str(i)) for i in np.arange(80, -81, -10) ]
# angle_ticks = [(np.radians(80), r"$80$"), 
#                (np.radians(40), r"$\frac{3}{4}\pi$"),
#                (1.*np.pi, r"$\pi$"),
#                (1.25*np.pi, r"$\frac{5}{4}\pi$"),
#                (1.5*np.pi, r"$\frac{3}{2}\pi$"),
#                (1.75*np.pi, r"$\frac{7}{4}\pi$")]

# set up ticks and spacing around the circle
  grid_locator1 = FixedLocator([v for v, s in angle_ticks])
  tick_formatter1 = DictFormatter(dict(angle_ticks))


# set up grid spacing along the 'radius'
  radius_ticks = [(1., '0'),
                # (1.5, '%i' % (MAX_R/2.)),
                (R/Rb, 'R')]

  grid_locator2 = FixedLocator([v for v, s in radius_ticks])
  tick_formatter2 = DictFormatter(dict(radius_ticks))

# define angle ticks around the circumference:
# set up axis:
# tr: the polar axis setup
# extremes: theta max, theta min, r max, r min
# the grid for the theta axis
# the grid for the r axis
# the tick formatting for the theta axis
# the tick formatting for the r axis
  grid_helper = floating_axes.GridHelperCurveLinear(tr,
                                                  extremes=(-np.pi/2, np.pi/2, R/Rb, 1),
                                                  grid_locator1=grid_locator1,
                                                  grid_locator2=grid_locator2,
                                                  tick_formatter1=tick_formatter1,
                                                  tick_formatter2=tick_formatter2)

  ax1 = floating_axes.FloatingSubplot(fig, num, grid_helper=grid_helper)
  fig.add_subplot(ax1)


  ax1.axis["left"].set_axis_direction("bottom")
  ax1.axis["right"].set_axis_direction("bottom")
  ax1.axis["left"].major_ticklabels.set_rotation(90)
# ax1.axis["right"].set_tick_direction('out')
# ax1.axis["left"].major_ticks.set_tick_out(True)
# # ax1.axis["bottom"].set_visible(False)
  ax1.axis["bottom"].set_axis_direction("right")
# ax1.axis["bottom"].major_ticklabels.locs_angles_labels = [0, 0, 0, 0,0,0,0] 
# print(ax1.axis["bottom"].major_ticklabels._locs_angles_labels)
#ax1.axis["bottom"].major_ticklabels.set_ha('right')

  majortick_iter, minortick_iter = ax1.axis['bottom']._axis_artist_helper.get_tick_iterators(ax1)
  tick_loc_angle, ticklabel_loc_angle_label = ax1.axis['bottom']._get_tick_info(majortick_iter)
#ax1.axis["bottom"].major_ticklabels.set_pad(1.0)
  aux_ax = ax1.get_aux_axes(tr)

  aux_ax.patch = ax1.patch # for aux_ax to have a clip path as in ax
  ax1.patch.zorder=0.9 # but this has a side effect that the patch is
                    #  drawn twice, and possibly over some other
                    #  artists. So, we decrease the zorder a bit to
                    #  prevent this.

  return ax1, aux_ax

lats = np.radians(np.linspace(-90, 90, len(EK[0])))
r = np.linspace(0, 1, len(EK)) + 1

ll, rr = np.meshgrid(lats, r)
data = 1.0*rr**2.0

# create a parasite axes whose transData in RA, cz

# plot your data:
# aux_ax.plot(theta, radius)

r_max = np.amax(RC)
r_min = np.amin(RC)

RC = (R/Rb - 1.0) * ((RC - r_min) / (r_max - r_min)) + 1.0

EK = EK[::-1, :]

ll2, rr2 = np.meshgrid(np.radians(lat), RC)
llon2, rr3 = np.meshgrid(np.radians(lon), RC)
print(R/Rb, RC)

ax1s = get_slice_axes(fig, 131)
ax2s = get_slice_axes(fig, 132)
ax3s = get_slice_axes(fig, 133)


c = ax1s[1].contourf(ll2, rr2, np.log10(np.mean(EK, axis=-1) ), 11, cmap=plt.cm.hot)
plt.colorbar(c, ax=ax1s[0], shrink=0.8, pad=0.1)#, orientation='horizontal')
#ax1s[1].contour(ll2, rr2, W[:,:,lon_slice], levels=[0], colors='k', linestyles='--', linewidths=0.8)
ax1s[0].set_ylabel("Longitude " + str(lon[lon_slice]) )
ax1s[0].yaxis.set_label_coords(10,10.02)

#c = ax2s[1].contourf(ll2, rr2, U[:,:,lon_slice], levels=level_u)

#plt.colorbar(c, ax=ax2s[0], shrink=0.8, pad=0.1)#, orientation='horizontal')

#c = ax3s[1].contourf(ll2, rr2, V[:,:,lon_slice], levels=level_v)
#ax3s[1].contour(ll2, rr2, V[:,:,lon_slice], levels=[0], colors='k', linestyles='--', linewidths=0.8)
#plt.colorbar(c, ax=ax3s[0], shrink=0.8, pad=0.1)#, orientation='horizontal')

#from sh_filter import sh_filter
#Usf = sh_filter(Us, XC, YC, lmax=12)

ax1s[0].set_title("Radial vel [cm/s]\n")
ax2s[0].set_title("Zonal vel [cm/s]\n")
ax3s[0].set_title("Merid vel [cm/s] \n")


fig.savefig("EK_map.pdf")



os.system("echo 'plot' | mailx -s 'plot from plot_mea_ke_map.py' -A EK_map.pdf hamish.hay@jpl.nasa.gov")
