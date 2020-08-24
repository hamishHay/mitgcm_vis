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
import os
from matplotlib.ticker import FormatStrFormatter


matplotlib.rcParams['axes.unicode_minus'] = False
matplotlib.rcParams['axes.linewidth'] = 0.5
matplotlib.rc('image', cmap='coolwarm')


# initialize figure:
fig = plt.figure(figsize=(8,5))
R = 1550e3
Rb = 1450e3
Rb= 0.9*R




def get_slice_axes(fig, num, R, Rb):

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
  radius_ticks = [(R/Rb, '0'),
                # (1.5, '%i' % (MAX_R/2.)),
                (R, 'R')]

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
                                                  extremes=(-np.pi/2, np.pi/2, Rb/R, 1),
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
  ax1.axis["bottom"].major_ticklabels.set_visible(False)
  ax1.axis["bottom"].set_axis_direction("left")

  ax1.axis['top'].major_ticklabels.set_visible("True")
  ax1.axis['top'].set_axis_direction('right')
# 
#ax1.axis["bottom"].major_ticklabels.locs_angles_labels = [0, 0, 0, 0,0,0,0] 
# print(ax1.axis["bottom"].major_ticklabels._locs_angles_labels)
#ax1.axis["bottom"].major_ticklabels.set_ha('right')

#  ax1.axis['bottom'].set_label_position('top') 

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

ax1s = get_slice_axes(fig, 121, R, Rb)
ax2s = get_slice_axes(fig, 122, R, Rb*0.7)


def zonal_sheet(lat, r):
   x = r*np.cos(np.radians(lat))
   return np.sin(3*2*np.pi*x)

lats = np.linspace(-90,90,200)
rs1 = np.linspace(Rb, R, 81)/R
rs2 = np.linspace(Rb*0.7, R, 81)/R

ll, rr = np.meshgrid(lats, rs1)
ll, rr2 = np.meshgrid(lats, rs2)

U = zonal_sheet(ll, rr)

print(U)

c = ax1s[1].pcolormesh(np.radians(ll), rr, U,linewidth=0, shading='gouraud')
c.set_edgecolor('face')
plt.colorbar(c, ax=ax1s[0], shrink=0.8, pad=0.1)#, orientation='horizontal')
#ax1s[1].contour(ll2, rr2, W[:,:,lon_slice], levels=[0], colors='k', linestyles='--', linewidths=0.8)

c = ax2s[1].pcolormesh(np.radians(ll), rr2, U,linewidth=0, shading='gouraud')
c.set_edgecolor('face')
plt.colorbar(c, ax=ax2s[0], shrink=0.8, pad=0.1)#, orientation='horizontal')


ax2s[0].set_title("Test data stretched\nto thicker ocean\n")
ax1s[0].set_title("Test data set\n")
fig.savefig("projection_test.pdf")



os.system("echo 'plot' | mailx -s 'plot from plot_slice.py' -A projection_test.pdf hamish.hay@jpl.nasa.gov")
