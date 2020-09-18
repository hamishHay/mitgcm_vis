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

iter = int(sys.argv[1])

var_list = ['U', 'V', 'XC', 'YC', 'RC', 'W', 'CS', 'SN', 'Eta', 'XG', 'YG']

ds = open_mdsdataset('.', prefix=var_list, iters=iter, geometry='curvilinear', ignore_unknown_vars=True, default_dtype=np.float)

XC = np.array(ds.XC)
YC = np.array(ds.YC)
XG = np.array(ds.XG)
YG = np.array(ds.YG)

Eta = np.array(ds.Eta)
RC = -mds.rdmds('RC').flatten()
#print(RC)

triang = tri.Triangulation(XC.flatten(), YC.flatten())

U = np.array(ds.U)[0,:,:,:]*100
V = np.array(ds.V)[0,:,:,:]*100
W = np.array(ds.W)[0,:,:,:]*100
Eta = Eta[0,:,:]
Us = U[0,:,:]
Vs = V[0,:,:]
print(Eta.shape, Us.shape)


CS = np.array(ds.CS)
SN = np.array(ds.SN)

hres = 1.0                                                                              
lon = np.arange(-180, 180.01, hres) #+ hres/2 
lat = np.arange(-90, 90.01, hres) #+ hres/2
lon_out, lat_out = np.meshgrid(lon, lat)

U, V = project_velocity(U, V, CS, SN)

Us, Vs = project_velocity(Us, Vs, CS, SN)

U = interpolate_data(U, XC, YC, lon_out, lat_out)
W = interpolate_data(W, XC, YC, lon_out, lat_out)
V = interpolate_data(V, XC, YC, lon_out, lat_out)
#Eta = interpolate_data(Eta, XC, YC, lon_out, lat_out)
r_slice = 0
lon_slice = 90
lat_slice = 90
#Us = U[r_slice, :, :]
#Vs = V[r_slice, :, :]
#Umap = interpolate_data(Umap, XC, YC, lon_out, lat_out)

print(U.shape)


#U = U[:, :, 0]
#W = W[:, :, 0]
#V = V[:, :, 0]


Umin = min( np.amin(Us), np.amin(U) )
Umax = max( np.amax(Us), np.amax(U) )
Vmin = min( np.amin(Vs), np.amin(V) )
Vmax = max( np.amax(Vs), np.amax(V) )
#Umin = min( np.amin(Us), np.amin(U) )
#Umax = max( np.amax(Us), np,amax(U) )
#Umin = min( np.amin(Us), np.amin(U) )
#Umax = max( np.amax(Us), np,amax(U) )

MAX_R = 40000.
radius = np.random.rand(100)*MAX_R
radius = radius/np.max(radius) + 1.

# initialize figure:
fig = plt.figure(figsize=(14,14))
R = 1550e3
Rb = 1450e3
Rb= R - 100e3

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
  radius_ticks = [(1., ''),
                # (1.5, '%i' % (MAX_R/2.)),
                (R/Rb, '')]

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
  #ax1.axis["right"].set_tick_direction('out')
  #ax1.axis
# ax1.axis["left"].major_ticks.set_tick_out(True)
# # ax1.axis["bottom"].set_visible(False)
  ax1.axis["bottom"].set_axis_direction("right")
  ax1.axis["bottom"].major_ticks.set_tick_out(True)
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

lats = np.radians(np.linspace(-90, 90, len(U[0])))
r = np.linspace(0, 1, len(U)) + 1

ll, rr = np.meshgrid(lats, r)
data = 1.0*rr**2.0

# create a parasite axes whose transData in RA, cz

# plot your data:
# aux_ax.plot(theta, radius)

r_max = np.amax(RC)
r_min = np.amin(RC)

RC = (R/Rb - 1.0) * ((RC - r_min) / (r_max - r_min)) + 1.0

W = W[::-1, :]
U = U[::-1, :]
V = V[::-1, :]

print(Us[:,0], U[-1,:])

ll2, rr2 = np.meshgrid(np.radians(lat), RC)
llon2, rr3 = np.meshgrid(np.radians(lon), RC)
print(R/Rb, RC)

ax1s = get_slice_axes(fig, 431)
ax2s = get_slice_axes(fig, 432)
ax3s = get_slice_axes(fig, 433)
ax4s = get_slice_axes(fig, 434)
ax5s = get_slice_axes(fig, 435)
ax6s = get_slice_axes(fig, 436)



Umin, Umax = [min(Umin, -Umax), max(-Umin, Umax)]

level_u = np.linspace(Umin, Umax, 11)
level_v = np.linspace(Vmin, Vmax, 11)

c = ax1s[1].contourf(ll2, rr2, W[:,:,lon_slice], 11)
plt.colorbar(c, ax=ax1s[0], shrink=0.8, pad=0.1)#, orientation='horizontal')
#ax1s[1].contour(ll2, rr2, W[:,:,lon_slice], levels=[0], colors='k', linestyles='--', linewidths=0.8)
ax1s[0].set_ylabel("Longitude " + str(lon[lon_slice]) )
ax1s[0].yaxis.set_label_coords(10,10.02)

c = ax2s[1].contourf(ll2, rr2, U[:,:,lon_slice], levels=level_u)

#W = W*0 + 0.1
sk=6
skl=2
vproj_x, vproj_y = [-V[::sk,::skl,lon_slice]*np.sin(ll2[::sk,::skl]), -V[::sk,::skl,lon_slice]*-np.cos(ll2[::sk,::skl]) ]
wproj_x, wproj_y = [W[::sk,::skl,lon_slice]*np.cos(ll2[::sk,::skl]), W[::sk,::skl,lon_slice]*np.sin(ll2[::sk,::skl]) ]

#ax2s[1].quiver(ll2[::sk,::skl],rr2[::sk,::skl],
#               vproj_x + wproj_x, 
#               vproj_y + wproj_y , zorder=10, pivot='mid')

#c = ax2s[1].contourf(llon2, rr3, U)
#ax2s[1].contour(ll2, rr2, U[:,:,lon_slice], levels=[0], colors='k', linestyles='--', linewidths=0.8)
plt.colorbar(c, ax=ax2s[0], shrink=0.8, pad=0.1)#, orientation='horizontal')

c = ax3s[1].contourf(ll2, rr2, V[:,:,lon_slice], levels=level_v)
#ax3s[1].contour(ll2, rr2, V[:,:,lon_slice], levels=[0], colors='k', linestyles='--', linewidths=0.8)
plt.colorbar(c, ax=ax3s[0], shrink=0.8, pad=0.1)#, orientation='horizontal')

c = ax4s[1].contourf(llon2, rr3, W[:,lat_slice,:], 11)
plt.colorbar(c, ax=ax4s[0], shrink=0.8, pad=0.1)#, orientation='horizontal')
#ax1s[1].contour(llon2, rr3, W[:,lat_slice,:], levels=[0], colors='k', linestyles='--', linewidths=0.8)
ax4s[0].set_ylabel("Latitude " + str(lat[lat_slice]) )

c = ax5s[1].contourf(llon2, rr3, U[:,lat_slice,:], levels=level_u)
#c = ax2s[1].contourf(llon2, rr3, U)
#ax2s[1].contour(llon2, rr3, U[:,lat_slice,:], levels=[0], colors='k', linestyles='--', linewidths=0.8)
plt.colorbar(c, ax=ax5s[0], shrink=0.8, pad=0.1)#, orientation='horizontal')
#ax5s[1].quiver(llon2[::5, :], rr3[::5, :], W[::5, lat_slice,:], U[::5, lat_slice, :])

c = ax6s[1].contourf(llon2, rr3, V[:,lat_slice,:], levels=level_v)
#ax3s[1].contour(llon2, rr3, V[:,lat_slice,:], levels=[0], colors='k', linestyles='--', linewidths=0.8)
plt.colorbar(c, ax=ax6s[0], shrink=0.8, pad=0.1)#, orientation='horizontal')

ax4 = fig.add_subplot(437)
ax5 = fig.add_subplot(438)
ax6 = fig.add_subplot(439)

ax7 = fig.add_subplot(4,3,10)
ax8 = fig.add_subplot(4,3,11)
ax9 = fig.add_subplot(4,3,12)

ax4.set_aspect('equal')
ax5.set_aspect('equal')
ax6.set_aspect('equal')

ax8.set_aspect('equal')
ax9.set_aspect('equal')

#filter = True
#if filter:
from sh_filter import sh_filter
Usf = sh_filter(Us, XG, YG, lmax=8)
Vsf = sh_filter(Vs, XG, YG, lmax=8)*-1
Etaf = sh_filter(Eta, XC, YC, lmax=8)
lonf = np.linspace(-180, 180, len(Etaf[0]))
latf = np.linspace(-90, 90, len(Etaf))

print(Etaf.shape)

c = ax4.tricontourf(triang, Eta.flatten(), cmap=plt.cm.viridis)
c4 = plt.colorbar(c, ax=ax4, shrink=0.8, orientation='horizontal')
#ax4.streamplot(lon_out, lat_out, Us, Vs)

c = ax7.contourf(lonf, latf, Etaf, cmap=plt.cm.viridis)
c4 = plt.colorbar(c, ax=ax7, shrink=0.8, orientation='horizontal')


#if not filter:
c = ax5.tricontourf(triang, Us.flatten(), levels=level_u)
c5 = plt.colorbar(c, ax=ax5, shrink=0.8, orientation='horizontal')

#else:
c = ax8.contourf(lonf, latf, Usf)
c5 = plt.colorbar(c, ax=ax8, shrink=0.8, orientation='horizontal')


#if not filter:
c = ax6.tricontourf(triang, Vs.flatten(), levels=level_v)
c6 = plt.colorbar(c, ax=ax6, shrink=0.8, orientation='horizontal')

#else:
c = ax9.contourf(lonf, latf, Vsf)
c6 = plt.colorbar(c, ax=ax9, shrink=0.8, orientation='horizontal')


skip = 5
#ax7.quiver(lonf[::skip], latf[::skip], Usf[::skip,::skip], Vsf[::skip,::skip], color='w')
#ax4.quiver(XC, YC, Us, Vs, color='w',scale=2000)

#ax4.set_xlim([-100, -75])
#ax4.set_ylim([15,45])

ax1s[0].set_title("Radial vel [cm/s]\n")
ax2s[0].set_title("Zonal vel [cm/s]\n")
ax3s[0].set_title("Merid vel [cm/s] \n")
ax4.set_title("Displacement [m]")
ax5.set_title("Zonal vel [cm/s]")
ax6.set_title("Merid vel [cm/s]")

ax7.set_title("Displacement [m]")
ax8.set_title("Zonal vel [cm/s]")
ax9.set_title("Merid vel [cm/s]")
#aux_ax.streamplot(ll2, rr, V, W, linewidth=0.9, density=0.5)

ax7.set_ylabel("Low-pass filter ($l\leq8$)")


c5.ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))


#fig, (ax1, ax2) = plt.subplots(ncols=2)

#c = ax1.contourf(lats, RC/1e3, U)
#c2 = ax2.contourf(Umap)
#plt.colorbar(c, ax=ax1, orientation='horizontal')
#plt.colorbar(c2,ax=ax2, orientation='horizontal')
fig.savefig("slice_test.png",dpi=300)



os.system("echo 'plot' | mailx -s 'plot from plot_slice.py' -A slice_test.png hamish.hay@jpl.nasa.gov")
