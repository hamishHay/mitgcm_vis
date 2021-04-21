import mpl_toolkits.axisartist.floating_axes as floating_axes
from matplotlib.projections import PolarAxes
from mpl_toolkits.axisartist.grid_finder import FixedLocator, \
     MaxNLocator, DictFormatter
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.transforms import Affine2D
import matplotlib

from matplotlib.patches import ConnectionPatch
from matplotlib.projections import get_projection_class
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

def get_polar_inset(fig, ax, th1, th2, rt, rb, x=0.0, y=0.0, scale="25%"):
    axins = inset_axes(ax, width=scale, height=scale, 
                       axes_class=get_projection_class("polar"), 
                       borderpad=0, 
                       bbox_to_anchor=(x,y,1,1), 
                       bbox_transform=ax2.transAxes)
    
    axins.set_thetamin(th1)
    axins.set_thetamax(th2)
    axins.set_zorder(20)
    axins.set_ylim( [rb, rt] )
    axins.set_rorigin(0)
    
    axins.grid(False,which='both')
    axins.tick_params(labelleft=False, labelright=False,
                labeltop=False, labelbottom=False)

    
    coords = [(th2, rt), (th1, rt), (th2, rb), (th1, rb)]
    for xl, yl in coords:
        xl = np.radians(xl)
        con = ConnectionPatch(xyA=(xl, yl), coordsA='data', axesA=ax,
                        xyB=(xl, yl), coordsB='data', axesB=axins, alpha=0.7)
        con.set_color([0, 0, 0])
        fig.add_artist(con)
        con.set_linewidth(0.5)

    angles = np.radians(np.linspace(th1, th2, 101))
    upper_r = np.ones(101)*rt
    lower_r = np.ones(101)*rb

    angles = np.concatenate( (angles, angles[::-1]) )
    rvals = np.concatenate(  (upper_r, lower_r) )

    ax.fill(angles, rvals, color='C1', alpha = 1, facecolor="none",zorder=1e8, lw=1.2)

    return axins



def get_slice_axes(fig, num, R, Rb, gs=None, lat1=-90.0, lat2=90.0, shrink=1.0):

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
  radius_ticks = [(Rb, ''),
                # (1.5, '%i' % (MAX_R/2.)),
                (R, '')]

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
  l1 = np.radians(lat1)
  l2 = np.radians(lat2)
  grid_helper = floating_axes.GridHelperCurveLinear(tr,
                                                  extremes=(l1, l2, R, Rb),
                                                  grid_locator1=grid_locator1,
                                                  grid_locator2=grid_locator2,
                                                  tick_formatter1=tick_formatter1,
                                                  tick_formatter2=tick_formatter2)


  if type(num) is tuple:
    ax1 = floating_axes.FloatingSubplot(fig, num[0], num[1], num[2], grid_helper=grid_helper)
  else:
    if gs is not None:
      ax1 = floating_axes.FloatingSubplot(fig, gs, grid_helper=grid_helper)
    else:
      ax1 = floating_axes.FloatingSubplot(fig, num, grid_helper=grid_helper)


  fig.add_subplot(ax1)


  ax1.axis["left"].set_axis_direction("bottom")
  ax1.axis["right"].set_axis_direction("bottom")
  ax1.axis["left"].major_ticklabels.set_rotation(90)
  ax1.axis["bottom"].set_axis_direction("right")
  ax1.axis["bottom"].major_ticks.set_tick_out(True)
#  ax1.axis["bottom"].major_ticks.set_ticks([])


  majortick_iter, minortick_iter = ax1.axis['bottom']._axis_artist_helper.get_tick_iterators(ax1)
  tick_loc_angle, ticklabel_loc_angle_label = ax1.axis['bottom']._get_tick_info(majortick_iter)

  # ax1.axis["bottom"].toggle(all=False)#, ticks=True)
  ax1.axis["top"].toggle(all=False)#, ticks=True)
  ax1.axis["left"].toggle(all=False)#, ticks=True)
  ax1.axis["right"].toggle(all=False)#, ticks=True)

  aux_ax = ax1.get_aux_axes(tr)

  aux_ax.patch = ax1.patch # for aux_ax to have a clip path as in ax
  ax1.patch.zorder=0.9 # but this has a side effect that the patch is
                    #  drawn twice, and possibly over some other
                    #  artists. So, we decrease the zorder a bit to
                    #  prevent this.

  return ax1, aux_ax

