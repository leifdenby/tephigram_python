import numpy as np
import matplotlib.pyplot as plot
import math

T1=-10.
T2=40.
theta1=-12.

rotation_angle = 45.
rotation_origin = [T1, theta1]

def rotation_matrix_2d(angle):
    return np.array([[np.cos(angle), -np.sin(angle)],[np.sin(angle), np.cos(angle)]])

from matplotlib.transforms import Affine2D

import mpl_toolkits.axisartist.floating_axes as floating_axes

import numpy as np
import  mpl_toolkits.axisartist.angle_helper as angle_helper
from matplotlib.projections import PolarAxes
from mpl_toolkits.axisartist.grid_finder import FixedLocator, MaxNLocator, \
     DictFormatter

from  mpl_toolkits.axisartist.grid_helper_curvelinear import GridHelperCurveLinear
from mpl_toolkits.axisartist import Subplot

plot.ion()

theta_min, theta_max = -20., 110.
T_min, T_max = -70., 40.

d_theta = 10.
d_T = 10.

def setup_axes1(fig, rect):
    """
    A simple one.
    """
    deg = -45.
    tr = Affine2D().rotate_deg(deg)

    T_ticks = np.arange(T_min, T_max, d_T)
    theta_ticks = np.arange(theta_min, theta_max, d_T)

    grid_helper = GridHelperCurveLinear(tr, grid_locator1=FixedLocator(T_ticks), grid_locator2=FixedLocator(theta_ticks))

    ax1 = Subplot(fig, 1, 1, 1, grid_helper=grid_helper, transform=tr)
    # ax1 will have a ticks and gridlines defined by the given
    # transform (+ transData of the Axes). Note that the transform of
    # the Axes itself (i.e., transData) is not affected by the given
    # transform.

    fig.add_subplot(ax1)

    tf = lambda x,y: tr.transform(np.array([x, y]).T).T

    # SW, SE, NE, NW
    corners = np.array([[-25., -20.], [30., 40.], [-40., 120.], [-105., 60.]])
    corners_t = tf(corners[:,0], corners[:,1])

    ax1.set_aspect(1.)
    ax1.set_xlim(-40., 60.0)
    ax1.set_ylim(-40., 120.)

    ax1.axis["t"]=ax1.new_floating_axis(0, 0.)
    T_axis = ax1.axis['t']
    theta_axis = ax1.axis["t2"]=ax1.new_floating_axis(1, 0.)
    
    plot.draw()
    plot.show()

    return ax1, tf


fig = plot.figure(figsize=(10,20))
import ipdb
with ipdb.launch_ipdb_on_exception():
    ax1, tf = setup_axes1(fig, '111')

# plot temperature lines and label them
for T in np.arange(T_min, T_max+0.1, d_T):
    ax1.plot(*tf([T, T], [theta_min, theta_max]  ), linestyle='-', color='red')

# plot potential temperature lines and label them
for theta in np.arange(theta_min, theta_max, d_theta):
    ax1.plot(*tf([T_min, T_max], [theta, theta]  ), linestyle='-', color='green')

raw_input()

