import numpy as np
import matplotlib.pyplot as plot
import math

from matplotlib.transforms import Affine2D

import mpl_toolkits.axisartist.floating_axes as floating_axes

import numpy as np
import  mpl_toolkits.axisartist.angle_helper as angle_helper
from matplotlib.projections import PolarAxes
from mpl_toolkits.axisartist.grid_finder import FixedLocator, MaxNLocator, \
     DictFormatter

from  mpl_toolkits.axisartist.grid_helper_curvelinear import GridHelperCurveLinear
from mpl_toolkits.axisartist import Subplot



def esat(T):
    """
    Based on http://www.iac.ethz.ch/staff/dominik/idltools/atmos_phys/esat.pro
    """

    e1=1013.250
    TK=273.15
    T0=0.
    esat= e1*10**(
            10.79586*(1-TK/T)
            -5.02808*np.log10(T/TK)
            +1.50474*1e-4*(1-10**(-8.29692*(T/TK-1)))
            +0.42873*1e-3*(10**(4.76955*(1-TK/T))-1)
            -2.2195983
            ) 

    return esat

plot.ion()

T1=-10.
T2=40.
theta1=-12.

rotation_angle = 45.
rotation_origin = [T1, theta1]

theta_min, theta_max, d_theta = -20., 110., 10.
T_min, T_max,d_T = -100., 40., 10.
P_min, P_max, d_P = 200., 1000., 100.
qs_ = 1.e-3*np.array([30., 20., 15., 10., 7., 5., 3., 2., 15., 1.0, 0.7,])

x_min = -40.
x_max = 60.0
        
        
T_ = np.arange(T_min, T_max+0.1, d_T)

class Tephigram:
    def __init__(self):
        fig = plot.figure(figsize=(10,20))
        self.setup_axes1(fig, '111')
        self.plot_temp_lines()
        self.plot_pot_temp_lines()
        self.plot_pressure_lines()
        self.plot_qs_lines()
        self.plot_sat_adiabats()

    def plot_sounding(self, P, T, r):
        """
        Input sounding defined by P, T and r
        """
        theta = self.f_theta(P=(P/100.), T=T)

        self.ax1.plot(*self._tf(T-273.15, theta-273.15), marker='o', color='pink')

    def setup_axes1(self, fig, rect):
        """
        A simple one.
        """
        deg = -45.
        self.tr = Affine2D().rotate_deg(deg)

        T_ticks = np.arange(T_min, T_max, d_T)
        theta_ticks = np.arange(theta_min, theta_max, d_T)

        grid_helper = GridHelperCurveLinear(self.tr, grid_locator1=FixedLocator(T_ticks), grid_locator2=FixedLocator(theta_ticks))

        ax1 = Subplot(fig, 1, 1, 1, grid_helper=grid_helper, transform=self.tr)
        # ax1 will have a ticks and gridlines defined by the given
        # transform (+ transData of the Axes). Note that the transform of
        # the Axes itself (i.e., transData) is not affected by the given
        # transform.

        fig.add_subplot(ax1)


        # SW, SE, NE, NW
        corners = np.array([[-25., -20.], [30., 40.], [-40., 120.], [-105., 60.]])
        corners_t = self._tf(corners[:,0], corners[:,1])

        ax1.set_aspect(1.)
        ax1.set_xlim(x_min, x_max)
        ax1.set_ylim(-10., 100.)

        ax1.axis["t"]=ax1.new_floating_axis(0, 0.)
        T_axis = ax1.axis['t']
        theta_axis = ax1.axis["t2"]=ax1.new_floating_axis(1, 0.)
        
        plot.draw()
        plot.show()
        self.ax1 = ax1

    def _tf(self, T, theta):
        """
        Rotate temp and potential temp into the cartesian basis.
        """
        return self.tr.transform(np.array([T, theta]).T).T

    def plot_temp_lines(self):
        for T in T_:
            self.ax1.plot(*self._tf([T, T], [theta_min, theta_max]  ), linestyle='-', color='red')

    def plot_pot_temp_lines(self):
        theta_ = np.arange(theta_min, theta_max, d_theta)
        for theta in theta_:
            self.ax1.plot(*self._tf([T_min, T_max], [theta, theta]  ), linestyle='-', color='green')

    def plot_pressure_lines(self):
        """
        theta = T*(1000/p)^0.286
        """
        for P in np.arange(P_min, P_max+0.1, d_P):
            T__ = np.linspace(T_min, T_max, 1000.)
            theta_constP = -273.15 + self.f_theta(P, T__+273.15)

            x, y = self._tf(T__, theta_constP)
            mask = x>(x_min + 0.1*(x_max-x_min))
            self.ax1.plot(x[mask], y[mask], linestyle='-', color='blue')

            k = np.argmin((x>x_min)*x)

            x0 = x[k]
            y0 = y[k]

            plot.text(x0, y0, '%dhPa' % P, color='blue')

    def f_theta(self, P, T):
        return T*(P/P_max)**-0.286

    def plot_qs_lines(self):
        """
        Make the saturated specific humidity curves
          q_s = 0.622*e_s/p
        """

        for qs in qs_:
            Tk = T_ + 273.15
            Pqs = 0.622*esat(Tk)/qs

            theta = -273.15 + self.f_theta(P=Pqs,T=Tk)

            self.ax1.plot(*self._tf(T_, theta), linestyle=':', color='grey')

    def plot_sat_adiabats(self):
        """
        Use dtheta = -L*theta/(Cp*T) dqs and integrate from T=-40 C
        """
        def f_dtheta_dq(theta, T):
            L = 2257*1.e3 # J/kg
            Cp = 1.004*1.e3 # J/kg/K
            return -L*theta/(Cp*T)

        for theta0 in np.arange(theta_min, theta_max, d_theta):
            T0 = T_min

            theta_const_qs_K = [theta0+273.15,]
            T_arr_K = [T0+273.15,]

            Tk = T_arr_K[0]
            P = 1000*( Tk/theta_const_qs_K[-1])**3.5
            qs0 = 0.622*esat(Tk)/P

            for i in range(int(T_max-T_min)):
                Tk += 1.0

                P = 1000*( Tk/theta_const_qs_K[-1])**3.5
                qs = 0.622*esat(Tk)/P

                dtheta = f_dtheta_dq(theta=theta_const_qs_K[-1], T=Tk)*(qs-qs0)
                theta_new_K = theta_const_qs_K[-1] + dtheta

                theta_const_qs_K.append(theta_new_K)
                T_arr_K.append(Tk)
                qs0 = qs

            T = np.array(T_arr_K) - 273.15
            theta_const_qs = np.array(theta_const_qs_K) - 273.15

            self.ax1.plot(*self._tf(T, theta_const_qs), linestyle='--', color='black')

tephigram = Tephigram()

sounding = np.loadtxt('sounding_example.dat', unpack=True)
P = sounding[0]*100.
T = sounding[2]+273.15
rh = sounding[4]

#from pycfd.reference.atmospheric_flow import stratification_profiles

#profile = stratification_profiles.getStandardIsentropicAtmosphere()
#z = [np.linspace(0., 10000., 100)]
#P = profile.p(z)
#T = profile.temp(z)

tephigram.plot_sounding(P=P, T=T, r=rh)
raw_input()


