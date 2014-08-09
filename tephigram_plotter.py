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

plot_sat_adiabat_labels = False

rotation_angle = 45.
rotation_origin = [T1, theta1]

theta_min, theta_max, d_theta = -20., 110., 10.
T_min, T_max,d_T = -100., 40., 10.
P_min, P_max, d_P = 200., 1000., 100.
qs_ = 1.e-3*np.array([30., 20., 15., 10., 7., 5., 3., 2., 15., 1.0, 0.7,])

x_min = -40.
x_max = 60.0
y_min = -10.
y_max = 80.

T_ = np.arange(T_min, T_max+0.1, d_T)

class Tephigram:
    def __init__(self, fig=None, subplotshape=None):
        if fig is None:
            fig = plot.figure(figsize=(10,10))
        self.fig = fig
        self.setup_axes1(fig, subplotshape)
        self.plot = plot
        self.lines = []

        self.plot_temp_lines()
        self.plot_pot_temp_lines()
        self.plot_pressure_lines()
        self.plot_qs_lines()
        self.plot_sat_adiabats()


    def savefig(self, filename):
        self.fig.savefig(filename)

    def plot_sounding(self, P, T, T_dp):
        """
        Input sounding defined by P, T and T_dp

        Expected units:
            T [C], P [hPa], T_dp [C]
        """
        theta = self.f_theta(P=P, T=(T+273.15))

        self.ax1.plot(*self._tf(T, theta-273.15), marker='o', color='red')

        theta = self.f_theta(P=P, T=(T_dp+273.15))
        self.ax1.plot(*self._tf(T_dp, theta-273.15), marker='o', color='green')

    def plot_temp(self, P, T, color='red'):
        """
        Input sounding defined by P and T

        Expected units:
            T [C], P [hPa]
        """
        theta = self.f_theta(P=P, T=(T+273.15))

        x, y = self._tf(T, theta-273.15)

        self.ax1.plot(*self._tf(T, theta-273.15), marker='o', color=color)

        return x, y

    def plot_RH(self, P, T, RH, color='green'):
        """
        dew-point equation source: http://andrew.rsmas.miami.edu/bmcnoldy/Humidity.html
        """
        T_dp = 243.04*(np.log(RH)+((17.625*T)/(243.04+T)))/(17.625-np.log(RH)-((17.625*T)/(243.04+T)))

        theta = self.f_theta(P=P, T=(T_dp+273.15))
        self.ax1.plot(*self._tf(T_dp, theta-273.15), marker='o', color='green')

    def setup_axes1(self, fig, subplotshape=None):
        """
        A simple one.
        """
        deg = -45.
        self.tr = Affine2D().rotate_deg(deg)

        T_ticks = np.arange(T_min, T_max, d_T)
        theta_ticks = [] #np.arange(theta_min, theta_max, d_T)

        grid_helper = GridHelperCurveLinear(self.tr, grid_locator1=FixedLocator(T_ticks), grid_locator2=FixedLocator(theta_ticks))

        if subplotshape is None:
            subplotshape = (1,1,1)

        ax1 = Subplot(fig, *subplotshape, grid_helper=grid_helper)
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
        ax1.set_ylim(y_min, y_max)
        ax1.set_xlabel('Temperature [C]')

        #ax1.axis["t"]=ax1.new_floating_axis(0, 0.)
        #T_axis = ax1.axis['t']
        #theta_axis = ax1.axis["t2"]=ax1.new_floating_axis(1, 0.)
        
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
            self.ax1.plot(*self._tf([T, T], [theta_min, theta_max]  ), linestyle=':', color='red')

    def plot_pot_temp_lines(self):
        theta_ = np.arange(theta_min, theta_max, d_theta)
        for theta in theta_:
            line = self.ax1.plot(*self._tf([T_min, T_max], [theta, theta]  ), linestyle=':', color='green', label='dry adiabat')

        self.lines += line

    def plot_pressure_lines(self):
        """
        theta = T*(1000/p)^0.286
        """
        for P in np.arange(P_min, P_max+0.1, d_P):
            T__ = np.linspace(T_min, T_max, 1000.)
            theta_constP = -273.15 + self.f_theta(P, T__+273.15)

            x, y = self._tf(T__, theta_constP)
            self.ax1.plot(x, y, linestyle=':', color='blue')

            k = np.argmax((x<x_max)*x)

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

        T_ = np.linspace(T_min, T_max, 1000)

        for qs in qs_:
            Tk = T_ + 273.15
            Pqs = 0.622*esat(Tk)/qs

            theta = -273.15 + self.f_theta(P=Pqs,T=Tk)

            x, y = self._tf(T_, theta)

            self.ax1.plot(x, y, linestyle='--', color='purple')

            k = np.argmin((y>y_min)*y)

            x0 = x[k] - 1.
            y0 = y[k] 

            plot.text(x0, y0, "%g" % (qs*1000.), color='purple')

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
            P_arr = []

            Tk = T_arr_K[0]
            P = 1000*( Tk/theta_const_qs_K[-1])**3.5
            P_arr.append(P)
            qs0 = 0.622*esat(Tk)/P

            for i in range(int(T_max-T_min)):
                Tk += 1.0

                P = 1000*( Tk/theta_const_qs_K[-1])**3.5
                qs = 0.622*esat(Tk)/P

                dtheta = f_dtheta_dq(theta=theta_const_qs_K[-1], T=Tk)*(qs-qs0)
                theta_new_K = theta_const_qs_K[-1] + dtheta

                theta_const_qs_K.append(theta_new_K)
                T_arr_K.append(Tk)
                P_arr.append(P)
                qs0 = qs

            T = np.array(T_arr_K) - 273.15
            theta_const_qs = np.array(theta_const_qs_K) - 273.15
            
            x, y = self._tf(T, theta_const_qs)
            line = self.ax1.plot(x, y, linestyle='-.', color='black', label='moist adiabat')

            k = np.argmin(np.abs(np.array(P_arr)-1000.))
            T_at_1000 = T[k]

            kk = np.argmin((x>x_min)*x - (y<y_max)*y)

            x0 = x[kk] + 1.
            y0 = y[kk] - 2.

            if plot_sat_adiabat_labels:
                self.ax1.text(x0, y0, "%gC" % T_at_1000, color='black')

        self.lines += line

    def plot_legend(self):
        labels = [l.get_label() for l in self.lines]
        plot.legend(self.lines, labels, loc='upper left', prop={'size': 10})
        

if __name__ == "__main__":
    tephigram = Tephigram()

    sounding = np.loadtxt('sounding_example.dat', unpack=True)
    P = sounding[0]
    T = sounding[2]
    T_dp = sounding[3]

    tephigram.plot_sounding(P=P, T=T, T_dp=T_dp)

    RH = sounding[4]
    T_dp2 = 243.04*(np.log(RH/100)+((17.625*T)/(243.04+T)))/(17.625-np.log(RH/100)-((17.625*T)/(243.04+T)))
    tephigram.plot_sounding(P=P, T=T, T_dp=T_dp2)

    tephigram.savefig('tephigram_example.png')


    plot.figure()
    z = sounding[1]

    plot.plot(T_dp2, z, label='calc')
    plot.plot(T_dp, z)
    plot.legend()
    raw_input()
