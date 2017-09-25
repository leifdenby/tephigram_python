import numpy as np
import matplotlib.pyplot as plot
import math

from matplotlib.transforms import Affine2D

from mpl_toolkits.axisartist.grid_finder import FixedLocator, MaxNLocator, \
     DictFormatter

from mpl_toolkits.axisartist.grid_helper_curvelinear import GridHelperCurveLinear
from mpl_toolkits.axisartist import Subplot

import profile_integration


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

REF_LINES_ALPHA = 0.4

T1=-10.
T2=40.
theta1=-12.

plot_sat_adiabat_labels = False

rotation_angle = 45.
rotation_origin = [T1, theta1]

theta_min, theta_max, d_theta = -20., 110., 10.
T_min, T_max,d_T = -100., 40., 10.
P_min, P_max, d_P = 200., 1000., 100.
qs_DEFAULT = 1.e-3*np.array([30., 20., 15., 10., 7., 5., 3., 2., 1.0, 0.7,])

T_ticks_DEFAULT = np.arange(T_min, T_max, d_T)

class Tephigram:
    def __init__(self, fig=None, subplotshape=None, plot_default_lines=True, y_range=(-10, 120), x_range=(-20, 60), with_labels=True, height_function=None, T_ticks=T_ticks_DEFAULT):
        self.y_range = y_range
        self.x_range = x_range

        self.height_function = height_function

        if fig is None:
            fig = plot.figure(figsize=(7,10))
        self.fig = fig
        self.setup_axes1(fig=fig, subplotshape=subplotshape, T_ticks=T_ticks)
        self.plot = plot
        self.lines = {}

        self.with_labels = with_labels

        self.plot_sat_adiabat_labels = self.with_labels
        self.plot_annotations = self.with_labels

        self.plotted_lines = []
        if plot_default_lines:
            self.plotted_lines += self.plot_temp_lines()
            self.plotted_lines += self.plot_pot_temp_lines()
            self.plotted_lines += self.plot_pressure_lines()
            self.plotted_lines += self.plot_qs_lines()
            self.plotted_lines += self.plot_sat_adiabats()

        self.ax1.set_yticklabels([])


    def _save_lines(self, line_type, lines):
        if not line_type in self.lines:
            self.lines[line_type] = []
        self.lines[line_type] += lines

    def savefig(self, filename):
        self.fig.savefig(filename)

    def plot_sounding(self, P, T, T_dp):
        """
        Input sounding defined by P, T and T_dp

        Expected units:
            T [C], P [hPa], T_dp [C]
        """
        lines = self.plot_temp(P=P, T=T)

        lines = self.plot_temp(P=P, T=T_dp, temp_type="moisture", label="moisture", marker='o', color='green')
        self._save_lines('moisture', lines)

    def plot_temp(self, P, T, color='red', marker='o', label='temperature', with_height_markers=[], marker_interval=None, temp_type="temp", **kwargs):
        """
        Input sounding defined by P and T

        Expected units:
            T [C], P [hPa]
        """
        theta = self.f_theta(P=P, T=(T+273.15))

        x, y = self._tf(T, theta-273.15)

        lines = self.ax1.plot(*self._tf(T, theta-273.15), marker=marker, color=color, label=label, **kwargs)

        if len(with_height_markers) != 0:
            if len(with_height_markers) != len(x):
                raise Exception("The array of height markers must have the same number of points as the pressure and temperature profiles")
            else:
                z = with_height_markers
                n_skip = np.count_nonzero(z < z.min() + marker_interval)
                for z_, x_, y_ in zip(with_height_markers, x, y)[::n_skip]:
                    t = 'z=%skm' % float('%.2g' % (float(z_)/1e3))
                    self.fig.gca().annotate(t,
                            xy = (x_, y_), xytext = (100, 0),
                            textcoords = 'offset points', ha = 'right', va = 'bottom',
                            bbox = dict(boxstyle = 'round,pad=0.5', fc = 'white', alpha = 0.8),
                            arrowprops = dict(arrowstyle = '->', connectionstyle = 'arc3,rad=0'))

        self._save_lines(temp_type, lines)

        return lines

    def plot_RH(self, P, T, RH, color='green', label='moisture', marker='o', moisture_type="moisture", **kwargs):
        """
        dew-point equation source: http://andrew.rsmas.miami.edu/bmcnoldy/Humidity.html
        """
        T_dp = 243.04*(np.log(RH)+((17.625*T)/(243.04+T)))/(17.625-np.log(RH)-((17.625*T)/(243.04+T)))

        theta = self.f_theta(P=P, T=(T_dp+273.15))

        lines = self.ax1.plot(*self._tf(T_dp, theta-273.15), marker=marker, color=color, label=label, **kwargs)

        self._save_lines(moisture_type, lines)

        return lines


    def setup_axes1(self, fig, T_ticks, subplotshape=None):
        """
        A simple one.
        """
        deg = -45.
        self.tr = Affine2D().rotate_deg(deg)

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

        # ax1.set_aspect(1.)
        x_min, x_max = self.x_range
        ax1.set_xlim(x_min, x_max)
        ax1.set_ylim(*self.y_range)
        ax1.set_xlabel('Temperature [C]')

        ax1.set_aspect(1)

        #ax1.axis["t"]=ax1.new_floating_axis(0, 0.)
        #T_axis = ax1.axis['t']
        #theta_axis = ax1.axis["t2"]=ax1.new_floating_axis(1, 0.)
        
        # plot.draw()
        # plot.show()
        self.ax1 = ax1

    def _tf(self, T, theta):
        """
        Rotate temp and potential temp into the cartesian basis.
        """
        return self.tr.transform(np.array([T, theta]).T).T

    def plot_temp_lines(self, alpha=REF_LINES_ALPHA, **kwargs):
        T_ = np.arange(T_min, T_max+0.1, d_T)
        lines = []
        for T in T_:
            lines += self.ax1.plot(*self._tf([T, T], [theta_min, theta_max]  ), linestyle=':', color='red', label='const. temp', alpha=alpha, **kwargs)
        return lines

    def plot_pot_temp_lines(self, alpha=REF_LINES_ALPHA, **kwargs):
        theta_ = np.arange(theta_min, theta_max, d_theta)

        lines = []
        for theta in theta_:
            lines += self.ax1.plot(*self._tf([T_min, T_max], [theta, theta]  ), linestyle=':', color='green', label='dry adiabat', alpha=alpha, **kwargs)

        self._save_lines('pot_temp_ref', lines)

        return lines

    def plot_pressure_lines(self, p=[], **kwargs):
        """
        theta = T*(1000/p)^0.286
        """
        if len(p) == 0:
            P_ = np.arange(P_min, P_max+0.1, d_P)
        else:
            P_ = np.array(p)/100.

        lines = []
        for P in P_:
            lines += self.plot_pressure_line(P=P, **kwargs)

        return lines

    def plot_pressure_line(self, P, p_label_format='{p:.0f}hPa', z_label_format=None, label_inside=False, label_fontsize=12, alpha=REF_LINES_ALPHA, **kwargs):
        """
        plot line at pressure `P` (in hPa)
        """
        x0 = None
        T__ = np.linspace(T_min, T_max+10., 1000.)
        theta_constP = -273.15 + self.f_theta(P, T__+273.15)

        x, y = self._tf(T__, theta_constP)
        lines = self.ax1.plot(x, y, linestyle=':', color='blue', alpha=alpha, **kwargs)

        xlim = plot.gca().get_xlim()
        x_min, x_max = xlim
        lx = x_max - x_min
        k = np.argmax((x<x_max)*x)

        if x0 is None:
            x0 = x_max + 1.0
        y0 = y[k]

        if (p_label_format is not None or z_label_format is not None) and y0 < self.y_range[1]:
            label_parts = [p_label_format.format(p=P)]

            if not z_label_format is None:
                if self.height_function is None:
                    raise Exception("To include the height on the pressure labels a `height_function` must be passed when creating a tephigram")
                else:
                    label_parts.append(z_label_format.format(z=float(self.height_function(P*100.))))

            label = "\n".join(label_parts)

            if label_inside:
                plot.text(x_max - 0.02*lx, y0, label , color='blue', horizontalalignment='right', size=label_fontsize)
            else:
                plot.text(x0, y0, label , color='blue', size=label_fontsize, verticalalignment='center')

        self._save_lines('p', lines)

        return lines

    def f_theta(self, P, T):
        return T*(P/P_max)**-0.286

    def plot_qs_lines(self, qs_=qs_DEFAULT, include_labels=True, alpha=REF_LINES_ALPHA, **kwargs):
        """
        Make the saturated specific humidity curves
          q_s = 0.622*e_s/p
        """

        T_ = np.linspace(T_min, T_max, 5000)

        lines = []
        for qs in qs_:
            Tk = T_ + 273.15
            Pqs = 0.622*esat(Tk)/qs

            theta = -273.15 + self.f_theta(P=Pqs,T=Tk)

            x, y = self._tf(T_, theta)

            lines += self.ax1.plot(x, y, linestyle='--', color='purple', label=r"$q_{sat}$ const.", alpha=alpha, **kwargs)

            x_min, x_max = plot.gca().get_xlim()
            y_min, y_max = plot.gca().get_ylim()
            k = np.argmin((y>y_min)*y)

            x0 = x[k] 
            y0 = y[k] + 1.5

            if include_labels and x0 > x_min:
                bbox = dict(facecolor='white', edgecolor='white', alpha=0.5)
                units_str = ""
                if qs == qs_[0]:
                    units_str = "g/kg"
                plot.text(x0, y0, "%g%s" % (qs*1000., units_str), color='purple', bbox=bbox)

        self._save_lines('qs_ref', lines)
        return lines


    def plot_sat_adiabats(self, include_labels=True, alpha=REF_LINES_ALPHA, **kwargs):
        """
        Use dtheta = -L*theta/(Cp*T) dqs and integrate from T=-40 C
        """
        def f_dtheta_dq(theta, T):
            L = 2500.08*1.e3 # J/kg
            Cp = 1.004*1.e3 # J/kg/K
            return -L*theta/(Cp*T)

        x_min, x_max = self.x_range

        lines = []
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
            lines += self.ax1.plot(x, y, linestyle='-.', color='black', label='moist adiabat', alpha=alpha, **kwargs)

            k = np.argmin(np.abs(np.array(P_arr)-1000.))
            T_at_1000 = T[k]

            kk = np.argmin((x>x_min)*x - (y<self.y_range[1])*y)

            x0 = x[kk] + 1.
            y0 = y[kk] - 2.

            if include_labels and x[kk] > x_min:
                self.ax1.text(x0, y0, "%gC" % T_at_1000, color='black')

        self._save_lines('sat_adiabat_ref', lines)

        return lines


    def plot_legend(self, lines=[], include_ref_lines=True, ncol=3):
        if len(lines) == 0:
            LINES_NAMES = ['temp', 'moisture', 'test_parcel_temp', 'test_parcel_moisture']
            lines = [self.lines[ln][0] for ln in LINES_NAMES if ln in self.lines]

        if include_ref_lines:
            LINES_NAMES = ['sat_adiabat_ref', 'qs_ref', 'pot_temp_ref', 'temp_ref']
            lines += [self.lines[ln][0] for ln in LINES_NAMES if ln in self.lines]

        # labels = [l.get_label() for l in self.lines]
        # plot.legend(self.lines, labels, loc='upper left', prop={'size': 10})

        plot.tight_layout()
        plot.subplots_adjust(bottom=0.17)

        plot.figlegend(lines, [l.get_label() for l in lines], loc='lower center', frameon=False, ncol=ncol)

    def plot_test_parcel(self, z, P, T, RH, dT0=0.1, dqv0=0.0001):
        p = P*100.
        Tk = T + 273.15
        parcel_integrator = profile_integration.ThermodynamicParcelIntegration(z=z, T=Tk, p=p, rel_humid=RH)
        Var = profile_integration.Var


        F_parcel, parcel_info = parcel_integrator.find_EL(dqv0=dqv0, dT0=dT0, raise_exception=False, include_info=True, z0=z[0])
        p_parcel, T_parcel, qv_parcel = F_parcel[:,Var.p], F_parcel[:,Var.T], F_parcel[:,Var.q_v]

        rh_parcel = qv_parcel/parcel_integrator.pv_sat.qv_sat(p=p_parcel, T=T_parcel)

        lines = []

        # l, = self.plot_temp(P=p__g_to_LCL/100., T=T__g_to_LCL-273.15, label="test parcel", temp_type="test_parcel_temp")
        lines += self.plot_temp(P=p_parcel/100., T=T_parcel-273.15, label="test parcel temp", 
                temp_type="test_parcel_temp", color='grey', linewidth=4, marker='',
        )
        lines += self.plot_RH(P=p_parcel/100., T=T_parcel-273.15, RH=rh_parcel, label="test parcel moisture", 
                moisture_type="test_parcel_moisture", color='grey', linewidth=4, marker='', linestyle=':'
        )

        parcel_info_str = r"$e_{{CIN}}={e_CIN:.0f}$J/kg".format(e_CIN=parcel_info.e_CIN)\
                + "\n" + r"$e_{{CAPE}}={e_CAPE:.0f}$J/kg".format(e_CAPE=parcel_info.e_CAPE)

        plot.text(x=1.0, y=1.0, s=parcel_info_str, horizontalalignment='right', 
                verticalalignment='top', transform=self.plot.gca().transAxes,
                bbox=dict(fc="white"))

        return parcel_info


if __name__ == "__main__":
    tephigram = Tephigram()

    sounding = np.loadtxt('sounding_example.dat', unpack=True)
    P = sounding[0]
    T = sounding[2]
    T_dp = sounding[3]
