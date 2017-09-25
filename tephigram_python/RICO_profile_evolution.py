"""
Looking at the ambient state evolution as a means to explain the variation in
the relationship between maximum cloud-base radius and maximum cloud-top
height.
"""
#import matplotlib
#matplotlib.use('Agg')

import tephigram_python
from cloud_tracking_analysis import CloudData, utils
import matplotlib.pyplot as plot
from matplotlib.gridspec import GridSpec

import numpy as np

# pretty hacky way of setting tephigram setup
tephigram_python.tephigram_plotter.d_theta = 5
tephigram_python.tephigram_plotter.d_T = 5
tephigram_python.tephigram_plotter.P_max = 1000.  # reference pressure in UCLALES is 1000hPa

cloud_data = CloudData('rico_gcss', 'track_1-1500')

# timesteps to extract from
tn_extract = (0 + 8*60*np.arange(3)).astype(int)


fig = plot.figure(figsize=(12, 8))
grids = iter(list(GridSpec(1, 3)))

for n, tn_ in enumerate(tn_extract):
    ax1 = plot.subplot(grids.next())
    ax1.set_xticklabels([])
    ax1.set_yticklabels([])

    profile_fh = cloud_data.get_profile_file_fh()

    T = utils.UCLALES.calc_temperature_profile(profile_fh=profile_fh, tn_=tn_)
    rel_humid = utils.UCLALES.calc_rh_profile(profile_fh=profile_fh, tn_=tn_)
    p = profile_fh.variables['p'][tn_,:]
    z = profile_fh.variables['zt'][:]

    def height_function(p_):
        return z[np.argmin(np.abs(p - p_))]

    tephigram = tephigram_python.Tephigram(fig=fig, subplotshape=(1, 3, n+1), plot_default_lines=False, height_function=height_function)
    plot.gca().set_xlim(25, 38)
    plot.gca().set_ylim(-1, 21)

    tephigram_lines = []
    lines = tephigram.plot_temp_lines()
    tephigram_lines.append(lines[0])

    lines = tephigram.plot_sat_adiabats()
    tephigram_lines.append(lines[0])
    
    lines = tephigram.plot_pot_temp_lines()
    tephigram_lines.append(lines[0])

    lines = tephigram.plot_qs_lines(1.0e-3 *np.array([14., 15., 16., 17., ]), include_labels=True)
    tephigram_lines.append(lines[0])

    for k in 1 + 20*np.arange(0, 7):
        tephigram.plot_pressure_line(p[k]/100, label_inside=True, label_format='{z:.1f}m')

    # add test-parcel plot
    from pycfd.reference.atmospheric_flow.stratification_profiles import DiscreteProfile
    profile = DiscreteProfile(
        T=T,
        p=p,
        rel_humidity=rel_humid,
        z=z,
    )
    from singlecloud.utils import profile_integration
    from pyclouds.common import Var
    dqv0 = -0.5e-3
    dT0 = -0.05

    p_min_parcel = p[140]
    F__g_to_EL = profile_integration.find_EL(profile, dqv0=dqv0, dT0=dT0, raise_exception=False, p_min=p_min_parcel)
    p__g_to_EL, T__g_to_EL, z_EL = F__g_to_EL[:,Var.p], F__g_to_EL[:,Var.T], F__g_to_EL[-1,Var.z]
    l, = tephigram.plot_temp(P=p__g_to_EL/100., T=T__g_to_EL-273.15)
    l.set_marker('.')
    l.set_color('grey')
    l.set_linewidth(4)
    l.set_label("Test-parcel temperature")

    # plot mean profile

    line_T, = tephigram.plot_temp(P=p/100., T=T-273.15, marker_interval=500.)
    line_T.set_marker('.')
    line_T.set_linestyle('')
    line_T.set_linewidth(3)
    line_T.set_label("horz. mean temperature")
    line_RH, = tephigram.plot_RH(P=p/100., T=T-273.15, RH=rel_humid)
    line_RH.set_marker('.')
    line_RH.set_linestyle('')
    line_RH.set_linewidth(3)
    line_RH.set_label("horz. mean water content")


    title = "t={time}min".format(time=tn_)
    if tn_ == 0:
        title += " (initial condition)"
    plot.title(title)

    # stupid matplotlib makes figlegends by column first...
    def alternate_lists(l1, l2):
        result = [None]*(len(l1)+len(l2))
        result[::2] = l1
        result[1::2] = l2
        return result

    lines = alternate_lists(tephigram_lines, [line_T, line_RH, l])

    F__g_to_LCL = profile_integration.find_LCL(profile, dqv0=dqv0, dT0=dT0)
    F__g_to_LFC = profile_integration.find_LFC(profile, dqv0=dqv0, dT0=dT0)
    print "z_LCL={}m, z_LFC={}m".format(F__g_to_LCL[-1,Var.z], F__g_to_LFC[-1,Var.z])

#plot.suptitle("RICO ambient profile evolution")
plot.tight_layout()
plot.subplots_adjust(bottom=0.17)

plot.figlegend(lines, [l.get_label() for l in lines], loc='lower center', frameon=False, ncol=4)

plot.savefig(__file__.replace('.py', '.pdf'))
