import numpy as np
import matplotlib.pyplot as plot
import sys
import os
import argparse

from pycfd.reference.atmospheric_flow import stratification_profiles
from tephigram_python import Tephigram

from singlecloud.utils.profile_integration import find_LFC, find_LCL, LCLNotFoundException, find_EL

from pyclouds.common import Var

from . import profile as default_profile

def plot_test_parcel(profile, dqv0=0., dT0=0.0):
    plot.ioff()

    F__g_to_LCL = find_LCL(profile, dqv0=dqv0, dT0=dT0)
    p__g_to_LCL, T__g_to_LCL, z_LCL = F__g_to_LCL[:,Var.p], F__g_to_LCL[:,Var.T], F__g_to_LCL[-1,Var.z]
    print 'z_LCL=', z_LCL, 'T_LCL=', F__g_to_LCL[-1,Var.T]

    F__g_to_LFC, e_CIN = find_LFC(profile, dqv0=dqv0, dT0=dT0, raise_exception=False, calc_CIN=True)
    p__g_to_LFC, T__g_to_LFC, z_LFC = F__g_to_LFC[:,Var.p], F__g_to_LFC[:,Var.T], F__g_to_LFC[-1,Var.z]
    print "CIN=", e_CIN, 'z_LFC=', z_LFC, "w(CIN)=", np.sqrt(e_CIN)

    F__g_to_EL, e_CAPE = find_EL(profile, dqv0=dqv0, dT0=dT0, raise_exception=False, calc_CAPE=True)
    p__g_to_EL, T__g_to_EL, z_EL = F__g_to_EL[:,Var.p], F__g_to_EL[:,Var.T], F__g_to_EL[-1,Var.z]

    z = np.linspace(0, 4e3, 100)
    T = profile.temp(z)
    p = profile.p(z)

    t = Tephigram(y_range=(-10, 70))
    t.plot_temp(P=p/100., T=T-273.15, with_height_markers=z, marker_interval=500.)
    
    RH = profile.rel_humidity(z)
    t.plot_RH(P=p/100., T=T-273.15, RH=RH)

    # l, = t.plot_temp(P=p__g_to_LCL/100., T=T__g_to_LCL-273.15)
    # l.set_marker('x')
    # l.set_color('grey')
    # l.set_linewidth(4)
    # l, = t.plot_temp(P=p__g_to_LFC/100., T=T__g_to_LFC-273.15)
    # l.set_marker('x')
    # l.set_color('grey')
    # l.set_linewidth(4)
    l, = t.plot_temp(P=p__g_to_EL/100., T=T__g_to_EL-273.15)
    l.set_marker('.')
    l.set_color('grey')
    l.set_linewidth(4)


    print "qv0=", profile.q_v0, "T0=", profile.T0
    plot.show()


if __name__ == "__main__":
    argparser = argparse.ArgumentParser(__doc__)
    argparser.add_argument('--dT', type=float, default=0.0)
    argparser.add_argument('--dqv', type=float, default=0.0)
    
    args = argparser.parse_args()

    profile = default_profile

    plot_test_parcel(profile=profile, dqv0=args.dqv, dT0=args.dT)
