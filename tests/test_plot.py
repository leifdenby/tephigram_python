import numpy as np
import matplotlib.pyplot as plot

from tephigram_python import Tephigram

def test_plot_dewpoint_temp():
    tephigram = Tephigram()

    sounding = np.loadtxt('examples/sounding_example.dat', unpack=True)
    P = sounding[0]
    T = sounding[2]
    T_dp = sounding[3]

    tephigram.plot_sounding(P=P, T=T, T_dp=T_dp)
    # tephigram.savefig('tephigram_example_dewpoint_temp.png')


def test_plot_rh():
    tephigram = Tephigram()

    sounding = np.loadtxt('examples/sounding_example.dat', unpack=True)
    P = sounding[0]
    T = sounding[2]
    RH = sounding[4]

    tephigram.plot_temp(P=P, T=T)
    tephigram.plot_RH(P=P, T=T, RH=RH)

    # tephigram.savefig('tephigram_example_rh.png')
