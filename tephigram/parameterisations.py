import numpy as np

from attrdict import AttrDict


class SaturationVapourPressure():
    """
    For calculating saturation vapour pressure and related values using Teten's formula
    """
    AttrDict()

    default_constants = {
        "pv_sat": AttrDict({
            "p0vs": 611.2,  # [Pa]
            "a0_lq": 17.67,
            "a1_lq": -32.19,
            "a0_ice": 22.587,
            "a1_ice": 0.7,
        }),
        "R_d": 287.05,
        "R_v": 461.51,
    }

    def __init__(self, constants=None):
        if constants is None:
            constants = self.default_constants

        self.constants = AttrDict(constants)

    def pv_sat_liquid(self, T):
        p0vs = self.constants.pv_sat.p0vs
        a0_lq = self.constants.pv_sat.a0_lq
        a1_lq = self.constants.pv_sat.a1_lq

        return p0vs*np.exp((a0_lq*(T-273.15)/(T+a1_lq)))

    def pv_sat_ice(self, T):
        p0vs = self.constants.pv_sat.p0vs
        a0_ice = self.constants.pv_sat.a0_ice
        a1_ice = self.constants.pv_sat.a1_ice

        return p0vs*np.exp((a0_ice*(T-273.15)/(T+a1_ice)))

    def pv_sat(self, T):
        v = np.zeros(np.array(T).shape)

        if v.shape == ():
            if T > 273.15:
                return self.pv_sat_liquid(T)
            else:
                return self.pv_sat_ice(T)
        else:
            T = np.array(T)
            idx_liquid = np.array(T) > 273.15
            idx_ice = idx_liquid == False
            v[idx_liquid] = self.pv_sat_liquid(T[idx_liquid])
            v[idx_ice] = self.pv_sat_ice(T[idx_ice])
            return v

    def qv_sat(self, T, p):
        R_v = self.constants.R_v
        R_d = self.constants.R_d

        pv_sat = self.pv_sat(T=T)
        epsilon = R_d/R_v
        qv_sat = (epsilon*pv_sat)/(p-(1.-epsilon)*pv_sat)

        return qv_sat

    def dpsat_dT(self, T):
        if T < 273.15:
            A, B = self.constants.pv_sat.a0_ice, self.constants.pv_sat.a1_ice
        else:
            A, B = self.constants.pv_sat.a0_lq, self.constants.pv_sat.a1_lq

        return self.pv_sat(T=T)*A*(273.15-B)/((T-B)**2.)

    def dqv_sat__dT(self, p, T):
        dpsat_dT = self.dpsat_dT(T=T)

        R_v = self.constants.R_v
        R_d = self.constants.R_d

        pv_sat = self.pv_sat(T=T)

        return R_d/R_v*p*dpsat_dT/((p-pv_sat)**2.)


    def __call__(self, T):
        return self.pv_sat(T)
