"""
Routines for calculating the dry and saturated ascent of a lifted parcel. Also
contains routines for CAPE and CIN calculation.
"""
import numpy as np
import scipy.interpolate
from collections import namedtuple

from attrdict import AttrDict

import parameterisations

default_constants = AttrDict({
    "R_d": 287.05,
    "R_v": 461.51,
    "L_v": 2.5008e6,
    "L_s": 2.8345e6,
    "cp_d": 1005.46,
    "cv_d": 717.60, # J/kg/K
    "cp_v": 1859.0,
    "cv_v": 1402.5, # J/kg/K
    "cp_l": 4183.0,
    "cv_l": 4183.0, # same as cp as liquid is assumed incompressible
    "rho_l": 1000.,
    "rho_i": 500.,
    "g": 9.80665,
    # constants for calculating saturation vapour pressure with Teten's formula
    "pv_sat": AttrDict({
        "p0vs": 611.2,
        "a0_lq": 17.67,
        "a1_lq": -32.19,
        "a0_ice": 22.587,
        "a1_ice": 0.7,
    })
})


def make_related_constants(constants):
    if 'cp_d' in constants and 'cp_v' in constants and not 'R_d' in constants:
        constants['R_d'] = constants['cp_d'] - constants['cv_d']
        constants['R_v'] = constants['cp_v'] - constants['cv_v']
    if 'cp_i' in constants and not 'cv_i' in constants:
        constants['cv_i'] = constants['cp_i']
    if 'cp_l' in constants and not 'cv_l' in constants:
        constants['cv_l'] = constants['cp_l']
    if all([v in constants for v in ('L_s', 'L_v')]) and not 'L_f' in constants:
        constants['L_f'] = constants['L_s'] - constants['L_v']
    return constants


class Var:
    """
    For indexing into array representing full thermodynamic state at given height
    """
    w = 0
    T = 1
    q_v = 2
    q_l = 3
    q_i = 4
    z = 5
    p = 6
    dq_pr = 7 # for storing increments in condensate when integrating pseudoadiabatically

    names = ['w', 'T', 'q_v', 'q_l', 'q_i', 'z', 'p', 'dq_pr']
    NUM = len(names)

    @staticmethod
    def print_formatted(v, formatting='%g'):
        print ",\t".join([("%s=" + formatting) % (Var.names[i], v[i]) for i in range(Var.NUM)])

    @staticmethod
    def repr(v, formatting='%g', skip=[]):
        units = { 'w': 'm/s', 'T': 'K', 'z': 'm', 'p': 'Pa',
                }
        return ", ".join([r"$%s=%g%s$" % (Var.names[i], v[i], units.get(Var.names[i], '')) for i in range(Var.NUM) if not Var.names[i] in skip])

    @staticmethod
    def make_state(**kwargs):
        s = np.zeros((Var.NUM))
        for k, v in kwargs.items():
            s[getattr(Var, k)] = v

        return s


class DiscreteProfile():
    """
    Wrapper for discretely defined ambient profile which exposes the same
    interface as the general profiles but uses interpolation internally
    """

    def __init__(self, z, **kwargs):
        self.vars = kwargs
        self.z = z

        if 'T' in kwargs and not 'temp' in kwargs:
            self.vars['temp'] = kwargs['T']

    def __getattr__(self, name):
        if name in self.__dict__:
            return self.__dict__[name]
        elif name in self.vars:
            y_discrete = self.vars[name]

            return lambda z: scipy.interpolate.interp1d(x=self.z, y=y_discrete)(z)
        else:
            raise AttributeError("Can't find variable `{}`".format(name))


class ThermodynamicParcelIntegration():
    class LCLNotFoundException(Exception):
        pass

    IntegrationInfo = namedtuple("IntegrationInfo", "z_LCL z_LFC z_EL e_CIN e_CAPE")

    def __init__(self, z, p, T, rel_humid=None, qv=None, constants=default_constants):
        self.constants = make_related_constants(constants)
        self.pv_sat = parameterisations.SaturationVapourPressure(constants)

        if sum([not qv is None, not rel_humid is None]) != 1:
            raise Exception("""The moisture profile must be defined through
                    either the relative humidity (`rel_humid`) or water vapour
                    specific concentration (`qv`)""")
        elif not rel_humid is None:
            qv = rel_humid*self.pv_sat.qv_sat(T=T, p=p)
        elif not qv is None:
            rel_humid = qv/self.pv_sat.qv_sat(T=T, p=p)
        else:
            raise NotImplementedError

        if not np.all(rel_humid < 1.0):
            raise NotImplementedError("Parcel integration can't be done with 100% relative humidity in column")

        profile = DiscreteProfile(temp=T, rel_humid=rel_humid, p=p, qv=qv, z=z)

        self.profile = profile

    def _calc_adjusted_state(self, F, iterations, constraint='isobaric'):
        """
        Calculate approximation saturation state using first-order Taylor
        expansion of the moist enthalpy equation at constant pressure:

        cp*dT = -L_v*dq

        Following suggestions in notes by Peter Bechtold
        """
        cp_d = self.constants.cp_d
        cp_v = self.constants.cp_v
        cv_d = self.constants.cv_d
        cv_v = self.constants.cv_v
        L_v = self.constants.L_v

        # only needed for isometric integration, where we use equation of state
        # to update pressure
        rho_l = self.constants.rho_l
        rho_i = self.constants.rho_i
        R_d = self.constants.R_d
        R_v = self.constants.R_v

        qv = F[Var.q_v]
        ql = F[Var.q_l]
        qi = F[Var.q_i]
        qd = 1. - qv - ql - qi

        T = F[Var.T]
        p = F[Var.p]
        qv = F[Var.q_v]

        for n in range(iterations):
            rho = 1.0/((qd*R_d + qv*R_v)*T/p + ql/rho_l + qi/rho_i)

            if constraint == 'isometric':
                c_m = cv_d*qd + (ql + qi + qv)*cv_v
            elif constraint == 'isobaric':
                c_m = cp_d*qd + (ql + qi + qv)*cp_v
            else:
                raise NotImplementedError("Model constraint mode '%s' not implemented" % constraint)

            qv_sat = self.pv_sat.qv_sat(T=T, p=p)
            dqv_sat__dT = self.pv_sat.dqv_sat__dT(T=T, p=p)

            dT = L_v/c_m*(qv - qv_sat)/(1 + L_v/c_m*dqv_sat__dT)

            T = T+dT
            qv = qv - c_m/L_v*dT
            ql = ql + c_m/L_v*dT

            if constraint == 'isometric':
                # rho is unchanged, update pressure
                p = T*(qd*R_d + qv*R_v)/(1./rho - ql/rho_l - qi/rho_i)

        Fs = np.copy(F)
        Fs[Var.q_v] = qv
        Fs[Var.q_l] = ql
        Fs[Var.T] = T
        Fs[Var.p] = p

        return Fs

    def _integrate(self, F0, adjust_function, stopping_criterion, dp=-100.):
        assert dp < 0.0

        F = [F0,]

        constants = self.constants

        g = constants.get('g')

        rho_l = constants.get('rho_l')
        rho_i = constants.get('rho_i')
        R_d = constants.get('R_d')
        R_v = constants.get('R_v')

        cp_v = constants.get('cp_v')
        cp_d = constants.get('cp_d')
        cp_l = constants.get('cp_l')

        qv0 = F0[Var.q_v]
        ql0 = F0[Var.q_l]
        qi0 = F0[Var.q_i]
        qd0 = 1. - ql0 - qi0 - qv0
        
        while not stopping_criterion(F_current=F[-1]):
            p  = F[-1][Var.p]
            T  = F[-1][Var.T]
            qv = F[-1][Var.q_v]
            ql = F[-1][Var.q_l]
            qi = F[-1][Var.q_i]
            qd = 1. - qv - ql - qi

            if qi != 0.0:
                raise NotImplementedError
            cp = qd*cp_d + qv*cp_v + ql*cp_l
            rho_parcel = self.calc_parcel_density(F_parcel=F[-1])

            dz = -dp/(rho_parcel*g)

            F_lifted = np.array(F[-1])
            F_lifted[Var.p] += dp
            F_lifted[Var.z] += dz

            # cool at adiabatic lapse rate, adiabatic expansion (this is not exact!)
            F_lifted[Var.T] -= g/cp*dz

            if adjust_function is not None:
                F_new = adjust_function(F_lifted=F_lifted)
            else:
                F_new = F_lifted

            F.append(F_new)

        return np.array(F)

    def calc_mixture_density(self, p, T, qv, ql, qi):
        R_d   = self.constants.R_d
        R_v   = self.constants.R_v
        rho_l = self.constants.rho_l
        rho_i = self.constants.rho_i

        qd = 1. - qv - ql - qi

        rho_parcel = 1.0/((qd*R_d + qv*R_v)*T/p + ql/rho_l + qi/rho_i)

        return rho_parcel

    def calc_parcel_density(self, F_parcel):
        T  = F_parcel[...,Var.T]
        p  = F_parcel[...,Var.p]
        qv = F_parcel[...,Var.q_v]
        ql = F_parcel[...,Var.q_l]
        qi = F_parcel[...,Var.q_i]

        return self.calc_mixture_density(p=p, T=T, qv=qv, ql=ql, qi=qi)


    def integration_along_moist_adiabat(self, p0, T0, p_min=500.0e2, dp=-10., parcel_sat_calc_num_iterations=6):
        """
        Start from saturation at p0, T0 and integrate along moist adiabat until the pressure drops to p_min
        """

        qv0 = self.pv_sat.qv_sat(T=T0, p=p0)

        F0 = Var.make_state(T=T0, p=p0, q_v=qv0, z=0.0)

        stopping_criterion = lambda F_current: F_current[Var.p] < p_min

        def adjust_function(F_lifted):
            # use moist adjustment at constant pressure to calculate new state
            F_new = self._calc_adjusted_state(F_lifted, iterations=parcel_sat_calc_num_iterations)

            assert F_new[Var.p] == F_lifted[Var.p]
            assert F_new[Var.q_v] < F_lifted[Var.q_v]

            return F_new

        F = self._integrate(F0=F0,
                       adjust_function=adjust_function,
                       stopping_criterion=stopping_criterion)

        return F


    def find_LCL(self, dqv0=0., dT0=0.0, z0=0.0, raise_exception=False):
        rh0 = self.profile.rel_humid(z0)
        T0 = self.profile.temp(z0) + dT0
        p0 = self.profile.p(z0)
        qv0 = rh0*self.pv_sat.qv_sat(T=T0, p=p0) + dqv0

        F0 = Var.make_state(p=p0, T=T0, q_v=qv0, z=0)

        def stopping_criterion(F_current):
            p = F_current[Var.p]

            if p < 500.0e2:
                if raise_exception:
                    raise self.LCLNotFoundException()
                else:
                    return True

            qv = F_current[Var.q_v]
            qv_sat = self.pv_sat.qv_sat(T=F_current[Var.T], p=F_current[Var.p])

            return qv > qv_sat

        return self._integrate(F0=F0, adjust_function=None, stopping_criterion=stopping_criterion)


    def find_LFC(self, dqv0=0.0, dT0=0.0, calc_CIN=False, z0=0.0, raise_exception=True, dp=-100., parcel_sat_calc_num_iterations=6):
        rh0 = self.profile.rel_humid(z0)
        T0 = self.profile.temp(z0) + dT0
        p0 = self.profile.p(z0)
        qv0 = rh0*self.pv_sat.qv_sat(T=T0, p=p0) + dqv0
        F0 = Var.make_state(p=p0, T=T0, q_v=qv0, z=0)

        def stopping_criterion(F_current):
            if F_current[Var.q_l] > 0.0:
                z = F_current[Var.z]

                rho_e = self.calc_mixture_density(
                            p=F_current[Var.p], T=self.profile.temp(z=z),
                            qv=self.profile.qv(z=z), ql=0.0, qi=0.0)

                rho_c = self.calc_parcel_density(F_parcel=F_current)

                # stop once the cloud becomes buoyant, i.e. reaches the level of free convection
                return rho_c < rho_e
            else:
                return False

        def adjust_function(F_lifted):
            T_lifted = F_lifted[Var.T]
            p_lifted = F_lifted[Var.p]

            qv_sat = self.pv_sat.qv_sat(T=T_lifted, p=p_lifted)
            
            if F_lifted[Var.q_v] > qv_sat:
                # use moist adjustment at constant pressure to calculate new state
                F_new = self._calc_adjusted_state(F_lifted, 6)

                assert F_new[Var.p] == F_lifted[Var.p]
                assert F_new[Var.q_v] < F_lifted[Var.q_v]
            else:
                F_new = F_lifted

            return F_new

        F = self._integrate(F0=F0,
                       adjust_function=adjust_function,
                       stopping_criterion=stopping_criterion,
                       dp=dp)

        
        if calc_CIN:
            rho_c = self.calc_parcel_density(F_parcel=F)

            rho_e = self.calc_mixture_density(
                        p=F_current[Var.p], T=self.profile.temp(z=z),
                        qv=self.profile.qv(z=z), ql=0.0, qi=0.0)

            is_buoyant = rho_c < rho_e
            has_condensate = F[:,Var.q_l] != 0.0

            CIN_region = np.logical_and(np.logical_not(is_buoyant), has_condensate)

            # de = drho/rho * g * dz, dp/dz = -rho * g => dz = - dp/(rho*g)
            #    = drho/rho * g * (-dp) / rho / g
            #    = drho/(rho**2) * (-dp)
            de = (rho_c - rho_e)/rho_e**2.*(-dp)

            e_CIN = np.sum(de[CIN_region])

            return F, e_CIN
        else:
            return F

    def find_EL(self, dqv0=0.0, dT0=0.0, pseudo_adiabatic=True, include_info=False, z0=0.0, raise_exception=True, dp=-100., p_min=None,):
        constants = self.constants

        rh0 = self.profile.rel_humid(z0)
        T0 = self.profile.temp(z0) + dT0
        p0 = self.profile.p(z0)
        qv0 = rh0*self.pv_sat.qv_sat(T=T0, p=p0) + dqv0
        F0 = Var.make_state(p=p0, T=T0, q_v=qv0, z=z0)

        class StoppingCriterion():
            def __init__(self, thermodynamic_parcel_integration):
                self.LFC_reached = False
                self.thermodynamic_parcel_integration = thermodynamic_parcel_integration

            def __call__(self, F_current):
                if p_min is not None:
                    if F_current[Var.p] < p_min:
                        return True

                if F_current[Var.q_l] > 0.0 or F_current[Var.dq_pr] > 0.0:
                    z = F_current[Var.z]
                    rho_c = self.thermodynamic_parcel_integration.calc_parcel_density(F_parcel=F_current)

                    rho_e = self.thermodynamic_parcel_integration.calc_mixture_density(
                                p=F_current[Var.p], T=self.thermodynamic_parcel_integration.profile.temp(z=z),
                                qv=self.thermodynamic_parcel_integration.profile.qv(z=z), ql=0.0, qi=0.0)

                    if rho_c < rho_e:
                        if not self.LFC_reached:
                            self.LFC_reached = True
                        return False
                    else:
                        return self.LFC_reached

        def adjust_function(F_lifted):
            T_lifted = F_lifted[Var.T]
            p_lifted = F_lifted[Var.p]

            qv_sat = self.pv_sat.qv_sat(T=T_lifted, p=p_lifted)
            
            if F_lifted[Var.q_v] > qv_sat:
                # use moist adjustment at constant pressure to calculate new state
                F_new = self._calc_adjusted_state(F_lifted, 6)

                assert F_new[Var.p] == F_lifted[Var.p]
                assert F_new[Var.q_v] < F_lifted[Var.q_v]
            else:
                F_new = F_lifted

            if pseudo_adiabatic:
                # put all condensate into precip variable which doesn't contribute to mixture density and heat capacity
                dq_pr = F_new[Var.q_l] + F_new[Var.q_i]
                F_new[Var.dq_pr] = dq_pr
                F_new[Var.q_i] = 0.0
                F_new[Var.q_l] = 0.0

            return F_new

        F = self._integrate(F0=F0,
                       adjust_function=adjust_function,
                       stopping_criterion=StoppingCriterion(self),
                       dp=dp)

        
        if include_info:
            return F, self._calc_parcel_info(F_parcel=F, dp=dp)
        else:
            return F

    def _calc_parcel_info(self, F_parcel, dp):
        rho_c = self.calc_parcel_density(F_parcel=F_parcel)

        z_parcel = F_parcel[:,Var.z]
        rho_e = self.calc_mixture_density(
                p=self.profile.p(z=z_parcel),
                T=self.profile.temp(z=z_parcel),
                qv=self.profile.qv(z=z_parcel), ql=0.0, qi=0.0)

        is_buoyant = rho_c < rho_e
        has_condensate = np.logical_or(F_parcel[:,Var.q_l] != 0.0, F_parcel[:,Var.dq_pr] != 0.0)

        CAPE_region = np.logical_and(is_buoyant, has_condensate)
        z_LFC = z_parcel[np.logical_and(has_condensate, is_buoyant)].min()

        below_CAPE_region = z_parcel < z_LFC

        CIN_region = np.logical_and(np.logical_and(np.logical_not(is_buoyant), has_condensate), below_CAPE_region)

        # de = drho/rho * g * dz, dp/dz = -rho * g => dz = - dp/(rho*g)
        #    = drho/rho * g * (-dp) / rho / g
        #    = drho/(rho**2) * (-dp)
        de = (rho_c - rho_e)/rho_e**2.*(-dp)

        e_CIN = np.sum(de[CIN_region])
        e_CAPE = np.sum(de[CAPE_region])

        z_LCL = z_parcel[has_condensate].min()
        z_EL = z_parcel.max()

        return self.IntegrationInfo(
                e_CIN=e_CIN, e_CAPE=e_CAPE,
                z_LCL=z_LCL, z_LFC=z_LFC, z_EL=z_EL,
        )
