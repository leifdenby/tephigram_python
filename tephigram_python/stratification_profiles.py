import warnings
import numpy as np
import scipy.constants
import scipy.optimize
import scipy.interpolate
from pycfd.reference.atmospheric_flow import gas_properties as ref_gas_properties

def getStandardIsothermalAtmosphere():
    gas_properties = ref_gas_properties.AtmosphericAir()
    rho0 = 1.205
    p0 = 101325.0
    dTdz = 0.0  # isothermal
    g = scipy.constants.g
    return HydrostaticallyBalancedAtmosphere(rho0=rho0, p0=p0, dTdz=dTdz, gas_properties=gas_properties, g=g)

def getStandardIsentropicAtmosphere():
    gas_properties = ref_gas_properties.AtmosphericAir()
    rho0 = 1.205
    p0 = 101325.0
    g = scipy.constants.g
    dTdz = -g/gas_properties.cp()
    return HydrostaticallyBalancedAtmosphere(rho0=rho0, p0=p0, dTdz=dTdz, gas_properties=gas_properties, g=g)

def getConstantDensityAtmosphere():
    gas_properties = ref_gas_properties.AtmosphericAir()
    rho0 = 1.205
    p0 = 101325.0
    g = scipy.constants.g
    dTdz = -g/gas_properties.R()
    return HydrostaticallyBalancedAtmosphere(rho0=rho0, p0=p0, dTdz=dTdz, gas_properties=gas_properties, g=g)

def getKleinIsentropicAtmosphere():
    gas_properties = ref_gas_properties.AtmosphericAir()
    rho0 = 1.0
    g = 10.0
    T0 = 300.0
    p0 = rho0*gas_properties.R()*T0
    dTdz = -g/gas_properties.cp()
    return HydrostaticallyBalancedAtmosphere(rho0=rho0, p0=p0, dTdz=dTdz, gas_properties=gas_properties, g=g)

def getIsentropicAtmosphere(theta0=300.0, g=9.81, p0=1.0e6):
    """
    Generate an isentropic atmopsheric profile.
    """
    gas_properties = ref_gas_properties.AtmosphericAir()
    T0 = theta0
    rho0 = p0/T0*1.0/gas_properties.R()
    dTdz = -g/gas_properties.cp()
    return HydrostaticallyBalancedAtmosphere(rho0=rho0, p0=p0, dTdz=dTdz, gas_properties=gas_properties, g=g)


class HydrostaticallyBalancedAtmosphere(object):
    """
    Class for setting a hydrostatically balanced atmosphere with
    an ideal gas and constant lapse rate given as dT/dz.

    p = rho*R/M*T

    R: Unified gas constant
    """
    def __init__(self, rho0, p0, dTdz, gas_properties, g = None):
        self.rho0 = rho0
        self.p0 = p0
        self.dTdz = dTdz
        self.gas_properties = gas_properties
        self.T0 = p0*gas_properties.M/(rho0*scipy.constants.R*1000.0)
        self.pot_temperature0 = self.T0  # p=p0 at surface
        if g is None:
            self.g = scipy.constants.g
        else:
            self.g = g

    def __str__(self):
        return "HydrostaticallyBalancedAtmosphere (rho0=%f, p0=%f, dTdz=%f) with %s" % (self.rho0, self.p0, self.dTdz, str(self.gas_properties))

    def temp(self, pos):
        p = np.array(pos)
        if len(p.shape) > 1:
            z = p[-1]
        else:
            z = p
        return self.T0 + self.dTdz*z

    def rho(self, pos):
        p = np.array(pos)
        if len(p.shape) > 1:
            z = p[-1]
        else:
            z = p

        if self.dTdz == 0.0:
            return self.rho0*np.exp(-z*self.g*self.gas_properties.M/(scipy.constants.R*1000.0*self.T0))
        else:
            alpha = self.g*self.gas_properties.M/(self.dTdz*scipy.constants.R*1000.0)
            return self.rho0*np.power(self.T0, alpha+1.0 )*np.power(self.temp(pos), -alpha - 1.0)

    def drho_dz(self, pos):
        if self.dTdz == 0.0:
            return -self.g*self.gas_properties.M/(scipy.constants.R*1000.0*self.T0)*self.rho(pos)
        else:
            alpha = self.g*self.gas_properties.M/(self.dTdz*scipy.constants.R*1000.0)
            return (-alpha-1.0)*np.power(self.temp(pos), -alpha - 2.0)*self.dTdz

    def p(self, pos):
        return self.rho(pos)*scipy.constants.R*1000.0/self.gas_properties.M*self.temp(pos)

    def pot_temperature(self, pos):
        """
        Calculate the potential temperature at pos.
        """
        return self.temp(pos)*np.power(self.p(pos)/self.p0, -self.gas_properties.kappa())

    def x_velocity(self, pos):
        return 0.0

    def y_velocity(self, pos):
        return 0.0

    def lapseRate(self):
        """
        Calculate lapse rate for using in the CNS-AMR compressible code
        """
        return self.dTdz*scipy.constants.R*1000.0/self.gas_properties.M

class PseudoHydrostaticallyBalancedMoistAtmosphere(HydrostaticallyBalancedAtmosphere):
    """
    Given an atmosphere with non-zero specific concentration of water vapour
    and a fixed lapse rate it is not strictly guaranteed that the atmosphere
    will be stable, and so the hydrostatic assumption is not strictly valid.
    """
    def __init__(self, rho0, p0, dTdz, RH0, dRHdz, gas_properties, g=None):
        super(PseudoHydrostaticallyBalancedMoistAtmosphere, self).__init__(rho0=rho0, p0=p0, dTdz=dTdz, gas_properties=gas_properties, g=g)

        self.dRHdz = dRHdz
        self.RH0 = RH0

    def rel_humidity(self, pos):
        z = pos[-1]
        return self.RH0 + self.dRHdz*z

a = 6.112
b = 12.62
c = 243.5
gamma = lambda T, RH : np.log(RH) + b*T/(c+T)
P_a = lambda T, RH: a*np.exp(gamma(T, RH))
T_dp = lambda T, RH: c*np.log(P_a(T, RH)/a)/(b-np.log(P_a(T, RH)/a))

class AttrDict:
    def __init__(self, d):
        self.__dict__ = d


class LayeredAtmosphere(object):
    def __init__(self, layers):
        self.layer_instances = {}

        # ground state
        z_min = 0.0
        for layer in layers:
            z_max = layer['z_max']
            z = (z_min, z_max)
            layer_instance = AttrDict(layer)
            self.layer_instances[z] = layer_instance

            # calculate the start values of the next layer, remember that this
            # layer is offset.
            z_offset = z_max - z_min
            z_min = z_max

    def _get_values_from_layer(self, variable, pos):

        if type(pos) in [float, np.float, np.float64 ]:
            z = pos

            for (z_min, z_max), layer in self.layer_instances.items():
                if z_min <= z and z <= z_max:
                    f = getattr(layer, variable)
                    return f(z - z_min)

        else:
            pos = np.array(pos)
            if len(pos.shape) > 1:
                z = pos[-1]
            else:
                z = pos

            values = np.zeros(z.shape)

            for (z_min, z_max), layer in self.layer_instances.items():
                idx_in_layer = np.logical_and(z_min <= z, z < z_max)
                f = getattr(layer, variable)
                values[idx_in_layer] = f([z[idx_in_layer] - z_min])

            return values

class LayeredDryAtmosphere(LayeredAtmosphere):
    def __init__(self, layers, rho0=None, p0=None, gas_properties=None):
        self.layers = layers
        if gas_properties is None:
            self.gas_properties = ref_gas_properties.AtmosphericAir()
        else:
            self.gas_properties = gas_properties

        # create an instance of HydroststaticallyBalancedAtmosphere for
        # each layer
        self.layer_instances = {}

        # ground state
        z_min = 0.0
        if rho0 is None:
            rho0 = 1.205
        if p0 is None:
            p0 = 101325.0
        for layer in self.layers:
            z_max = layer['z_max']
            z = (z_min, z_max)
            layer_instance = HydrostaticallyBalancedAtmosphere(rho0=rho0,
                                                               p0=p0,
                                                               dTdz=layer['dTdz'],
                                                               gas_properties=self.gas_properties,
                                                               )
            self.layer_instances[z] = layer_instance

            # calculate the start values of the next layer, remember that this
            # layer is offset.
            z_offset = z_max - z_min
            z_min = z_max
            rho0 = layer_instance.rho([z_offset])
            p0 = layer_instance.p([z_offset])

    def temp(self, pos):
        return self._get_values_from_layer('temp', pos)

    def p(self, pos):
        return self._get_values_from_layer('p', pos)

    def rho(self, pos):
        return self._get_values_from_layer('rho', pos)

class NearIsentropic(HydrostaticallyBalancedAtmosphere):
    """
    This profile is forced a little more stable that isentropic (neutral)
    stability, since ATHAM can't run with a profile that is rigth on neutral.
    """
    def __init__(self, dTdz_offset=1.0e-3):
        gas_properties = ref_gas_properties.AtmosphericAir()
        rho0 = 1.205
        p0 = 101325.0
        g = scipy.constants.g
        dTdz = -g/gas_properties.cp() + dTdz_offset
        super(NearIsentropic, self).__init__(rho0=rho0, p0=p0, dTdz=dTdz, gas_properties=gas_properties, g=g)

        self.dTdz_offset = dTdz_offset

    def __str__(self):
        return "Near-isentropic, dTdz_offset={offset}K/km (dry)".format(offset=self.dTdz_offset*1.e3)

class LayeredStable(LayeredDryAtmosphere):

    def __init__(self):
        gas_properties = ref_gas_properties.AtmosphericAir()
        g = scipy.constants.g
        dTdz = -g/gas_properties.cp()
        rho0 = 1.205
        p0 = 101325.0

        layers = []
        layers.append({'z_max': 2.0e3, 'dTdz': dTdz+1.0e-3,})
        layers.append({'z_max': np.finfo('f').max, 'dTdz': dTdz+5.0e-3})
        super(LayeredStable, self).__init__(layers=layers, rho0=rho0, p0=p0)

    def __str__(self):
        return "Very stable two-layered atmosphere"

class LayeredMoistAtmosphere(object):
    def __init__(self, layers, RH0, RH_min=None, rho0=None, p0=None):
        self.layers = layers
        self.RH_min = RH_min
        self.RH0 = RH0
        self.gas_properties = ref_gas_properties.AtmosphericAir()

        # create an instance of HydrostHydrostaticallyBalancedAtmosphere for
        # each layer
        self.layer_instances = {}

        # ground state
        z_min = 0.0
        if rho0 is None:
            rho0 = 1.205
        if p0 is None:
            p0 = 101325.0
        RH0 = self.RH0
        for layer in self.layers:
            z_max = layer['z_max']
            z = (z_min, z_max)
            layer_instance = PseudoHydrostaticallyBalancedMoistAtmosphere(rho0=rho0,
                                                                    p0=p0,
                                                                    dTdz=layer['dTdz'],
                                                                    dRHdz=layer['dRHdz'],
                                                                    RH0=RH0,
                                                                    gas_properties=self.gas_properties,
                                                                    )
            self.layer_instances[z] = layer_instance

            # calculate the start values of the next layer, remember that this
            # layer is offset.
            z_offset = z_max - z_min
            z_min = z_max
            rho0 = layer_instance.rho([z_offset])
            p0 = layer_instance.p([z_offset])
            RH0 = layer_instance.rel_humidity([z_offset])

    def temp(self, pos):
        return self._get_values_from_layer('temp', pos)

    def p(self, pos):
        return self._get_values_from_layer('p', pos)

    def rho(self, pos):
        return self._get_values_from_layer('rho', pos)

    def _get_values_from_layer(self, variable, pos):

        if len(np.array(pos).shape) > 1:
            z = pos[-1]
        else:
            z = pos

        try:
            values = np.zeros(z.shape)

            for (z_min, z_max), layer in self.layer_instances.items():
                idx_in_layer = np.logical_and(z_min <= z, z < z_max)
                f = getattr(layer, variable)
                values[idx_in_layer] = f([z[idx_in_layer] - z_min])

            return values
        except AttributeError:
            for (z_min, z_max), layer in self.layer_instances.items():
                if z_min <= z and z <= z_max:
                    f = getattr(layer, variable)
                    return f([z - z_min])


    def rel_humidity(self, pos):
        rel_humidity = self._get_values_from_layer('rel_humidity', pos)
        try:
            rel_humidity[rel_humidity < self.RH_min] = self.RH_min
            return rel_humidity
        except TypeError:
            if rel_humidity < self.RH_min:
                return self.RH_min
            else:
                return rel_humidity

    def dew_point(self, pos):
        temp = self.temp(pos)
        rel_humidity = self.rel_humidity(pos)

        return T_dp(temp-273.15, rel_humidity) + 273.15

class Soong1973(LayeredMoistAtmosphere):
    def __init__(self, cloud_base_height=None):
        self.cloud_base_height = cloud_base_height

        layers = []
        dTdz_dry = -10.0e-3  # K/m
        dTdz_moist = -6.0e-3  # K/m

        RH0 = 0.70

        if self.cloud_base_height is not None:
            layer_thickness_0 = self.cloud_base_height # m
            # TODO: This is completely arbitrary, the RH needs lowering that's
            # for sure, not sure how much at this point
            if cloud_base_height == 1600.:
                RH0 = 0.50
            elif cloud_base_height == 1100.:
                RH0 = 0.59
            else:
                raise NotImplementedError
        else:
            layer_thickness_0 = 800.0 # m

        RH_LCL = 0.90

        dRHdz_0 = (RH_LCL - RH0)/layer_thickness_0  # %/m
        dRHdz_1 = -0.075e-3 # %/m

        layers.append({'z_max': layer_thickness_0, 'dTdz': dTdz_dry, 'dRHdz': dRHdz_0})
        layers.append({'z_max': 12800., 'dTdz': dTdz_moist, 'dRHdz': dRHdz_1})
        layers.append({'z_max': np.finfo('f').max, 'dTdz': 0.0, 'dRHdz': 0.0})

        T0 = 25.0 + 273.15 #  [K], from Soong 1973 paper
        p0 = 101325.0 #  [Pa], default value used, Soong 1973 doesn't give a value

        gas_properties = ref_gas_properties.AtmosphericAir()
        rho0 = p0/T0*1.0/gas_properties.R()

        super(Soong1973, self).__init__(layers=layers, RH0=RH0, RH_min=0.3, rho0=rho0, p0=p0)

    def __str__(self):
        if self.cloud_base_height is None:
            return "Soong 1973 layered moist atmosphere"
        else:
            return "Soong 1973 layered moist atmosphere (with modified cloud base at %s)" % self.cloud_base_height

class Soong1973Dry(LayeredDryAtmosphere):
    def __init__(self):
        layers = []
        dTdz_dry = -10.0e-3  # K/m
        dTdz_moist = -6.0e-3  # K/m

        self.g = 10.0

        layer_thickness_0 = 800.0 # m

        layers.append({'z_max': 800.0, 'dTdz': dTdz_dry, })
        layers.append({'z_max': 12800., 'dTdz': dTdz_moist, })
        layers.append({'z_max': np.finfo('f').max, 'dTdz': 0.0, })

        super(Soong1973Dry, self).__init__(layers=layers, )

    def __str__(self):
        return "Soong 1973 layered dry atmosphere"

class SimpleMoistStable(LayeredMoistAtmosphere):
    def __init__(self, cloud_base_height=None, dRHdz=None, RH_LCL=None):
        self.cloud_base_height = cloud_base_height

        layers = []
        dTdz_dry = -8.0e-3  # K/m
        dTdz_moist = -6.0e-3  # K/m

        if self.cloud_base_height is not None:
            layer_thickness_0 = self.cloud_base_height # m
        else:
            layer_thickness_0 = 800.0 # m

        if dRHdz is None:
            dRHdz = -0.2e-3 # %/m

        self.dRHdz = dRHdz

        RH0 = 0.70
        if RH_LCL is None:
            RH_LCL = 0.90

        dRHdz_0 = (RH_LCL - RH0)/layer_thickness_0  # %/m
        dRHdz_1 = self.dRHdz # %/m

        RH_min = 0.3

        layers.append({'z_max': layer_thickness_0, 'dTdz': dTdz_dry, 'dRHdz': dRHdz_0})
        layers.append({'z_max': 12800., 'dTdz': dTdz_moist, 'dRHdz': dRHdz_1})
        layers.append({'z_max': np.finfo('f').max, 'dTdz': 0.0, 'dRHdz': 0.0})

        super(SimpleMoistStable, self).__init__(layers=layers, RH0=RH0, RH_min=RH_min)

    def __str__(self):
        return "Simple stable moist atmosphere based on Soong1973 (dRHdz=%g%%/km)" % (self.dRHdz*1.e3)


class RICO:
    """
    Based on KNMI's synthesis of the RICO field compaign for a LES intercomparison study

    http://projects.knmi.nl/rico/setup3d.html

    OBS: The sea surface temperature is different from the temperature at z=0m,
    which is necessary to have surface fluxes.
    """
    def __init__(self, include_wind=False):
        from pyclouds import parameterisations

        self.include_wind = include_wind
        if include_wind:
            self.u_wind = self._u_wind
            self.v_wind = self._v_wind
        else:
            self.u_wind = lambda z: np.zeros_like(z)
            self.v_wind = lambda z: np.zeros_like(z)

        # surface conditions
        self.ps = 101540.  # [Pa], surface pressure
        self.p0 = 100000. # [Pa], reference pressure
        self.Ts = 299.8  # [K], sea surface temperature

        # constants
        self.Lv = 2.5e6  # [J/kg]
        self.c_p = 1005. # [J/kg/K]
        self.g = 9.81    # [m/s^2]
        self.R_d = 287.  # [J/kg/K]

        # XXX: R_v and cp_v are not given in the RICO test definition on the
        # the KNMI site I will use what I believe are standard values here
        self.R_v = parameterisations.common.default_constants.get('R_v')
        self.cp_v = parameterisations.common.default_constants.get('cp_v')

        self.z_max = 4e3

        self._create_profile()

    def _create_profile(self):
        """ Create a vertical profile that we can interpolate into later. 
        Integrating with the hydrostatic assumption.
        """
        from pyclouds import parameterisations

        parameterisation = parameterisations.SaturationVapourPressure()

        dz = 100.
        R_d = self.R_d

        R_v = self.R_v
        R_d = self.R_d
        cp_d = self.c_p
        cp_v = self.cp_v

        z = 0.0
        p = self.ps

        # Cathy suggested using the liquid water potential temperature as the
        # temperature in the first model level
        T = self.theta_l(0.0)

        profile = []

        n = 0
        while z < self.z_max:
            qt = self.q_t(z)

            # assume no liquid water
            ql = 0.0
            qv = qt

            qd = 1.0 - qt

            theta_l = self.theta_l(z)

            R_l = R_d*qd + R_v*qv
            c_l = cp_d*qd + cp_v*qv

            T = theta_l/((self.p0/p)**(R_l/c_l))
            # T = self.iteratively_find_temp(theta_l=theta_l, p=p, q_t=qt, q_l=ql, T_initial=T)

            rho = 1.0/((qd*R_d + qv*R_v)*T/p) # + 1.0/(ql/rho_l), ql = 0.0

            profile.append((z, rho, p, T))

            # integrate pressure
            z += dz
            p += - rho * self.g * dz

            n += 1

        self._profile = np.array(profile)

    def iteratively_find_temp(self, theta_l, p, q_t, q_l, T_initial):
        R_v = self.R_v
        R_d = self.R_d
        cp_d = self.c_p
        cp_v = self.cp_v

        p0 = self.p0
        omega_l = 1.0
        L_v = self.Lv

        R_l = R_d*(1.0 - q_t) + R_v*q_t
        c_l = cp_d*(1.0 - q_t) + cp_v*q_t

        C = R_l/c_l*np.log(p0/p) + np.log(omega_l) - np.log(theta_l)

        theta_l_f = lambda T: T*(p0/p)**(R_l/c_l)*omega_l*np.exp(- L_v * q_l /( c_l * T))

        f = lambda T: C + np.log(T) - L_v*q_l/(c_l*T)
        T = scipy.optimize.brentq(f, 1., 410.)

        if np.isnan(T):
            raise Exception("Integration failed")

        if not np.abs(theta_l_f(T) - theta_l) < 1.0e-10:
            print theta_l_f(T), theta_l, np.abs(theta_l_f(T) - theta_l) 

        return T

    def _get_value_from_precomputed_profile(self, pos, var_indx):
        z = self._profile[:,0]
        if np.any(z > self.z_max):
            raise Exception("RICO test case is only defined for z < 7km")
        T = self._profile[:,var_indx]
        return np.interp(pos, z, T)

    def rho(self, pos):
        return self._get_value_from_precomputed_profile(pos, 1)

    def p(self, pos):
        return self._get_value_from_precomputed_profile(pos, 2)

    def temp(self, pos):
        return self._get_value_from_precomputed_profile(pos, 3)

    def rel_humidity(self, z):
        q_v = self.q_t(z)
        p = self.p(z)
        T = self.temp(z)

        from pyclouds import parameterisations
        parameterisation = parameterisations.SaturationVapourPressure()

        qv_sat = parameterisation.qv_sat(T=T, p=p)

        return q_v/qv_sat

    def q_t(self, z):
        """ Total water specific concentration [kg/kg]"""
        return self._q_t(z)/1000.

    @np.vectorize
    def _q_t(z):
        if 0 <= z < 740:
            return 16.0 + (13.8 - 16.0) / (740) * z
        if 740 < z < 3260:
            return 13.8 + (2.4 - 13.8) / (3260 - 740)*(z - 740)
        else:
            return 2.4 + (1.8 - 2.4)/(4000 - 3260)*(z - 3260) 

    @np.vectorize
    def theta_l(z):
        """ Liquid water potential temperature [K]"""
        if 0 <= z < 740:
            return 297.9
        else:
            return 297.9 + (317.0 - 297.9)/(4000 - 740) *(z - 740)

    @np.vectorize
    def _u_wind(z):
        """
        u-component of wind
        """
        if z > 0.0:
            return -9.9 + 2.0e-3*z
        else:
            return 0.0


    @np.vectorize
    def _v_wind(z):
        """
        v-component of wind
        """
        if z > 0.0:
            return -3.8 
        else:
            return 0.0

    @np.vectorize
    def ddt_theta_l__ls(z):
        """
        Large Scale Horizontal Liq. Water Pot. Temperature Advection combined
        with Radiative Cooling [K/s] 

        NB: Initial profile contains no liquid water so `temp = pot. temp`
        """
        if z > 0:
            return -2.5 / 86400
        else:
            return 0.0

    @np.vectorize
    def ddt_qv_ls(z):
        """
        Large Scale Horizontal Moisture Advection [(kg/kg)/s]

        NB: not exactly as the KNMI website because we want to return tendencies
        in kg/kg/s, not g/kg/s
        """
        if 0 < z <= 2980:
            return (-1.0 / 86400 + (1.3456/ 86400) * z / 2980) * 1.0e-3
        elif z > 2980:
            return 4.*1.0e-6  * 1.0e-3
        else:
            return 0.0

    @np.vectorize
    def w_subsidence(z):
        """
        Large Scale Subsidence w [m/s] Apply the subsidence on the prognostic fields of q_t, theta_l.
        """
        if 0 < z < 2260:
            return - (0.005/2260) * z
        else:
            return - 0.005

    @np.vectorize
    def tke(z):
        """Initial subgrid profile of subgrid TKE"""
        return  1 - z/4000


    def __str__(self):
        return "RICO, LES test case from KNMI (%s wind)" % ['without', 'with'][self.include_wind]

class RICO_SCM:
    @np.vectorize
    def temp(self, z):
        if 0.0 <= z < 740.:
            return 299.2 + (292.0 - 299.2) / (740) * z 
        elif 740 < z < 4000:
            return 292.0 + (278.0 - 292.0) / (4000 - 740) * (z - 740)
        elif 4000 < z < 15000:
            return 278.0 + (203.0 - 278.0) / (15000 - 4000) * (z - 4000)
        elif 15000 < z < 17500:
            return 203.0 + (194.0 - 203.0) / (17500 - 15000)* (z - 15000)
        elif 17500 < z < 20000:
            return 194.0 + (206.0 - 194.0) / (20000 - 17500)* (z - 17500)
        elif 20000 < z < 60000:
            return 206.0 + (270.0 - 206.0) / (60000 - 20000)* (z - 20000) 
        else:
            raise Exception("z out of range")

    def q_v(self, z):
        """ Water vapour specific concentration in [kg/kg]"""
        return self._q_v(z)/1000.

    @np.vectorize
    def _q_v(z):
        if 0 <= z < 740:
            return 16.0 + (13.8 - 16.0) / (740) * z
        elif 740 < z < 3260:
            return 13.8 + (2.4 - 13.8) / (3260 - 740) * (z - 740)
        elif 3260 < z < 4000:
            return 2.4 + (1.8 - 2.4) / (4000 - 3260) * (z - 3260)
        elif 4000 < z < 9000:
            return 1.8 + (0 - 1.8) / (10000 - 4000) * (z - 4000)
        else:
            return 0.0

    def _create_profile():
        pass


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

            return scipy.interpolate.interp1d(x=self.z, y=y_discrete)

        else:
            raise AttributeError("Can't find variable `{}`".format(name))

class TwoLayerMoistIsentropicPBL():
    """
    Moist sub-saturated well-mixed (isentropic and constant water vapour
    concentration) boundary layer with a layer above with higher lapse-rate"""
    def __init__(self, z_BL, RH0, T0, z_INV, dRHdz_1=-0.1e-3, p0=101325., dTdz_mid=-6.0e-3, constants=None):
        self.z_BL = z_BL
        self.T0 = T0
        self.p0 = p0
        self.RH0 = RH0

        self.z_INV = z_INV

        if constants is None:
            from pyclouds.common import default_constants
            constants = default_constants
            warnings.warn("Using default constants from pyclouds")

        from attrdict import AttrDict
        constants = AttrDict(constants)
        cp_v = constants.cp_v
        cp_d = constants.cp_d
        R_v = constants.R_v
        R_d = constants.R_d
        g = constants.g

        self.constants = constants

        layers = []

        qv_sat_0 = self.qv_sat(z=0.)
        q_v = RH0*qv_sat_0
        q_d = 1. - q_v

        self.q_v0 = q_v

        cp = q_v*cp_v + q_d*cp_d
        R = q_v*R_v + q_d*R_d

        rho0 = p0/(R*T0)
        self.rho0 = rho0

        self.dTdz_BL = -g/cp
        self.dTdz_2 = dTdz_mid  # K/m

        self.dRHdz_1 = dRHdz_1

        # TODO: would be nice to abstract this away to have a better storage
        # for gas-mixture properties
        gas_properties = AttrDict(M=scipy.constants.R*1000./R)

        self.BL_profile = HydrostaticallyBalancedAtmosphere(
            rho0=rho0, p0=p0, dTdz=self.dTdz_BL, gas_properties=gas_properties)

        p_BL_top = self.BL_profile.p(z_BL)
        T_BL_top = self.BL_profile.temp(z_BL)

        rho_BL_top = p_BL_top/(R*T_BL_top)

        # self.MIDDLE_profile = HydrostaticallyBalancedAtmosphere(
            # rho0=rho_BL_top,
            # p0=p_BL_top,
            # dTdz=self.dTdz_2,
            # gas_properties=gas_properties,
        # )

        self.dRHdz_2 = -0.4e-3
        self.__init_profile(p_min=500e2)

        rho_top = self.rho(self.z_INV)
        p_top = self.p(self.z_INV)

        # self.TOP_profile = HydrostaticallyBalancedAtmosphere(
            # rho0=rho_top,
            # p0=p_top,
            # dTdz=0.0,
            # gas_properties=gas_properties,
        # )

    def __init_profile(self, p_min):
        from pyclouds import parameterisations
        qv_sat__f = parameterisations.ParametersationsWithSpecificConstants(self.constants).pv_sat.qv_sat
        R_v = self.constants.get('R_v')
        R_d = self.constants.get('R_d')
        g = self.constants.get('g')

        def rho_f(T, p, qv):
            qd = 1.0 - qv
            rho_inv = (qd*R_d + qv*R_v)*T/p
            return 1.0/rho_inv

        # do numerical integration to take into account that heat capacity changes
        dp = -100. # [Pa]

        p = [self.p0,]
        rho = [self.rho0,]
        z = [0.0, ]

        while p[-1] > p_min:
            p.append(p[-1] + dp)

            # dp/dz = - rho * g
            dz = - dp/(rho[-1]*g)

            z.append(z[-1] + dz)

            T = self.temp(z[-1])
            if z[-1] <= self.z_BL:
                qv = self.q_v0
            else:
                qv = self.rel_humidity(z[-1])*qv_sat__f(T=T, p=p[-1])
            rho.append(rho_f(T=T, p=p[-1], qv=qv))

            self._p = np.array(p)
            self._rho = np.array(rho)
            self._z = np.array(z)


    def qv_sat(self, z):
        from pyclouds import parameterisations
        T = self.temp(z)
        p = self.p(z)
        return parameterisations.ParametersationsWithSpecificConstants(self.constants).pv_sat.qv_sat(T=T, p=p)

    def temp(self, z):
        @np.vectorize
        def f(z):
            if z == 0.0:
                return self.T0
            elif z <= self.z_BL:
                return self.T0 + z*self.dTdz_BL
            elif z <= self.z_INV:
                T_BL_top = self.temp(self.z_BL)

                return T_BL_top + (z - self.z_BL)*self.dTdz_2
            else:
                T_top = self.temp(self.z_INV)
                return T_top

        return f(z)

    def p(self, z):
        @np.vectorize
        def f(z):
            if z == 0.0:
                return self.p0
            else:
                return scipy.interp(z, self._z, self._p)
        return f(z)


        # @np.vectorize
        # def f(z):
            # if z == 0.:
                # return self.p0
            # elif z <= self.z_BL:
                # return self.BL_profile.p(z)
            # elif z <= self.z_INV:
                # return self.MIDDLE_profile.p(z - self.z_BL)
            # else:
                # return self.TOP_profile.p(z - self.z_INV)
        # return f(z)

    def rel_humidity(self, z):
        @np.vectorize
        def f(z):
            if z <= self.z_BL:
                qv_sat = self.qv_sat(z)

                return self.q_v0/qv_sat
            elif z <= self.z_INV:
                RH_BL_top = self.rel_humidity(self.z_BL)

                return RH_BL_top + (z - self.z_BL)*self.dRHdz_1
            else:
                RH_top = self.rel_humidity(self.z_INV)

                return RH_top + (z - self.z_INV)*self.dRHdz_2

        return f(z)

    def q_v(self, z):
        return self.rel_humidity(z)*self.qv_sat(z)


    def rho(self, z):
        @np.vectorize
        def f(z):
            if z == 0.0:
                return self.rho0
            else:
                return scipy.interp(z, self._z, self._rho)
        return f(z)

        # from pyclouds import parameterisations
        # T = self.temp(z)
        # p = self.p(z)
        # qv_sat__f = parameterisations.ParametersationsWithSpecificConstants(self.constants).pv_sat.qv_sat
        # qv = self.rel_humidity(z)*qv_sat__f(T=T, p=p)

        # R_v = self.constants.get('R_v')
        # R_d = self.constants.get('R_d')

        # qd = 1.0 - qv

        # rho_inv = (qd*R_d + qv*R_v)*T/p
        # return 1.0/rho_inv

        # @np.vectorize
        # def f(z):
            # if z < self.z_BL:
                # return self.BL_profile.rho(z)
            # else:
                # return self.MIDDLE_profile.rho(z - self.z_BL)
        # return f(z)

    def __str__(self):
        return "Idealised Two-layer atmosphere ($RH_0={RH0:.0f}\%$, "\
               "$T_0={T0}K$, $z_{{BL}}={z_bl}m$, $z_{{INV}}={z_INV}m$)".format(
               RH0=self.RH0*100., T0=self.T0, z_bl=self.z_BL, z_INV=self.z_INV)


if __name__ == "__main__":
    from matplotlib import pyplot as plot

    profile = Soong1973()

    plot.ion()
    z = np.linspace(0.0, 20000.0, 100)
    temp = profile.temp([z])

    plot.subplot(131)
    plot.plot(profile.temp([z]), z)
    #plot.plot(profile.dew_point([z]), z)
    plot.xlabel("Temperature [K]")
    plot.ylabel("Height [m]")
    plot.grid(True)

    plot.subplot(132)
    plot.plot(profile.rel_humidity([z]), z)
    plot.xlabel("Relative humidity [%]")
    plot.ylabel("Height [m]")
    plot.xlim(0.0, 1.0)
    plot.grid(True)

    plot.subplot(133)
    plot.plot(profile.p([z]), z)
    plot.xlabel("Pressure [Pa]")
    plot.ylabel("Height [m]")
    plot.xlim(0.0, None)
    plot.grid(True)

    plot.suptitle("Atmospheric stratification profile from Soong 1973")

    plot.draw()
    raw_input()


    r = RICO()

    z = np.linspace(0., 4e3, 100)

    print r.p(z)

    import matplotlib.pyplot as plot
    plot.plot(r.p(z), z)
    plot.draw()
    plot.show()
