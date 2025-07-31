from _pymcvspec import _cataclysmic_variable, _mass_to_radius, _luminosity_to_mdot
import numpy as np
import astropy.units as u
from astropy.constants import G

cgs = [(u.statC, ((u.g*u.cm**3)**0.5)/u.s, lambda x:x, lambda x:x),
       (u.G, ((u.g/u.cm)**0.5/u.s), lambda x:x, lambda x:x)]

class cataclysmic_variable:
    @u.quantity_input
    def __init__(
        self,
        mass: u.Quantity[u.M_sun],
        b_field: u.Quantity[u.MG],
        accretion_rate: u.Quantity[u.g/u.s],
        accretion_area: u.Quantity[u.cm**2],
        magnetospheric_radius: u.Quantity[u.cm] = 0*u.cm,
        metalicity=1,
        shock_ratio=0.75,
        area_exponent=0,
        cos_inclination_angle=0.5,
        distance: u.Quantity[u.pc] = 1*u.pc
    ) -> None:
        self.mass = mass.to(u.M_sun)
        self.radius = (_mass_to_radius(mass.to_value(u.g))*u.cm).to(u.R_sun)
        self.b_field = b_field.to(u.MG)
        self.accretion_rate = accretion_rate.to(u.g/u.s)
        self.accretion_area = accretion_area.to(u.cm**2)
        self.mdot = (accretion_rate/accretion_area).to(u.g/u.cm**2/u.s)
        self.magnetospheric_radius = magnetospheric_radius.to(u.cm)
        self.metalicity = 1
        self.shock_ratio = shock_ratio
        self.n = area_exponent
        self.cos_incl = cos_inclination_angle
        self.distance = distance.to(u.pc)
        if magnetospheric_radius==0:
            irm = 0./u.cm
        else:
            irm = 1/magnetospheric_radius
        self.cpp_impl = _cataclysmic_variable(mass.to_value(u.g), self.radius.to_value(u.cm), b_field.to_value(u.G), accretion_rate.to_value(u.g/u.s),
                                              irm.to_value(1/u.cm), metalicity, accretion_area.to_value(u.cm**2),
                                              cos_inclination_angle, area_exponent, distance.to_value(u.cm), 1)

        self.altitude = self.cpp_impl.get_altitude()*u.cm
        self.electron_temperature = self.cpp_impl.get_electron_temperature()*u.keV
        self.ion_temperature = self.cpp_impl.get_ion_temperature()*u.keV
        self.electron_density = self.cpp_impl.get_electron_density()/u.cm**3
        self.ion_density = self.cpp_impl.get_ion_density()/u.cm**3
        self.total_pressure = self.cpp_impl.get_total_pressure()*u.dyne/u.cm**2
        self.electron_pressure = self.cpp_impl.get_electron_pressure()*u.dyne/u.cm**2

class polar(cataclysmic_variable):
    @u.quantity_input
    def __init__(
        self,
        mass: u.Quantity[u.M_sun],
        b_field: u.Quantity[u.MG],
        luminosity: u.Quantity[u.erg/u.s],
        accretion_area: u.Quantity[u.cm**2]=None,
        fractional_area=1e-3,
        metalicity=1,
        shock_ratio=0.75,
        area_exponent=0,
        cos_inclination_angle=0.5,
        distance: u.Quantity[u.pc] = 1*u.pc
    ) -> None:
        radius = (_mass_to_radius(mass.to_value(u.g))*u.cm).to(u.R_sun)
        if accretion_area is None:
            accretion_area = fractional_area*4*np.pi*(radius**2)
        mdot = _luminosity_to_mdot(luminosity.to_value(u.erg/u.s), mass.to_value(u.g), radius.to_value(u.cm), 0)
        cataclysmic_variable.__init__(self, mass, b_field, mdot, accretion_area, 0*u.cm, metalicity, shock_ratio, area_exponent, cos_inclination_angle, distance)

class intermediate_polar(cataclysmic_variable):
    @u.quantity_input
    def __init__(
        self,
        mass: u.Quantity[u.M_sun],
        spin_period: u.Quantity[u.s],
        luminosity: u.Quantity[u.erg/u.s],
        accretion_area: u.Quantity[u.cm**2]=None,
        fractional_area=1e-3,
        metalicity=1,
        shock_ratio=0.75,
        area_exponent=0,
        cos_inclination_angle=0.5,
        distance: u.Quantity[u.pc] = 1*u.pc,
        mag_radius_ratio=1
    ) -> None:
        radius = (_mass_to_radius(mass.to_value(u.g))*u.cm).to(u.R_sun)
        r_m = mag_radius_ratio*np.cbrt(G*mass*(spin_period**2)/(4*np.pi*np.pi))
        if accretion_area is None:
            accretion_area = fractional_area*4*np.pi*(radius**2)
        mdot = _luminosity_to_mdot(luminosity.to_value(u.erg/u.s), mass.to_value(u.g), radius.to_value(u.cm), 1/r_m.to_value(u.cm))*u.g/u.s
        b_field = (np.sqrt(32*mdot*np.sqrt(G*mass*(r_m**7)))/(radius**3)).to(u.G, equivalencies=cgs)
        cataclysmic_variable.__init__(self, mass, b_field, mdot, accretion_area, r_m, metalicity, shock_ratio, area_exponent, cos_inclination_angle, distance)
