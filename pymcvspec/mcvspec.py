"""Python implementation of MCVSPEC.

This module contains a user-friendly interface to the C++ implementation of
MCVSPEC (using pybind11). It uses astropy's units package for input and
provides spectral output using pyatomdb. Currently there is no reflection
model implemented directly in python.
"""

import numpy as np
import pyatomdb
import astropy.units as u
from astropy.constants import G

from _pymcvspec import _cataclysmic_variable, _dipole, _white_dwarf
from _pymcvspec import _mass_to_radius, _luminosity_to_mdot
from _pymcvspec import _atomic_charges, _atomic_masses

cgs = [(u.statC, ((u.g*u.cm**3)**0.5)/u.s, lambda x: x, lambda x: x),
       (u.G, ((u.g/u.cm)**0.5/u.s), lambda x: x, lambda x: x)]

atomic_charges = np.array(_atomic_charges)
atomic_masses = np.array(_atomic_masses)*u.u

@u.quantity_input
def mass_to_radius(mass: u.Quantity[u.M_sun]) -> u.Quantity[u.cm]:
    """Finds the WD radius corresponding to an input mass.

    The white dwarf mass-radius relationship is specified by a header file
    which is supplied at compile-time. The header contains two C-arrays
    which contain samples of the mass-radius relationship. Linear
    interpolation is used for arbitrary sampling.

    Parameters
    ----------
    mass : `~astropy.units.Quantity`
        White dwarf mass

    Returns
    -------
    radius : `~astropy.units.Quantity`
        Corresponding white dwarf radius (in cm)

    References
    ----------
    ..[1] E. E. Salpeter, "Energy and Pressure of a Zero-Temperature
          Plasma", Astrophysical Journal, vol. 134, p.669, 1961
          (https://ui.adsabs.harvard.edu/abs/1961ApJ...134..669S)
    ..[2] T. Hamada, and E. E. Salpeter, "Models for Zero-Temperature
          Stars", Astrophysical Journal, vol. 134, p.683, 1961
          (https://ui.adsabs.harvard.edu/abs/1961ApJ...134..683H)
    """
    if not mass.unit.is_equivalent(u.g):
        raise u.UnitTypeError("mass must have units of mass")
    return _mass_to_radius(mass.to_value(u.g))*u.cm

@u.quantity_input
def luminosity_to_mdot(luminosity : u.Quantity[u.erg/u.s],
                       mass : u.Quantity[u.M_sun],
                       radius : u.Quantity[u.R_sun],
                       mag_radius : u.Quantity[u.R_sun] = 0*u.cm
                      ) -> u.Quantity[u.g/u.s]:
    """Computes the accretion rate of an mCV based on its luminosity.

    Parameters
    ----------
    luminosity : `~astropy.units.Quantity`
        Bolometric luminosity of the mCV

    mass : `~astropy.units.Quantity`
        White dwarf mass

    radius : `~astropy.units.Quantity`
        White dwarf radius

    mag_radius: `~astropy.units.Quantity`
        Magnetospheric radius

    Returns
    -------
    mdot : `~astropy.units.Quantity`
        Accretion rate in g/s

    Notes
    -----
    The accretion rate is computed by equating the bolometric luminosity
    of the mCV to the gravitational energy avialable is the source :
    L = G*M*mdot*(1/R_wd - 1/R_m). Where mdot is the mass accretion rate
    (mass/time)
    """

    if not luminosity.unit.is_equivalent(u.erg/u.s):
        raise u.UnitTypeError("luminosity must have units of energy/time")
    if not mass.unit.is_equivalent(u.g):
        raise u.UnitTypeError("mass must have units of mass")
    if not radius.unit.is_equivalent(u.cm):
        raise u.UnitTypeError("radius must have units of length")
    if not mag_radius.unit.is_equivalent(u.cm):
        raise u.UnitTypeError("mag_radius must have units of length")
    irm = 0
    if mag_radius != 0:
        irm = 1/mag_radius.to_value(u.cm)
    lum = luminosity.to_value(u.erg/u.s)
    m = mass.to_value(u.g)
    r = radius.to_value(u.cm)
    return _luminosity_to_mdot(lum, _white_dwarf(m, r, inv_mag_rad=irm))*u.g/u.s

class dipole(_dipole):
    """Public interface to _dipole.

    Store information about the geometry of the accretion column and
    provides a method to compute relevant properties of the column
    coordiante system. As the name suggests, this class specifically
    defines a dipolar coordinate system. It considers coordinate
    variability only along magnetic field lines.

    Parameters
    ----------
    u : float
        `u` coordinate of the column, the "choice" of field line.
        u=sin(theta)**2 where theta is the magnetic colatitude of the
        column at the WD surface

    Attributes
    __________
    w_0 : float
        `w` coordinate of the column at the WD surface. w=sin(theta) with
        theta defined as above
    a_0 : float
        cross-sectional area proxy at the WD surface. With h_i as the
        metric scaling factors a_0 = h_u*h_phi with h evaluated at w=w_0.
        The product h_u*h_phi is proportional to the column's
        cross-sectional area.

    References
    ----------
    ..[1] J. B. G. Canalle, "Accretion in dipole magnetic fields: flow
        structure and X-ray emission of accreting white dwarfs", Astronomy
        and Astrophysics, Volume 440, Issue 1, September II 2005,
        pp.185-198 (https://ui.adsabs.harvard.edu/abs/2005A&A...440..185C)
    ..[2] H. Wang, "A General Curvilinear Magnetic Field-Line-Following
        Coordinate System for Ionosphere-Plasmasphere Modeling",
        Journal of Geophysical Research: Space Physics, Volume 127,
        Issue 3, article id. e30017
        (https://ui.adsabs.harvard.edu/abs/2022JGRA..12730017W)
    """

    def __init__(self, u : float) -> None:
        super().__init__(u)
    """
    update_coordinates(w)

    computes geometric terms for an input coordinate "w"

    Parameters
    ----------
    w : float
        `w` coordinate

    Returns
    -------
    r : float
        corresponding radial coordinate
    proj_r_w : float
        projection of `w` onto `r`
    convergance : float
        (1/a)*da/dw where a = h_u*h_phi
    scale_factors : `~numpy.ndarray`
        [h_u, h_w, h_phi] where h_i are the scale factors assocaited with
        the subscripted coordinates
    """


class cataclysmic_variable:
    """Base class for all mCV objects.

    The `cataclysmic_variable` class is a base class that can be constructed
    from, often unobservable, physical parameters. It calls an underlying
    C++ object to solve for the thermal profile of the mCV accretion column
    and exposes the profile as a selection of attributes.

    Parameters
    ----------
    mass : `~astropy.units.Quantity`
        white dwarf mass
    b_field : `~astropy.units.Quantity`
        surface magnetic field strength at base of accretion column (in G,
        using gaussian cgs units)
    mdot : `~astropy.units.Quantity`
        specific accretion rate (in g/cm2/s)
    accretion_area : `~astropy.units.Quantity`
        area of accretion column at WD surface
    magnetospheric_radius : `~astropy.units.Quantity`
        magnetospheric radius of WD
    corotation_radius : `~astropy.units.Quantity`
        corotation radius, defined as cbrt(G*M*P_spin**2/(4*pi**2))
    metallicity : float
        metallicity of accretion column relative to solar abundance
    shock_ratio : float
        ratio of electron and ion partial pressures at the shock front
    orbital_inclination : `~astropy.units.Quantity`
        inclination angle of WD to observer (only relevant for reflection)
    column_magnetic_colatitude : `~astropy.units.Quantity`
        magnetic colatitude of the column
    distance : `~astropy.units.Quantity`
        distance to the source

    Attributes
    __________
    accretion_rate : `~astropy.units.Quantity`
        total accretion rate (in g/s)
    geometry : `~pymcvspec.dipole`
        class which defines the geometry of the column
    shock_height : `~astropy.units.Quantity`
        radial altitude of the shock (in cm)
    mbar : `~astropy.units.Quantity`
        average mass of ions in the accretion column (in g)
    zbar : float
        average atomic charge of ions in the accretion column
    density_const : float
        constant which defines the relationship between mass density and
        electron number density. n_e = density_const*density/mbar and
        density_const = zbar/(1 + zbar*m_e/mbar)
    exchange_const : float
        prefactor used in computing electron-ion energy exchange rate in
        dimensionless unit system
    bremss_const : float
        prefactor used in computing bremsstrahlug radiation rate in
        dimensionless unit system
    cyclotron_const : float
        prefactor used in computing cyclotron to bremsstrahlug ratio in
        dimensionless unit system
    length_unit : `~astropy.units.Quantity`
        conversion factor between length in dimensionless unit system and
        cgs
    mass_unit : `~astropy.units.Quantity`
        conversion factor between mass in dimensionless unit system and
        cgs
    time_unit : `~astropy.units.Quantity`
        conversion factor between time in dimensionless unit system and
        cgs
    velocity_unit : `~astropy.units.Quantity`
        conversion factor between velocity in dimensionless unit system
        and cgs
    volume_unit : `~astropy.units.Quantity`
        conversion factor between volume in dimensionless unit system
        and cgs
    energy_unit : `~astropy.units.Quantity`
        conversion factor between energy in dimensionless unit system
        and cgs
    density_unit : `~astropy.units.Quantity`
        conversion factor between density in dimensionless unit system
        and cgs
    altitude : `~astropy.units.Quantity`
        altitude of plasma element in thermal profile
    volume : `~astropy.units.Quantity`
        volume of plasma element in thermal profile
    velocity : `~astropy.units.Quantity`
        bulk flow velocity of plasma element in thermal profile
    density : `~astropy.units.Quantity`
        mass density of plasma element in thermal profile
    total_pressure : `~astropy.units.Quantity`
        total pressure of plasma element in thermal profile
    electron_pressure : `~astropy.units.Quantity`
        electron partial pressure of plasma element in thermal profile
    electron_density : `~astropy.units.Quantity`
        electron number density of plasma element in thermal profile
    electron_temperature : `~astropy.units.Quantity`
        electron temperature (in keV) of plasma element in thermal profile
    ion_pressure : `~astropy.units.Quantity`
        ion partial pressure of plasma element in thermal profile
    ion_density : `~astropy.units.Quantity`
        ion number density of plasma element in thermal profile
    ion_temperature : `~astropy.units.Quantity`
        ion temperature (in keV) of plasma element in thermal profile
    """

    @u.quantity_input
    def __init__(
        self,
        mass: u.Quantity[u.M_sun],
        b_field: u.Quantity[u.MG],
        mdot: u.Quantity[u.g/u.cm**2/u.s],
        accretion_area: u.Quantity[u.cm**2],
        magnetospheric_radius: u.Quantity[u.cm] = 0*u.cm,
        corotation_radius: u.Quantity[u.cm] = 0*u.cm,
        metallicity: u.Quantity[u.dimensionless_unscaled] = 1,
        shock_ratio: u.Quantity[u.dimensionless_unscaled] = 0.75,
        orbital_inclination: u.Quantity[u.deg] = 45*u.deg,
        column_magnetic_colatitude: u.Quantity[u.deg] = 1*u.deg,
        distance: u.Quantity[u.pc] = 1*u.pc,
    ) -> None:
        self.mass = mass.to(u.M_sun)
        self.radius = mass_to_radius(mass).to(u.R_sun)
        self.b_field = b_field.to(u.MG)
        self.mdot = mdot.to(u.g/u.cm**2/u.s)
        self.accretion_rate = (mdot*accretion_area).to(u.g/u.s)
        self.accretion_area = accretion_area.to(u.cm**2)
        self.magnetospheric_radius = magnetospheric_radius.to(u.cm)
        self.corotation_radius = corotation_radius.to(u.cm)
        self.metallicity = metallicity
        self.shock_ratio = shock_ratio
        self.orbital_inclination = orbital_inclination.to(u.deg)
        self.magnetic_colatitude = column_magnetic_colatitude.to(u.deg)
        self.distance = distance.to(u.pc)
        if magnetospheric_radius.value == 0:
            magnetospheric_radius = np.inf*u.cm
            corotation_radius = 1*u.cm
        irm = 1/magnetospheric_radius
        cos_incl = np.cos(self.orbital_inclination.to_value(u.radian))
        u_coord = np.sin(self.magnetic_colatitude.to_value(u.radian))**2
        self._cpp_impl = _cataclysmic_variable(
            mass=self.mass.to_value(u.g),
            radius=self.radius.to_value(u.cm),
            b_field=self.b_field.to_value(u.G),
            mdot=self.mdot.to_value(u.g/u.cm**2/u.s),
            area=accretion_area.to_value(u.cm**2),
            inv_r_m=irm.to_value(1/u.cm),
            corot_radius=corotation_radius.to_value(u.cm),
            metallicity=metallicity,
            cos_incl_angle=cos_incl,
            shock_ratio=self.shock_ratio,
            column_coord=u_coord,
            src_distance=distance.to_value(u.cm)
        )
        self.geometry = dipole(self._cpp_impl.column_coord)
        status = self._cpp_impl.solve()
        if status == -1:
            raise RuntimeError("No valid solution found: Column does not reach WD surface for any shock height")

    @property
    def shock_height(self):
        return self._cpp_impl.shock_height*u.cm
    @property
    def abundance(self):
        return self._cpp_impl.abundance
    @property
    def mbar(self):
        return self._cpp_impl.average_ion_mass*u.g
    @property
    def zbar(self):
        return self._cpp_impl.average_ion_charge
    @property
    def density_const(self):
        return self._cpp_impl.density_const
    @property
    def exchange_const(self):
        return self._cpp_impl.exchange_const
    @property
    def bremss_const(self):
        return self._cpp_impl.bremss_const
    @property
    def cyclotron_const(self):
        return self._cpp_impl.cyclotron_const
    @property
    def length_unit(self):
        return self._cpp_impl.length_converter*u.cm
    @property
    def mass_unit(self):
        return self._cpp_impl.mass_converter*u.g
    @property
    def time_unit(self):
        return self._cpp_impl.time_converter*u.s
    @property
    def velocity_unit(self):
        return self._cpp_impl.velocity_converter*u.cm/u.s
    @property
    def volume_unit(self):
        return self._cpp_impl.volume_converter*u.cm**3
    @property
    def energy_unit(self):
        return self._cpp_impl.energy_converter*u.erg
    @property
    def density_unity(self):
        return self._cpp_impl.density_converter*u.g/u.cm**3
    @property
    def altitude(self):
        return self._cpp_impl.altitude*u.cm
    @property
    def volume(self):
        return self._cpp_impl.volume*u.cm**3
    @property
    def velocity(self):
        return self._cpp_impl.velocity*u.cm/u.s
    @property
    def density(self):
        return self._cpp_impl.density*u.g/u.cm**3
    @property
    def total_pressure(self):
        return self._cpp_impl.total_pressure*u.erg/u.cm**3
    @property
    def electron_pressure(self):
        return self._cpp_impl.electron_pressure*u.erg/u.cm**3
    @property
    def electron_density(self):
        return self._cpp_impl.electron_density/u.cm**3
    @property
    def electron_temperature(self):
        return self._cpp_impl.electron_temperature*u.keV
    @property
    def ion_temperature(self):
        return self._cpp_impl.ion_temperature*u.keV
    @property
    def ion_density(self):
        return self.electron_density/self.zbar
    @property
    def ion_pressure(self):
        return self.total_pressure-self.electron_pressure

    @u.quantity_input
    def spectrum(self,
                 energy_bins: u.Quantity[u.keV]
                 ) -> u.Quantity[1/u.s/u.keV/u.cm**2]:
        """Produces spectrum from the thermal profile.

        Uses `pyatomdb` to compute the expected optically-thin
        bremsstrahlung spectrum from the mCV thermal profile

        Parameters
        ----------
        energy_bins : `~astropy.units.Quantity`
            bounds of energy bins for the output spectrum

        Returns
        -------
        flux : `~astropy.units.Quantity`
            photon flux in each bin (units of photons/s/cm**2/keV)
        """

        session = pyatomdb.spectrum.CIESession()
        session.set_abund(atomic_charges, self.metallicity)
        session.set_response(energy_bins.to_value(u.keV), raw=True)
        flux = np.zeros(len(energy_bins)-1)/(u.s*u.keV*u.cm**2)
        apec_unit = ((u.cm**3)/u.s/energy_bins.unit)
        for kT, n_e, n_i, vol in zip(
            self.electron_temperature,
            self.electron_density,
            self.electron_density/self.zbar,
            self.volume,
        ):
            # The units of a pyatomdb "spectrum" are photons*cm**3/s/keV
            # This is normalized to a flux (units of photons/s/cm**2/keV)
            # by multiplying by the plasma "emission measure" divided by
            # surface area of a sphere with radius "d". If an arf is set
            # pyatomdb will instead compute a spectrum with units
            # photons*cm**5/s/keV and the arf would need to be factored out.
            norm = n_e*n_i*vol/(4*np.pi*self.distance**2)
            flux += session.return_spectrum(kT.to_value(u.keV))*apec_unit*norm
        return flux.to(1/u.s/u.keV/u.cm**2)


class polar(cataclysmic_variable):
    """Subclass of `cataclysmic_variable` optimized for polars.

    This class uses "observational" inputs to infer the values of the
    physical parameters required by the `cataclysmic_variable` class.
    It assumes that the magnetospheric radius -> infinity due to the
    strong B-field found in polars.

    Parameters
    ----------
    mass : `~astropy.units.Quantity`
        white dwarf mass
    b_field : `~astropy.units.Quantity`
        surface magnetic field strength at base of accretion column (in G,
        using gaussian cgs units)
    luminosity : `~astropy.units.Quantity`
        bolometric luminosity of mCV
    accretion_area : `~astropy.units.Quantity`
        area of accretion column at WD surface
    metallicity : float
        metallicity of accretion column relative to solar abundance
    shock_ratio : float
        ratio of electron and ion partial pressures at the shock front
    orbital_inclination : `~astropy.units.Quantity`
        inclination angle of WD to observer (only relevant for reflection)
    distance : `~astropy.units.Quantity`
        distance to the source

    Other Parameters
    ----------------
    fractional_area : float
        if accretion_area is 0 cm `polar` uses this instead to calculate
        the accretion area based on the surface area of the WD.

    See Also
    --------
    cataclysmic_variable
    """

    @u.quantity_input
    def __init__(
        self,
        mass: u.Quantity[u.M_sun],
        b_field: u.Quantity[u.MG],
        luminosity: u.Quantity[u.erg/u.s],
        accretion_area: u.Quantity[u.cm**2] = 0*u.cm**2,
        fractional_area: u.Quantity[u.dimensionless_unscaled] = 1e-3,
        metallicity: u.Quantity[u.dimensionless_unscaled] = 1,
        shock_ratio: u.Quantity[u.dimensionless_unscaled] = 0.75,
        orbital_inclination: u.Quantity[u.deg] = 45*u.deg,
        distance: u.Quantity[u.pc] = 1*u.pc,
    ) -> None:
        radius = mass_to_radius(mass)
        if accretion_area == 0*u.cm**2:
            accretion_area = fractional_area*4*np.pi*(radius**2)
        mdot = luminosity_to_mdot(luminosity, mass, radius)/accretion_area
        cataclysmic_variable.__init__(
            self,
            mass,
            b_field,
            mdot,
            accretion_area,
            metallicity=metallicity,
            shock_ratio=shock_ratio,
            orbital_inclination=orbital_inclination,
            distance=distance,
        )

class intermediate_polar(cataclysmic_variable):
    """Subclass of `cataclysmic_variable` optimized for IPs.

    This class uses "observational" inputs to infer the values of the
    physical parameters required by the `cataclysmic_variable` class.
    It derives the magnetic field strength of the mCV based on the
    magnetospheric radius.

    Parameters
    ----------
    mass : `~astropy.units.Quantity`
        white dwarf mass
    spin_period : `~astropy.units.Quantity`
        spin period of the WD
    luminosity : `~astropy.units.Quantity`
        bolometric luminosity of mCV
    accretion_area : `~astropy.units.Quantity`
        area of accretion column at WD surface
    metallicity : float
        metallicity of accretion column relative to solar abundance
    shock_ratio : float
        ratio of electron and ion partial pressures at the shock front
    orbital_inclination : `~astropy.units.Quantity`
        inclination angle of WD to observer (only relevant for reflection)
    distance : `~astropy.units.Quantity`
        distance to the source
    mag_radius_ratio : float
        ratio of the magnetospheric radius to the corotation radius

    Other Parameters
    ----------------
    fractional_area : float
        if accretion_area is 0 cm `polar` uses this instead to calculate
        the accretion area based on the surface area of the WD.

    See Also
    --------
    cataclysmic_variable
    """

    @u.quantity_input
    def __init__(
        self,
        mass: u.Quantity[u.M_sun],
        spin_period: u.Quantity[u.s],
        luminosity: u.Quantity[u.erg/u.s],
        accretion_area: u.Quantity[u.cm**2] = 0*u.cm**2,
        fractional_area: u.Quantity[u.dimensionless_unscaled] = 1e-3,
        metallicity: u.Quantity[u.dimensionless_unscaled] = 1,
        shock_ratio: u.Quantity[u.dimensionless_unscaled] = 0.75,
        orbital_inclination: u.Quantity[u.deg] = 45*u.deg,
        distance: u.Quantity[u.pc] = 1*u.pc,
        mag_radius_ratio: float = 1
    ) -> None:
        radius = mass_to_radius(mass)
        corotation_radius = np.cbrt(G*mass*(spin_period**2)/(4*np.pi*np.pi))
        r_m = mag_radius_ratio*corotation_radius
        if accretion_area == 0*u.cm**2:
            accretion_area = fractional_area*4*np.pi*(radius**2)
        mdot = luminosity_to_mdot(luminosity, mass, radius, r_m)/accretion_area
        b_field = (np.sqrt(32*mdot*accretion_area*np.sqrt(G*mass*(r_m**7)))/(radius**3))
        b_field = b_field.to(u.G, equivalencies=cgs)
        cataclysmic_variable.__init__(
            self,
            mass,
            b_field,
            mdot,
            accretion_area,
            magnetospheric_radius=r_m,
            corotation_radius=corotation_radius,
            metallicity=metallicity,
            shock_ratio=shock_ratio,
            orbital_inclination=orbital_inclination,
            distance=distance,
        )
