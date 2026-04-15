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

from _pymcvspec import _cataclysmic_variable, _saxton_cv, _cropper_cv, _wu_cv
from _pymcvspec import _dipole, _white_dwarf, _accretion_column, _tolerance
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

class white_dwarf:
    """Public interface to _white_dwarf

    container for storing information about the white dwarf in an mCV

    Parameters
    __________
    mass : `~astropy.units.Quantity`
        white dwarf mass
    b_field : `~astropy.units.Quantity`
        surface magnetic field strength at base of accretion column (in G,
        using gaussian cgs units)
    orbital_inclination : `~astropy.units.Quantity`
        inclination angle of WD to observer (only relevant for reflection)
    magnetospheric_radius : `~astropy.units.Quantity`
        magnetospheric radius of WD
    corotation_radius : `~astropy.units.Quantity`
        corotation radius, defined as cbrt(G*M*P_spin**2/(4*pi**2))
    distance : `~astropy.units.Quantity`
        distance to the source
    """

    @u.quantity_input
    def __init__(
        self,
        mass: u.Quantity[u.M_sun],
        b_field: u.Quantity[u.MG],
        orbital_inclination: u.Quantity[u.deg] = 45*u.deg,
        magnetospheric_radius: u.Quantity[u.cm] = 0*u.cm,
        corotation_radius: u.Quantity[u.cm] = 1*u.cm,
        distance: u.Quantity[u.pc] = 1*u.pc,
    ) -> None:
        radius = mass_to_radius(mass)
        self._orbital_inclination = orbital_inclination
        self._magnetospheric_radius = magnetospheric_radius

        cosi = np.cos(orbital_inclination.to_value(u.radian))
        if magnetospheric_radius.value == 0:
            magnetospheric_radius = np.inf*u.cm
            corotation_radius = 1*u.cm
        irm = 1/magnetospheric_radius
        self._cpp_impl = _white_dwarf(
            mass=mass.to_value(u.g),
            radius=radius.to_value(u.cm),
            b_field=b_field.to_value(u.gauss),
            cos_incl=cosi,
            inv_mag_rad=irm.to_value(1./u.cm),
            corot_rad=corotation_radius.to_value(u.cm),
            distance=distance.to_value(u.cm)
        )

    @property
    def mass(self):
        """get or set mass

        Setting the WD mass will automatically adjust the WD radius

        """
        return self._cpp_impl.mass*u.g
    @mass.setter
    @u.quantity_input
    def mass(self, m: u.Quantity[u.M_sun]):
        self._cpp_impl.mass = m.to_value(u.g)
        self._cpp_impl.radius = mass_to_radius(m).to_value(u.cm)

    @property
    def radius(self):
        """get or set radius

        Setting the WD radius directly will not adjust the WD mass. The radius
        will be automatically adjusted by setting the mass. This allows the
        user to play around with different WD mass/radius relations.

        """
        return self._cpp_impl.radius*u.cm
    @radius.setter
    @u.quantity_input
    def radius(self, r: u.Quantity[u.R_sun]):
        self._cpp_impl.radius = r.to_value(u.cm)

    @property
    def b_field(self):
        return self._cpp_impl.b_field*u.gauss
    @b_field.setter
    @u.quantity_input
    def b_field(self, b_field: u.Quantity[u.MG]):
        self._cpp_impl.b_field = b_field.to_value(u.gauss)

    @property
    def orbital_inclination(self):
        """get or set orbital inclination
        This method will also compute cos(i) to update the c++ implementation

        """
        return self._orbital_inclination
    @orbital_inclination.setter
    @u.quantity_input
    def orbital_inclination(self, i: u.Quantity[u.deg]):
        self._orbital_inclination = i
        self._cpp_impl.cos_inclination = np.cos(i.to_value(u.radian))

    @property
    def magnetospheric_radius(self):
        """get or set the magnetoshperic radius
        This method will also compute 1/r_m which is stored in the c++
        implementation. It does not have bounds checking.

        """
        return self._magnetospheric_radius
    @magnetospheric_radius.setter
    @u.quantity_input
    def magnetospheric_radius(self, r_m: u.Quantity[u.cm]):
        self._magnetospheric_radius = r_m
        self._cpp_impl.inverse_mag_radius = (1./r_m).to_value(1./u.cm)

    @property
    def corotation_radius(self):
        return self._cpp_impl.corotation_radius*u.cm
    @corotation_radius.setter
    @u.quantity_input
    def corotation_radius(self, r_co: u.Quantity[u.cm]):
        self._cpp_impl.corotation_radius = r_co.to_value(u.cm)

    @property
    def distance(self):
        return self._cpp_impl.distance*u.cm
    @distance.setter
    @u.quantity_input
    def distance(self, d: u.Quantity[u.pc]):
        self._cpp_impl.distance = d.to_value(u.cm)

class accretion_column:
    """Public interface to _accretion_column

    container for storing information about the accretion column in an mCV

    Parameters
    __________
    mdot : `~astropy.units.Quantity`
        specific accretion rate (g/cm2/s)
    area : `~astropy.units.Quantity`
        area of the base of the accretion column
    metallicity : float
        metallicity (relative to solar) of the material in the column
    shock_ratio : float
        ratio of electron to ion partial pressure at the shock
    colatitude : `~astropy.units.Quantity`
        magnetic colatitude of the base of the accretion column
    """

    @u.quantity_input
    def __init__(
        self,
        mdot: u.Quantity[u.g/u.s/u.cm**2],
        area: u.Quantity[u.cm**2],
        metallicity: u.Quantity[u.dimensionless_unscaled] = 1,
        shock_ratio: u.Quantity[u.dimensionless_unscaled] = 0.75,
        colatitude: u.Quantity[u.deg] = 0*u.deg
    ) -> None:
        self._colatitude = colatitude
        self._cpp_impl = _accretion_column(
            mdot=mdot.to_value(u.g/u.s/u.cm**2),
            area=area.to_value(u.cm**2),
            metallicity=metallicity,
            shock_ratio=shock_ratio,
            sin_colat=np.sin(colatitude.to_value(u.radian))
        )

    @property
    def mdot(self):
        return self._cpp_impl.mdot*u.g/u.s/u.cm**2
    @mdot.setter
    @u.quantity_input
    def mdot(self, md: u.Quantity[u.g/u.s/u.cm**2]):
        self._cpp_impl.mdot = md.to_value(u.g/u.s/u.cm**2)

    @property
    def area(self):
        return self._cpp_impl.area*u.cm**2
    @area.setter
    @u.quantity_input
    def area(self, a: u.Quantity[u.cm**2]):
        self._cpp_impl.area = a.to_value(u.cm**2)

    @property
    def metallicity(self):
        return self._cpp_impl.metallicity
    @metallicity.setter
    @u.quantity_input
    def metallicity(self, abund: u.Quantity[u.dimensionless_unscaled]):
        self._cpp_impl.metallicity = abund

    @property
    def shock_ratio(self):
        return self._cpp_impl.shock_ratio
    @shock_ratio.setter
    @u.quantity_input
    def shock_ratio(self, r: u.Quantity[u.dimensionless_unscaled]):
        self._cpp_impl.shock_ratio = r

    @property
    def magnetic_colatitude(self):
        """get or set the magnetic colatitiude of the column base
        This method will also compute cos(b) to update the c++ implementation

        """
        return self._colatitude
    @magnetic_colatitude.setter
    @u.quantity_input
    def magnetic_colatitude(self, b: u.Quantity[u.deg]):
        self._colatitude = b
        self._cpp_impl.sin_colat=np.sin(b.to_value(u.radian))

class tolerance:
    """Public interface to _tolerance

    container for storing numerical methods parameters
    Parameters
    __________
    abs_err : float
        absolute tolerance for numerical integration and root finding
    rel_err : float
        relative tolerance for numerical integration
    delta_kT : `~astropy.units.Quantity`
        minimum grid spacing in temperature
    delta_z : float
        minimum grid spacing in altitude (fractional)
    """

    @u.quantity_input
    def __init__(
        self,
        abs_err: float = 1e-8,
        rel_err: float = 1e-6,
        delta_kT: u.Quantity[u.keV] = 0.5*u.keV,
        delta_z: u.Quantity[u.dimensionless_unscaled] = 0.1,
    ) -> None:
        self._cpp_impl = _tolerance(abserr=abs_err,
                                    relerr=rel_err,
                                    dkT=delta_kT.to_value(u.keV),
                                    dz=delta_z)

    @property
    def abs_err(self):
        return self._cpp_impl.abs_err
    @abs_err.setter
    def abs_err(self, abserr: float):
        self._cpp_impl.abs_err = abserr

    @property
    def rel_err(self):
        return self._cpp_impl.rel_err
    @rel_err.setter
    def rel_err(self, relerr: float):
        self._cpp_impl.rel_err = relerr

    @property
    def delta_kT(self):
        return self._cpp_impl.dkT*u.keV
    @delta_kT.setter
    @u.quantity_input
    def delta_kT(self, dkT: u.Quantity[u.keV]):
        self._cpp_impl.dkT = dkT.to_value(u.keV)

    @property
    def delta_z(self):
        return self._cpp_impl.dz
    @delta_z.setter
    @u.quantity_input
    def delta_z(self, dz: float):
        self._cpp_impl.dz = dz

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
    model: `str`
        choice of model for accretion column. One of: mcvspec, saxton, cropper, or wu
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
    abs_err : float
        absolute tolerance for integration and root finding
    rel_err :
        relative tolerance for integration
    delta_kT : `~astropy.units.Quantity`
        maximum allowed temperature grid spacing
    delta_z: float
        maximum allowed fractional altitude grid spacing

    Attributes
    __________
    white_dwarf : `~pymcvspec.white_dwarf`
        container for white dwarf properties
    accretion_column : `~pymcvspec.accretion_column`
        container for accretion column properties
    tolerance : `~pymcvspec.tolerance`
        container for numeric tolerances
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
        model: str,
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
        abs_err: float = 1e-8,
        rel_err: float = 1e-6,
        delta_kT: u.Quantity[u.keV] = 0.5*u.keV,
        delta_z: float = 0.1,
        solve: bool = True
    ) -> None:
        self.model = model
        self.white_dwarf = white_dwarf(
            mass,
            b_field,
            orbital_inclination,
            magnetospheric_radius,
            corotation_radius,
            distance
        )
        self.accretion_column = accretion_column(
            mdot,
            accretion_area,
            metallicity,
            shock_ratio,
            column_magnetic_colatitude
        )
        self.tolerance = tolerance(
            abs_err,
            rel_err,
            delta_kT,
            delta_z
        )

        if self.model == "mcvspec":
            self._cpp_impl = _cataclysmic_variable(
                self.white_dwarf._cpp_impl,
                self.accretion_column._cpp_impl,
                self.tolerance._cpp_impl
            )
            self.geometry = dipole(self._cpp_impl.column_coord)
        elif self.model == "saxton":
            self._cpp_impl = _saxton_cv(
                self.white_dwarf._cpp_impl,
                self.accretion_column._cpp_impl,
                self.tolerance._cpp_impl
            )
        elif self.model == "cropper":
            self._cpp_impl = _cropper_cv(
                self.white_dwarf._cpp_impl,
                self.accretion_column._cpp_impl,
                self.tolerance._cpp_impl
            )
        elif self.model == "wu":
            self._cpp_impl = _wu_cv(
                self.white_dwarf._cpp_impl,
                self.accretion_column._cpp_impl,
                self.tolerance._cpp_impl
            )
        else:
            raise RuntimeError("model must be one of: mcvspec, saxton, cropper, or wu")

        if solve:
            status = self._cpp_impl.solve()
            if status == -1:
                raise RuntimeError(("No valid solution found: Column does not "
                                    "reach WD surface for any shock height")
                                )

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
        return self._cpp_impl.density_to_ne
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
    def density_unit(self):
        return self._cpp_impl.density_converter*u.g/u.cm**3
    @property
    def position(self):
        return self._cpp_impl.position*u.dimensionless_unscaled
    @property
    def altitude(self):
        return self._cpp_impl.altitude*u.cm
    @property
    def volume_element(self):
        return self._cpp_impl.volume_element*u.cm**3
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

        apec_unit = ((u.cm**3)/u.s)
        dE = np.diff(energy_bins)
        continuum = pyatomdb.spectrum.CIESession()
        continuum.dolines = False
        continuum.dopseudo = False
        continuum.set_abund(
            atomic_charges[2:],
            self.accretion_column.metallicity
        )
        continuum.set_response(energy_bins.to_value(u.keV), raw=True)

        lines = pyatomdb.spectrum.CIESession()
        lines.docont = False
        lines.set_abund(atomic_charges[2:], self.accretion_column.metallicity)
        lines.set_response(energy_bins.to_value(u.keV), raw=True)

        flux_integrand = np.zeros((len(self.position), len(energy_bins)-1))
        flux_integrand /= (u.s*u.keV*u.cm**2)
        flux = np.zeros(len(energy_bins)-1)/(u.s*u.keV*u.cm**2)

        profile = zip(self.electron_temperature,
                      self.ion_temperature,
                      self.electron_density,
                      self.ion_density,
                      self.volume_element)

        for i, (kT_e,kT_i,n_e,n_i,dV) in enumerate(profile):
            emissivity = continuum.return_spectrum(kT_e.to_value(u.keV))
            emissivity += lines.return_spectrum(kT_i.to_value(u.keV))
            norm = n_e*n_i*dV/(4*np.pi*self.white_dwarf.distance**2)
            flux_integrand[i] = norm*emissivity*apec_unit/dE

        for i, dw in enumerate(np.diff(self.position)):
            flux += 0.5*dw*(flux_integrand[i]+flux_integrand[i+1])

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
