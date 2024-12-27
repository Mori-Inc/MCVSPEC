#ifndef MCVSPEC_VARS_H
#define MCVSPEC_VARS_H

#include <cmath>
#include <valarray>

using std::pow;
using std::valarray;

// physical constants
// these never change within the model, they can only be modified by the code itself
// these are mostly physical constants of the universe
namespace  mcvspec_physical_constants
{
    const double pi = 3.14159265358979323846264338327950; // pi
    const double solar_mass = 1.989100e+33; // solar mass in grams
    const double solar_radius = 6.9599e10; // solar radius in cm
    const double grav_const = 6.672590e-8; // Newton's gravitational constant in cgs units
    const double wd_mol_mass = 2.0; // mean molecular mass of White dwarf
    const double col_mol_mass = 0.615; // mean molecular mass of acretion column (see Cropper 1999)
    const double bremss_const = 6.99e16; // bremsstrahlung constant (see Saxton: this coefficient gives the flux norm close to Suleimanov's model for the same WD mass input - 18% difference)
    const double boltz_const = 1.380658e-16; // Boltzmann constant in cgs
    const double boltz_const_kev = 8.617333262e-8;
    const double hydrg_mass = 1.672623e-24; // Mass of hydrogen in grams
    const double alpha = 2.0; // Thermodynamic Constant
    const double beta = 3.85; // Thermodynamic Constant
    const double proton_molar_mass = 1.0072764665789; // molar mass of a proton
    const double electron_molar_mass = 5.485799090441e-4; // molar mass of an electron
    const double mass_limit = 5.816*solar_mass/(pow(wd_mol_mass,2.0));
    const double helium_ratio = (2*col_mol_mass-electron_molar_mass-proton_molar_mass)/(2*electron_molar_mass+4*proton_molar_mass-6*col_mol_mass);
    const double electron_ion_ratio = (1+2*helium_ratio)/(1+helium_ratio);
    const double apec_norm = 2.62511E-34; // (1/4)*(1 km/1 kpc)^2 normalization needed for apec
}

// constants of the model
// these never change within the model, they can only be modified by the code itself
// these are semi-arbitrary decisions about how the model should operate
namespace mcvspec_model_constants
{
    const double shock_velocity = 0.25; // downstream velocity at shock, normalized to upstream free fall velocity at shock
    const double initial_veloicty = 0; // velocity at WD surface, normalized to upstream velocity at shock
    const double initial_height = 0;
    const double absolute_err = 1e-3;
    const double relative_err = 1e-2;
}

// model values
namespace mcvspec
{
    extern double mass, b_field, p_spin, luminosity, col_abund, wd_abund, mag_ratio, cos_incl;
    extern double accretion_area, accretion_rate, shock_height, velocity_at_shock, shock_temperature;
    extern double free_fall_velocity, wd_radius, mag_radius, shock_electron_dens, b_free_shock_height;
    extern double epsilon_s, epsilon_zero; // ratio of brems cooling time to cytclotron cooling time
    extern int reflection_sel;

    extern const int max_num_grid_points;
    extern int num_grid_points;
    extern valarray<double> density;
    extern valarray<double> pressure;
    extern valarray<double> temperature;
    extern valarray<double> altitude;
    extern valarray<double> velocity;
    extern valarray<double> electron_dens;
    extern valarray<double> offset_altitude;
}

#endif
