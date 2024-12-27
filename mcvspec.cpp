#include "mcvspec_core.hh"
#include <limits>
using std::cbrt;
using std::numeric_limits;

extern "C"
void IPspec(const RealArray& energy, const RealArray& params, int spectrum_num, RealArray& flux, RealArray& err, const string& init_string)
{
    flux.resize(energy.size()-1,0);
    err.resize(energy.size()-1,0);

    reflection_sel = params[0]; // how to apply reflection 0 = off, 1 = on, apply at each layer, 2 = on, apply once to whole spectrum
    p_spin = params[1]; // spin period [s]; assumes spin equilibrium
    accretion_rate = params[2]; //fractional accretion area
    luminosity = params[3]*1e33; // luminosity [ergs/s]
    mass = params[4]*solar_mass; // WD mass [grams]
    col_abund = params[5]; // accretion column abundance [solar abundances]
    wd_abund = params[6]; // wd surface abundance [solar abundances]
    cos_incl = params[7]; // cos inclination angle

    wd_radius = Calculate_White_Dwarf_Radius();
    mag_radius = cbrt(4*pi*pi*grav_const*mass*p_spin*p_spin);
    accretion_rate = Calculate_Accretion_Rate();
    b_field = Calculate_Magnetic_Field();

    if(wd_radius < mag_radius){
        free_fall_velocity = sqrt(2.*grav_const*mass*((1./wd_radius)-(1./mag_radius)));
    }
    else{
        free_fall_velocity = sqrt(2.*grav_const*mass/wd_radius);
    }

    shock_height = Estimate_Shock_Height();
    epsilon_s = Root_Finder(1e5, 100000, 1e-6);
    Shock_Height_Shooting(0.005, 100);
    MCVspec_Spectrum(energy, spectrum_num, flux, init_string);
}

extern "C"
void Deqspec(const RealArray& energy, const RealArray& params, int spectrum_num, RealArray& flux, RealArray& err, const string& init_string)
{
    flux.resize(energy.size()-1,0);
    err.resize(energy.size()-1,0);

    reflection_sel = params[0]; // how to apply reflection 0 = off, 1 = on, apply at each layer, 2 = on, apply once to whole spectrum
    mag_ratio = params[1]; // magnetospheric to wd radius ratio
    accretion_area = params[2]; //fractional accretion area
    luminosity = params[3]*1e33; // luminosity [ergs/s]
    mass = params[4]*solar_mass; // WD mass [grams]
    col_abund = params[5]; // accretion column abundance [solar abundances]
    wd_abund = params[6]; // wd surface abundance [solar abundances]
    cos_incl = params[7]; // cos inclination angle

    wd_radius = Calculate_White_Dwarf_Radius();
    mag_radius = wd_radius*mag_ratio;
    accretion_rate = Calculate_Accretion_Rate();
    b_field = Calculate_Magnetic_Field();

    if(wd_radius < mag_radius){
        free_fall_velocity = sqrt(2.*grav_const*mass*((1./wd_radius)-(1./mag_radius)));
    }
    else{
        free_fall_velocity = sqrt(2.*grav_const*mass/wd_radius);
    }

    shock_height = Estimate_Shock_Height();
    epsilon_s = Root_Finder(1e5, 100000, 1e-6);
    Shock_Height_Shooting(0.005, 100);
    MCVspec_Spectrum(energy, spectrum_num, flux, init_string);
}

extern "C"
void Polarspec(const RealArray& energy, const RealArray& params, int spectrum_num, RealArray& flux, RealArray& err, const string& init_string)
{
    flux.resize(energy.size()-1,0);
    err.resize(energy.size()-1,0);

    reflection_sel = params[0]; // how to apply reflection 0 = off, 1 = on, apply at each layer, 2 = on, apply once to whole spectrum
    b_field = params[1]*1e6; // surface b field [gauss]
    accretion_area = params[2]; //fractional accretion area
    luminosity = params[3]*1e33; // luminosity [ergs/s]
    mass = params[4]*solar_mass; // WD mass [grams]
    col_abund = params[5]; // accretion column abundance [solar abundances]
    wd_abund = params[6]; // wd surface abundance [solar abundances]
    cos_incl = params[7]; // cos inclination angle

    wd_radius = Calculate_White_Dwarf_Radius();
    mag_radius = numeric_limits<double>::max();
    accretion_rate = Calculate_Accretion_Rate();
    b_field = Calculate_Magnetic_Field();

    if(wd_radius < mag_radius){
        free_fall_velocity = sqrt(2.*grav_const*mass*((1./wd_radius)-(1./mag_radius)));
    }
    else{
        free_fall_velocity = sqrt(2.*grav_const*mass/wd_radius);
    }

    shock_height = Estimate_Shock_Height();
    epsilon_s = Root_Finder(1e5, 100000, 1e-6);
    Shock_Height_Shooting(0.005, 100);
    MCVspec_Spectrum(energy, spectrum_num, flux, init_string);
}
