#include "mcvspec.hh"

extern "C"
void Polarspec(const RealArray& energy, const RealArray& params, int spectrum_num, RealArray& flux, RealArray& err, const string& init_string)
{
    flux.resize(energy.size()-1,0);
    err.resize(energy.size()-1,0);

    reflection_sel = params[0]; // how to apply reflection 0 = off, 1 = on, apply at each layer, 2 = on, apply once to whole spectrum
    b_field = params[1]*1e6; // surface b field [gauss]
    fractional_area = params[2]; //fractional accretion area
    luminosity = params[3]*1e33; // luminosity [ergs/s]
    mass = params[4]*solar_mass; // WD mass [grams]
    col_abund = params[5]; // accretion column abundance [solar abundances]
    wd_abund = params[6]; // wd surface abundance [solar abundances]
    cos_incl = params[7]; // cos inclination angle

    wd_radius = Calculate_White_Dwarf_Radius(mass);
    accretion_area = fractional_area*4*pi*wd_radius*wd_radius;
    accretion_rate = Calculate_Accretion_Rate(mass, luminosity, wd_radius);
    specific_accretion = accretion_rate/accretion_area;
    free_fall_velocity = sqrt(2.*grav_const*mass*((1./wd_radius)-(1./mag_radius)));
    shock_height = Calculate_B_Free_Shock_Height(free_fall_velocity, accretion_rate);
    shock_electron_dens = Calculate_Electron_Density(specific_accretion, b_free_shock_height);
    shock_temperature = Calculate_Shock_Temperature(wd_radius, shock_height, mag_radius);
    epsilon_zero = Calculate_Epsilon(b_field, shock_temperature, shock_electron_dens, shock_height);
    epsilon_shock = Root_Finder(Epsilon_Diff, Epsilon_Diff_Derivative, &epsilon_zero, 1e5, 100000, 1e-6);;

    Shock_Height_Shooting(0.005, 100);
    MCVspec_Spectrum(energy, spectrum_num, flux, init_string);
}
